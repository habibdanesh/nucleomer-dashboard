"""Stats: per-model distributions (significance per k, p-values, values per k, edges) and a cross-model table."""
from functools import lru_cache

import numpy as np
import plotly.graph_objects as go
from dash import Input, Output

from .. import store as S
from ..kmers import OFF
from ..theme import TAB10, empty
from ..ui import box, fig, page, row, table


def layout():
    return page([],
                row(box("Significant k-mers", fig("st-sig", 300),
                        tip="Fraction of all k-mers of each length with node p below the threshold (top bar)."),
                    box("Node p-values", fig("st-p", 300),
                        tip="Distribution of node p-values for each k >= 4 (darker = longer k-mers)."),
                    box("Value per k", fig("st-v", 300),
                        tip="All k-mers of each length. Box = 25-75%, line = median, whiskers = 0.1-99.9%.")),
                row(box("Edge p-values", fig("st-e", 300), tip="Edge p-values of insertion and deletion edges in the graph."),
                    box("Background predictions", fig("st-bg", 300),
                        tip="Model prediction for each background sequence without any insertion. Dashed = median.")),
                box("All models", table("st-tab", page_size=20), tip="One row per model in the bundle."))


@lru_cache(64)
def per_k(m, kind, pthr):
    v, p = S.value(m, kind), np.asarray(S.arr(m, "pval"))
    out = []
    for k in range(1, S.meta(m)["maxk"] + 1):
        sl = S.OFF_K(k)
        vk, pk = v[sl], p[sl]
        q = np.nanpercentile(vk, [0.1, 25, 50, 75, 99.9])
        out.append(dict(k=k, n=len(pk), nsig=int((pk < pthr).sum()), q=q,
                        phist=np.histogram(pk, bins=np.linspace(0, 1, 51))[0]))
    return out


@lru_cache(16)
def edge_stats(m, pthr):
    E = S.edges(m)
    p, d = E["pval"], E["is_del"]
    bins = np.linspace(0, 1, 51)
    return dict(ins=np.histogram(p[~d], bins)[0], dele=np.histogram(p[d], bins)[0], n=len(p),
                nsig=int((p < pthr).sum()), nsig_ins=int((p[~d] < pthr).sum()), nsig_del=int((p[d] < pthr).sum()))


def register(app):
    @app.callback(Output("st-sig", "figure"), Output("st-p", "figure"), Output("st-v", "figure"),
                  Output("st-e", "figure"), Output("st-bg", "figure"),
                  Input("model", "value"), Input("value", "value"), Input("pthr", "value"))
    def update(m, kind, pthr):
        if not m:
            return (empty(),) * 5
        pthr = pthr or S.PVAL
        st = per_k(m, kind, pthr)
        ks = [s["k"] for s in st]
        f1 = go.Figure(go.Bar(x=ks, y=[s["nsig"] / s["n"] for s in st], marker_color=TAB10[0],
                              customdata=[[s["nsig"], s["n"]] for s in st],
                              hovertemplate="k %{x}: %{customdata[0]:,} / %{customdata[1]:,}<extra></extra>"))
        f1.update_layout(xaxis=dict(title="k", dtick=1), yaxis=dict(title=f"fraction with p < {pthr:g}"))
        f2 = go.Figure()
        x = np.linspace(0, 1, 51)
        for s in st:
            if s["k"] < 4:
                continue
            f2.add_scatter(x=x[:-1], y=s["phist"] / s["n"] * 50, mode="lines", line_shape="hv", name=f"{s['k']}",
                           line=dict(width=1, color=f"rgba(31,119,180,{0.25 + 0.75 * s['k'] / 11:.2f})"))
        f2.update_layout(xaxis_title="p", yaxis=dict(title="density", type="log"),
                         legend=dict(title_text="k", x=1.02, xanchor="left", y=1, font_size=10, tracegroupgap=0))
        f3 = go.Figure()
        for s in st:
            q = s["q"]
            f3.add_box(x=[s["k"]], q1=[q[1]], median=[q[2]], q3=[q[3]], lowerfence=[q[0]], upperfence=[q[4]],
                       marker_color=TAB10[0], line_width=1, showlegend=False, hoverinfo="skip")
        f3.update_layout(xaxis=dict(title="k", dtick=1), yaxis_title=S.VALUE_LABEL[kind] + " (0.1-99.9%)")
        es = edge_stats(m, pthr)
        f4 = go.Figure()
        for name, h, c in [(f"insertion ({es['nsig_ins']:,} sig.)", es["ins"], TAB10[0]),
                           (f"deletion ({es['nsig_del']:,} sig.)", es["dele"], TAB10[1])]:
            f4.add_scatter(x=x[:-1], y=h, mode="lines", line_shape="hv", name=name, line=dict(color=c, width=1.2))
        f4.update_layout(xaxis_title="edge p", yaxis=dict(title="edges", type="log"),
                         legend=dict(x=0.98, xanchor="right", y=0.98))
        pb = np.asarray(S.meta(m)["pred_before"])
        f5 = go.Figure(go.Histogram(x=pb, nbinsx=30, marker_color="#999", name="backgrounds"))
        f5.add_vline(x=np.median(pb), line_dash="dash", line_width=1)
        f5.update_layout(title=f"{len(pb)} backgrounds")
        f5.update_layout(xaxis_title="prediction (no insertion)", yaxis_title="backgrounds", showlegend=False)
        return f1, f2, f3, f4, f5

    @app.callback(Output("st-tab", "data"), Output("st-tab", "columns"), Input("pthr", "value"))
    def tab(pthr):
        pthr = pthr or S.PVAL
        rows = []
        for m in S.models():
            md = S.meta(m)
            p = np.asarray(S.arr(m, "pval"))
            top = sorted(md["top_cores"].items(), key=lambda t: -t[1][1])[:3]
            rows.append(dict(model=m, family=md["family"], assay=md["assay"], maxk=md["maxk"], backgrounds=md["n_backgrounds"],
                             sig=int((p < pthr).sum()), nodes=md["n_graph_nodes"], edges=md["n_edges"],
                             core=", ".join(dict.fromkeys(v[0] for _, v in top))))
        cols = [{"name": n, "id": i} for n, i in [("model", "model"), ("family", "family"), ("assay", "assay"), ("max k", "maxk"),
                                                  ("backgrounds", "backgrounds"), (f"k-mers p < {pthr:g}", "sig"),
                                                  ("graph nodes", "nodes"), ("graph edges", "edges"), ("top cores", "core")]]
        return rows, cols

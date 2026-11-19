"""Compare: two models k-mer by k-mer, all models against each other, and top k-mers x models."""
from functools import lru_cache

import numpy as np
import plotly.graph_objects as go
from dash import Input, Output, State, ctx, no_update
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.stats import rankdata

from .. import store as S
from ..kmers import OFF, revcomp_gid
from ..theme import TAB10, empty, key
from ..ui import box, check, ctl, dd, fig, num, page, row, slider, table


def layout():
    ms = S.models()
    side = [
        ctl("Model A", dd("cp-a", ms, ms[0])),
        ctl("Model B", dd("cp-b", ms, "ATAC" if "ATAC" in ms else ms[1])),
        ctl("k", slider("cp-k", 4, 11, 8), "all comparisons use every k-mer of this length"),
        check("cp-rc", "merge reverse complements", True, "a k-mer and its reverse complement count once (the stronger strand)"),
        ctl("Top per model (heatmap)", num("cp-top", 10, min=1, max=100)),
        ctl("Top of A (table)", num("cp-ntab", 200, min=10, max=5000, step=10)),
    ]
    return page(side,
                row(box("A vs B", fig("cp-ab", 420, ),
                        tip="Every k-mer of length k: value in A (x) vs B (y). Grey = density of all k-mers; colours = "
                            "significant in A only, B only or both. Spearman correlation over all k-mers in the title."),
                    box("All models (Spearman)", fig("cp-corr", 420),
                        tip="Rank correlation of values over all k-mers of length k, for every pair of models; "
                            "ordered by hierarchical clustering.")),
                row(box("Top k-mers of A, looked up in B", table("cp-tab", page_size=14),
                        tip="A's strongest significant k-mers with their value, p and percentile in B's complete results "
                            "(all k-mers, not only B's pruned graph)."),
                    box("Top k-mers x models", fig("cp-heat", 560),
                        tip="The top significant k-mers of every model (union), coloured by value in each model; "
                            "rows clustered.")),
                about="Two models k-mer by k-mer, and all models side by side.")


def vals(m, kind, k, merge):
    """Value of every k-mer of length k (canonical ids if merged: max over the two strands)."""
    g = np.arange(OFF[k], OFF[k + 1], dtype=np.int64)
    v = S.value(m, kind)
    p = np.asarray(S.arr(m, "pval"))
    if not merge:
        return g, v[g], p[g]
    r = revcomp_gid(g)
    keep = g <= r
    return g[keep], np.maximum(v[g], v[r])[keep], np.minimum(p[g], p[r])[keep]


@lru_cache(32)
def ranks(m, kind, k, merge):
    _, v, _ = vals(m, kind, k, merge)
    return rankdata(v, nan_policy="omit")


@lru_cache(16)
def corr(models, kind, k, merge):
    R = np.stack([ranks(m, kind, k, merge) for m in models])
    return np.corrcoef(R)


def register(app):
    @app.callback(Output("cp-ab", "figure"), Output("cp-tab", "data"), Output("cp-tab", "columns"),
                  Input("cp-a", "value"), Input("cp-b", "value"), Input("cp-k", "value"), Input("cp-rc", "value"),
                  Input("cp-ntab", "value"), Input("value", "value"), Input("pthr", "value"))
    def ab(a, b, k, merge, ntab, kind, pthr):
        if not a or not b:
            return empty(), [], []
        if k > min(S.meta(a)["maxk"], S.meta(b)["maxk"]):
            return empty(f"k > maxk of {a} or {b}"), [], []
        pthr = pthr or S.PVAL
        g, va, pa = vals(a, kind, k, bool(merge))
        _, vb, pb = vals(b, kind, k, bool(merge))
        lab = S.VALUE_LABEL[kind]
        sa, sb = pa < pthr, pb < pthr
        f = go.Figure()
        if len(g) > 20000:                      # density for the bulk, points only for significant k-mers
            H, xe, ye = np.histogram2d(va, vb, bins=150)
            f.add_heatmap(z=np.log10(1 + H.T), x=(xe[1:] + xe[:-1]) / 2, y=(ye[1:] + ye[:-1]) / 2, colorscale="Greys",
                          showscale=False, hoverinfo="skip")
            bulk = np.zeros(len(g), bool)
        else:
            bulk = ~(sa | sb)
        seqs_all = None
        if len(g) > 20000:
            key(f, "all k-mers (density)", "#888")
        groups = [("neither", bulk, "#c8c8c8"), (f"{a} only", sa & ~sb, TAB10[0]), (f"{b} only", sb & ~sa, TAB10[1]),
                  ("both", sa & sb, TAB10[3])]
        for name, msk, col in groups:
            idx = np.where(msk)[0]
            if len(idx) > 20000:
                idx = idx[np.argsort(-np.maximum(va[idx], vb[idx]))[:20000]]
            if not len(idx):
                continue
            f.add_scattergl(x=va[idx], y=vb[idx], mode="markers", name=f"{name} ({msk.sum():,})",
                            marker=dict(size=3 if name == "neither" else 4, color=col, opacity=0.8),
                            customdata=S.decode_np(g[idx]), hovertemplate="%{customdata}<br>%{x:.3f}, %{y:.3f}<extra></extra>")
        rho = np.corrcoef(rankdata(va), rankdata(vb))[0, 1]
        f.update_layout(xaxis_title=f"{lab} ({a})", yaxis_title=f"{lab} ({b})",
                        title=f"k = {k}, {len(g):,} k-mers, Spearman {rho:.2f}",
                        legend=dict(x=0.01, y=0.99, bgcolor="rgba(255,255,255,0.7)"))
        # table: top significant k-mers of A, with their values in B's complete (unpruned) results
        o = np.where(sa)[0]
        o = o[np.argsort(-va[o])][:ntab or 200]
        pct = 100 * rankdata(vb)[o] / len(vb)
        data = [dict(seq=s, a=round(float(x), 4), pa=S.fmt_p(y), b=round(float(z), 4), pb=S.fmt_p(w), pct=round(float(q), 1))
                for s, x, y, z, w, q in zip(S.decode_np(g[o]), va[o], pa[o], vb[o], pb[o], pct)]
        cols = [{"name": n, "id": i} for n, i in [("k-mer", "seq"), (f"{a}", "a"), (f"p {a}", "pa"), (f"{b}", "b"),
                                                  (f"p {b}", "pb"), (f"percentile in {b}", "pct")]]
        return f, data, cols

    @app.callback(Output("cp-corr", "figure"), Output("cp-heat", "figure"),
                  Input("cp-k", "value"), Input("cp-rc", "value"), Input("cp-top", "value"), Input("value", "value"),
                  Input("pthr", "value"))
    def allm(k, merge, top, kind, pthr):
        ms = tuple(m for m in S.models() if S.meta(m)["maxk"] >= k)
        if len(ms) < 2:
            return empty(), empty()
        C = corr(ms, kind, k, bool(merge))
        o = leaves_list(linkage(1 - C[np.triu_indices(len(ms), 1)], "average"))
        mo = [ms[i] for i in o]
        f1 = go.Figure(go.Heatmap(z=C[np.ix_(o, o)], x=mo, y=mo, colorscale="RdBu_r", zmin=-1, zmax=1, xgap=1, ygap=1,
                                  colorbar=dict(thickness=10, title=dict(text="Spearman", side="right")), hovertemplate="%{x} / %{y}: %{z:.2f}<extra></extra>"))
        f1.update_layout(yaxis_autorange="reversed", xaxis_tickangle=-45, title=f"k = {k}")
        # top k-mers per model (significant), union, rows clustered
        pthr = pthr or S.PVAL
        rows = set()
        V = {}
        for m in mo:
            g, v, p = vals(m, "score", k, bool(merge))
            V[m] = v
            i = np.where(p < pthr)[0]
            rows.update(i[np.argsort(-v[i])][:top or 10].tolist())
        rows = np.array(sorted(rows))
        g, _, _ = vals(mo[0], kind, k, bool(merge))
        M = np.stack([vals(m, kind, k, bool(merge))[1][rows] for m in mo], 1)
        Mz = np.nan_to_num(M)
        ro = leaves_list(linkage(Mz, "average", metric="cosine")) if len(rows) > 2 else np.arange(len(rows))
        seqs = S.decode_np(g[rows])[ro]
        f2 = go.Figure(go.Heatmap(z=M[ro], x=mo, y=seqs, colorscale="Viridis", xgap=0.5, ygap=0.5,
                                  colorbar=dict(title=dict(text=S.VALUE_LABEL[kind], side="right"), thickness=10),
                                  hovertemplate="%{y}<br>%{x}: %{z:.3f}<extra></extra>"))
        f2.update_layout(yaxis=dict(autorange="reversed", tickfont=dict(family="Menlo, monospace", size=8),
                                    showticklabels=len(rows) <= 120), xaxis_tickangle=-45, title=f"{len(rows)} k-mers")
        return f1, f2

    @app.callback(Output("kmer", "value", allow_duplicate=True),
                  Input("cp-ab", "clickData"), Input("cp-heat", "clickData"), Input("cp-tab", "active_cell"),
                  State("cp-tab", "derived_viewport_data"), prevent_initial_call=True)
    def pick(c1, c2, cell, rows):
        t = ctx.triggered_id
        if t == "cp-tab":
            return rows[cell["row"]]["seq"] if cell and rows else no_update
        c = c1 if t == "cp-ab" else c2
        if not c or not c["points"]:
            return no_update
        pt = c["points"][0]
        return pt.get("customdata") or (pt.get("y") if t == "cp-heat" else no_update)

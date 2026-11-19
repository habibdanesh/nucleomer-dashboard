"""K-mers: every k-mer of one model; density of all, points + table for the ones passing the filters."""
import numpy as np
import plotly.graph_objects as go
from dash import Input, Output, State, no_update

from .. import store as S
from ..kmers import OFF, revcomp_gid
from ..theme import TAB10, empty, key, legend_top
from ..ui import box, check, ctl, dd, fig, page, rng, row, table, txt

SEL = "selected (table rows)"
TIPS = {"seq": "k-mer (with 'merge reverse complements': the stronger strand)", "rc": "reverse complement",
        "value": "value chosen in the top bar (score, lift or prediction)",
        "pred": "mean model prediction over the backgrounds with the k-mer inserted",
        "p": "node p-value: bootstrap over k-mers with the same dinucleotide content (1000 samples, so p < 0.001 is the floor)",
        "core": "shortest k-mer reached by deletions that keeps >= 80% of this k-mer's prediction (Nucleomer core)"}


def layout():
    side = [
        ctl("k", rng("ex-k", 1, 11, [6, 11]), "k-mer lengths to include"),
        ctl("Pattern", txt("ex-pat", placeholder="regex, either strand", ),
            "regular expression matched on either strand, e.g. GATA or GA.TAA"),
        check("ex-sig", "significant only", True, "only k-mers with node p below the p threshold in the top bar"),
        check("ex-rc", "merge reverse complements", True, "one row per k-mer / reverse-complement pair: the stronger strand"),
        ctl("Rows", dd("ex-n", [100, 500, 2000, 10000], 500), "how many k-mers (strongest first) go into the table and the blue points"),
    ]
    return page(side,
                row(box("Value by k-mer length", fig("ex-land", 330),
                        tip="Grey: density of all k-mers of each length (log scale, each row scaled to its own maximum). "
                            "Blue: the k-mers in the table. Click a point to select it."),
                    box("Value vs significance", fig("ex-volc", 330),
                        tip="Grey: density of all k-mers. Blue: the k-mers in the table. Dashed: p threshold. "
                            "-log10 p stops at 3 because p comes from 1000 bootstrap samples.")),
                box("K-mers passing the filters", table("ex-table", page_size=14, tips=TIPS),
                    tip="Sorted by value. Click a row to select the k-mer in every tab."),
                about="Every k-mer of the model chosen in the top bar. Filter on the left; click a point or row to select a k-mer.")


def selected(m, kind, krange, sig_only, merge, pattern, pthr, n):
    g = np.arange(OFF[krange[0]], OFF[min(krange[1], S.meta(m)["maxk"]) + 1], dtype=np.int64)
    v = S.value(m, kind)
    p = np.asarray(S.arr(m, "pval"))
    keep = g
    if sig_only:
        keep = keep[p[keep] < pthr]
    if merge:                                                # one row per k-mer/RC pair: the stronger strand
        keep = keep[np.argsort(-v[keep], kind="stable")]
        keep = keep[np.unique(S.canonical(keep), return_index=True)[1]]
    keep = S.search(keep, pattern)
    keep = keep[np.argsort(-v[keep])][:n]
    return g, keep


def register(app):
    @app.callback(Output("ex-land", "figure"), Output("ex-volc", "figure"), Output("ex-table", "data"),
                  Output("ex-table", "columns"),
                  Input("model", "value"), Input("value", "value"), Input("pthr", "value"), Input("ex-k", "value"),
                  Input("ex-sig", "value"), Input("ex-rc", "value"), Input("ex-pat", "value"), Input("ex-n", "value"))
    def update(m, kind, pthr, krange, sig, merge, pattern, n):
        if not m:
            return empty(), empty(), [], []
        try:
            g, keep = selected(m, kind, krange, bool(sig), bool(merge), pattern, pthr or S.PVAL, n)
        except Exception as e:  # bad regex
            return empty(str(e)), empty(), [], []
        v, nlp, k = S.value(m, kind), S.value(m, "nlp"), None
        lab = S.VALUE_LABEL[kind]
        from ..kmers import k_of
        kg = k_of(g)
        vg = v[g]
        ok = np.isfinite(vg)
        # density of all k-mers per k (each row scaled to its own max so short k are visible)
        lo, hi = np.nanpercentile(vg, [0.01, 100])
        bins = np.linspace(lo, hi, 161)
        ks = np.arange(krange[0], min(krange[1], S.meta(m)["maxk"]) + 1)
        H = np.array([np.histogram(vg[ok & (kg == kk)], bins)[0] for kk in ks], float)
        H = np.log10(1 + H) / np.log10(1 + H).max(1, keepdims=True).clip(1e-9)
        f1 = go.Figure(go.Heatmap(z=H, x=(bins[1:] + bins[:-1]) / 2, y=ks, colorscale="Greys", showscale=False,
                                  hoverinfo="skip", zmin=0, zmax=1.4))
        seqs = S.decode_np(keep)
        kk = k_of(keep)
        jitter = np.random.default_rng(0).uniform(-0.3, 0.3, len(keep))
        f1.add_scattergl(x=v[keep], y=kk + jitter, mode="markers", marker=dict(size=4, color=TAB10[0], opacity=0.7),
                         customdata=seqs, hovertemplate="%{customdata}<br>%{x:.3f}<extra></extra>", name=SEL)
        key(f1, "all k-mers (density)", "#888")
        f1.update_layout(xaxis_title=lab, yaxis_title="k", yaxis_dtick=1, clickmode="event")
        legend_top(f1)

        # value vs -log10 p (density + selection)
        H2, xe, ye = np.histogram2d(vg[ok], nlp[g][ok], bins=[bins, np.linspace(0, 3.05, 62)])
        f2 = go.Figure(go.Heatmap(z=np.log10(1 + H2.T), x=(xe[1:] + xe[:-1]) / 2, y=(ye[1:] + ye[:-1]) / 2,
                                  colorscale="Greys", showscale=False, hoverinfo="skip"))
        f2.add_scattergl(x=v[keep], y=nlp[keep] + jitter * 0.1, mode="markers",
                         marker=dict(size=4, color=TAB10[0], opacity=0.7), customdata=seqs,
                         hovertemplate="%{customdata}<br>%{x:.3f}<extra></extra>", name=SEL)
        f2.add_hline(y=-np.log10(pthr or S.PVAL), line_dash="dash", line_width=1, line_color="#888")
        key(f2, "all k-mers (density)", "#888")
        key(f2, f"p = {pthr or S.PVAL:g}", "#888", dash="dash")
        f2.update_layout(xaxis_title=lab, yaxis_title="-log10 p")
        legend_top(f2)

        pred, pv, core = S.arr(m, "pred"), S.arr(m, "pval"), S.arr(m, "core")
        cs = S.decode_np(np.maximum(np.asarray(core[keep]), 0))
        data = [dict(seq=s, rc=S.decode_np(revcomp_gid([g_]))[0], k=int(k_), value=round(float(v[g_]), 4),
                     pred=round(float(pred[g_]), 4), p=S.fmt_p(pv[g_]),
                     core=(c if core[g_] >= 0 else ""))
                for s, g_, k_, c in zip(seqs, keep, kk, cs)]
        cols = [{"name": {"value": lab, "seq": "k-mer"}.get(c, c), "id": c} for c in ["seq", "rc", "k", "value", "pred", "p", "core"]]
        return f1, f2, data, cols

    @app.callback(Output("kmer", "value", allow_duplicate=True),
                  Input("ex-land", "clickData"), Input("ex-volc", "clickData"), Input("ex-table", "active_cell"),
                  State("ex-table", "derived_viewport_data"), prevent_initial_call=True)
    def pick(c1, c2, cell, rows):
        from dash import ctx
        t = ctx.triggered_id
        if t == "ex-table" and cell and rows:
            return rows[cell["row"]]["seq"]
        c = c1 if t == "ex-land" else c2
        if c and c["points"] and "customdata" in c["points"][0]:
            return c["points"][0]["customdata"]
        return no_update

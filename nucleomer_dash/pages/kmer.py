"""K-mer: one k-mer across models, per background, under single-base edits, and its graph neighbours."""
import numpy as np
import plotly.graph_objects as go
from dash import Input, Output, State, ctx, html, no_update
from plotly.subplots import make_subplots

from .. import store as S
from ..kmers import KMAX, encode, k_of, rc, revcomp_gid
from ..theme import FAMILY_COLOR, TAB10, empty, key, legend_top
from ..ui import box, fig, hint, page, row, table

NB_TIPS = {"nbr": "k-mer one insertion (k+1) or one deletion (k-1) away", "dir": "out: edge from the selected k-mer; in: edge into it",
           "edge_p": "edge p-value: how often a random edit of the same source k-mer (bootstrap over its graph children) gains as much",
           "dvalue": "value of the neighbour minus value of the selected k-mer", "node_p": "node p-value of the neighbour"}


def layout():
    side = [html.Div(id="km-info", className="info")]
    return page(side,
                row(box("Across models", fig("km-models", 300),
                        tip="Value of this k-mer (dark) and its reverse complement (light) in every model. "
                            "Bar colour = model family; black outline = significant (p below the threshold)."),
                    box("Per background", fig("km-bg", 300),
                        tip="Prediction with the k-mer inserted minus prediction without it, one point per background sequence. "
                            "Grey = the median k-mer of the same length, for reference. Stored for all k <= 6 and the top "
                            "25,000 k-mers of each model (and their reverse complements).")),
                row(box("Single-base edits", fig("km-scan", 360),
                        tip="Change in value when one base is substituted, deleted, or inserted before a position. "
                            "Red = the edited k-mer scores higher, blue = lower."),
                    box("Graph neighbours", table("km-nbrs", page_size=12, tips=NB_TIPS),
                        tip="K-mers connected to this one in the model's Nucleomer graph. Click a row to move to it.")),
                about="One k-mer (type it in the top bar, or click one in any tab).")


def seq_edits(s):
    """All substitutions, deletions and insertions of s: list of (kind, pos, base, new_seq)."""
    out = []
    for i, c in enumerate(s):
        out += [("sub", i, b, s[:i] + b + s[i + 1:]) for b in "ACGT"]
        if len(s) > 1:
            out.append(("del", i, "", s[:i] + s[i + 1:]))
    if len(s) < KMAX:
        out += [("ins", i, b, s[:i] + b + s[i:]) for i in range(len(s) + 1) for b in "ACGT"]
    return out


def register(app):
    @app.callback(Output("km-info", "children"), Output("km-models", "figure"), Output("km-bg", "figure"),
                  Output("km-scan", "figure"), Output("km-nbrs", "data"), Output("km-nbrs", "columns"),
                  Input("kmer", "value"), Input("model", "value"), Input("value", "value"), Input("pthr", "value"))
    def update(s, m, kind, pthr):
        g = S.parse_kmer(s)
        if g is None or not m:
            msg = "enter a k-mer above" if not s else "not a valid k-mer (ACGT, length 1-11)"
            return html.Div(msg, className="muted"), empty(), empty(), empty(), [], []
        s = s.strip().upper()
        k, pthr = len(s), pthr or S.PVAL
        r = int(revcomp_gid([g])[0])
        lab = S.VALUE_LABEL[kind]
        ms = S.models()
        if k > S.meta(m)["maxk"]:
            return html.Div(f"{m} was run up to k = {S.meta(m)['maxk']}", className="muted"), empty(), empty(), empty(), [], []

        # info panel
        core = int(S.arr(m, "core")[g])
        cseq = S.decode_np([core])[0] if core >= 0 else None
        pos = s.find(cseq) if cseq else -1
        pos_rc = rc(s).find(cseq) if cseq and pos < 0 else -1
        v = S.value(m, kind)
        info = [html.Div(s, className="seq big"), html.Div(["RC ", rc(s)], className="seq muted"),
                html.Table([html.Tr([html.Td(a), html.Td(b)]) for a, b in [
                    ("k", k), (lab, f"{v[g]:.3f}"), ("RC " + lab, f"{v[r]:.3f}"),
                    ("p", S.fmt_p(S.arr(m, "pval")[g])), ("RC p", S.fmt_p(S.arr(m, "pval")[r])),
                    (html.Span(["core", hint("shortest k-mer reached by deletions that keeps >= 80% of this k-mer's prediction (Nucleomer core)")]),
                     cseq or "-")]], className="kv")]
        if cseq:
            al = (" " * pos + cseq) if pos >= 0 else (" " * pos_rc + cseq + "  (RC)")
            info.append(html.Pre(f"{s}\n{al}", className="align"))
        tt = S.meta(m)["tomtom"].get(s) or S.meta(m)["tomtom"].get(rc(s))
        if tt:
            info.append(html.Div([html.Div(["TOMTOM", hint("best JASPAR matches of this core (TOMTOM)")], className="box-title")] +
                                 [html.Div(f"{t['target']}  {t['p']:.1e}", className="small") for t in tt]))

        # across models
        f1 = go.Figure()
        for strand, gg, op in [("k-mer", g, 1.0), ("RC", r, 0.45)]:
            vals, ps, cols = [], [], []
            for mm in ms:
                ok = k <= S.meta(mm)["maxk"]
                vals.append(float(S.value(mm, kind)[gg]) if ok else None)
                ps.append(float(S.arr(mm, "pval")[gg]) if ok else None)
                cols.append(FAMILY_COLOR.get(S.meta(mm)["family"], TAB10[7]))
            f1.add_bar(x=ms, y=vals, name=strand, showlegend=False, marker=dict(color=cols, opacity=op,
                       line=dict(width=[2 if (p is not None and p < pthr) else 0 for p in ps], color="#222")),
                       customdata=ps, hovertemplate="%{x}<br>" + lab + " %{y:.3f}<br>p %{customdata:.3g}<extra>" + strand + "</extra>")
        for fam in dict.fromkeys(S.meta(mm)["family"] for mm in ms):
            key(f1, fam, FAMILY_COLOR.get(fam, TAB10[7]))
        key(f1, "significant", "#222", symbol="square-open")
        f1.update_layout(barmode="group", yaxis_title=lab, xaxis_tickangle=-45, margin=dict(b=70))
        legend_top(f1)

        # per-background effect: k-mer, RC, and the median k-mer of the same length (reference)
        med = S.median_kmer(m, k)
        f2 = go.Figure()
        for name, gg, col in [(s, g, TAB10[0]), ("RC", r, TAB10[1]), ("median " + str(k) + "-mer", med, "#888")]:
            after, before = S.bg_preds(m, gg)
            if after is None:
                continue
            d = after - before
            f2.add_box(y=d, name=name, boxpoints="all", jitter=0.5, pointpos=0, marker=dict(size=3, color=col),
                       line=dict(color=col, width=1), fillcolor="rgba(0,0,0,0)",
                       hovertemplate="%{y:.3f}<extra>background %{pointNumber}</extra>")
        f2.add_hline(y=0, line_width=1, line_color="#bbb")
        if S.bg_preds(m, g)[0] is None:
            f2.add_annotation(text="per-background values not stored for this k-mer", showarrow=False, x=0.5, y=1,
                              xref="paper", yref="paper", font_color="#999")
        f2.update_layout(yaxis_title="with k-mer - without", showlegend=False)

        # single-base edits (difference in value from the k-mer itself)
        E = seq_edits(s)
        ok = [e for e in E if len(e[3]) <= S.meta(m)["maxk"]]
        gv = encode([e[3] for e in ok])
        dv = v[gv] - v[g]
        zmax = float(np.nanmax(np.abs(dv))) or 1
        nrow = 3 if k < S.meta(m)["maxk"] else 2
        f3 = make_subplots(rows=nrow, cols=1, shared_xaxes=False, vertical_spacing=0.13,
                           row_heights=[4, 1, 4][:nrow], subplot_titles=["substitution", "deletion", "insertion"][:nrow])
        sub = np.full((4, k), np.nan)
        dele = np.full((1, k), np.nan)
        ins = np.full((4, k + 1), np.nan)
        for (kind_, i, b, ns), d in zip(ok, dv):
            if kind_ == "sub":
                sub["ACGT".index(b), i] = d
            elif kind_ == "del":
                dele[0, i] = d
            else:
                ins["ACGT".index(b), i] = d
        hm = dict(colorscale="RdBu_r", zmin=-zmax, zmax=zmax, xgap=1, ygap=1,
                  hovertemplate="%{y} at %{x}: %{z:.3f}<extra></extra>")
        f3.add_heatmap(z=sub, x=[f"{i + 1}{c}" for i, c in enumerate(s)], y=list("ACGT"), row=1, col=1,
                       colorbar=dict(title=dict(text="Δ " + lab, side="right"), thickness=10), **hm)
        f3.add_heatmap(z=dele, x=[f"{i + 1}{c}" for i, c in enumerate(s)], y=["del"], row=2, col=1, showscale=False, **hm)
        if nrow == 3:
            f3.add_heatmap(z=ins, x=[f"before {i + 1}" for i in range(k + 1)], y=list("ACGT"), row=3, col=1,
                           showscale=False, **hm)
            f3.update_xaxes(showticklabels=False, row=3, col=1)
        for i_ in range(1, nrow + 1):
            f3.update_yaxes(autorange="reversed", ticks="", row=i_, col=1)
            f3.update_xaxes(ticks="", row=i_, col=1, showticklabels=i_ == 2)
        f3.update_annotations(font_size=11, x=0, xanchor="left")
        f3.update_layout(margin=dict(t=24))

        # neighbours in the Nucleomer graph
        nb = S.neighbours(m, g)
        data, cols = [], []
        if len(nb):
            nb["seq"] = S.decode_np(nb.nbr.to_numpy())
            nb["k"] = k_of(nb.nbr.to_numpy())
            nb["dv"] = v[nb.nbr.to_numpy()] - v[g]
            nb["np"] = np.asarray(S.arr(m, "pval"))[nb.nbr.to_numpy()]
            nb = nb.sort_values(["pval", "dv"], ascending=[True, False])
            data = [dict(nbr=x.seq, k=int(x.k), dir=x.dir, edge_p=S.fmt_p(x.pval), dvalue=round(float(x.dv), 4),
                         node_p=S.fmt_p(x.np)) for x in nb.itertuples()]
            cols = [{"name": n, "id": i} for n, i in [("neighbour", "nbr"), ("k", "k"), ("dir", "dir"), ("edge p", "edge_p"),
                                                      ("Δ " + lab, "dvalue"), ("node p", "node_p")]]
        return info, f1, f2, f3, data, cols

    @app.callback(Output("kmer", "value", allow_duplicate=True),
                  Input("km-nbrs", "active_cell"), State("km-nbrs", "derived_viewport_data"), prevent_initial_call=True)
    def pick(cell, rows):
        return rows[cell["row"]]["nbr"] if cell and rows else no_update

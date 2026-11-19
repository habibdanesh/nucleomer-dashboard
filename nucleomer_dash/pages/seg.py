"""Segmentation: eval summary, PR/ROC curves, and a chr1 locus browser (eval tracks + live k-mer annotation)."""
import re
import numpy as np
import plotly.graph_objects as go
from dash import Input, Output, State, ctx, html, no_update
from plotly.subplots import make_subplots

from .. import store as S
from ..kmers import encode
from ..theme import BASE_COLOR, METHOD_COLOR, TAB10, empty, key
from ..ui import box, check, ctl, dd, fig, page, radio, row, table, txt

METHODS = ["Nucleomer", "FIMO", "Gradient", "DeepLIFT", "ISM", "Model"]


def layout():
    sm = S.seg_meta()
    tfs = list(sm["tfs"])
    side = [
        ctl("Metric", radio("sg-metric", [{"label": "AUPR", "value": "aupr"}, {"label": "AUROC", "value": "auc_roc_genome"},
                                          {"label": "AUPR / Model", "value": "rel"}], "rel"),
            "base-pair level, against ChIP-seq peaks on chr1. AUPR / Model: relative to the model's own predictions"),
        ctl("Methods", dcc_check("sg-methods", [m for m in METHODS if m in set(S.seg_summary().get("method", []))],
                                 ["Nucleomer", "FIMO", "Gradient", "DeepLIFT"])),
        html.Hr(),
        ctl("TF", dd("sg-tf", tfs, tfs[0] if tfs else None), "TF for the curves and the locus browser"),
        ctl("Locus", txt("sg-locus", "", placeholder="chr1:start-end"), "type a locus, e.g. chr1:1,000,000-1,001,000 (chr1 only)"),
        html.Div([html.Button("‹", id="sg-left", className="btn sm", title="move left"),
                  html.Button("−", id="sg-out", className="btn sm", title="zoom out"),
                  html.Button("+", id="sg-in", className="btn sm", title="zoom in"),
                  html.Button("›", id="sg-right", className="btn sm", title="move right"),
                  html.Button("peak", id="sg-peak", className="btn sm", title="jump to a random ChIP-seq peak")], className="btns"),
        ctl("Annotation k", html.Div([dcc_num("sg-kmin", 6), dcc_num("sg-kmax", 11)], className="pair"),
            "lengths of the k-mers drawn in the k-mers track"),
        check("sg-ann", "annotate significant k-mers", True,
              "k-mers track: significant k-mers of the TF's model in the window, strongest first, non-overlapping (windows <= 5 kb)"),
    ]
    return page(side,
                row(box("Segmentation accuracy (chr1)", fig("sg-bars", 300),
                        tip="Best run of each method per TF (eval v7). Dashed = random baseline (fraction of chr1 in peaks "
                            "for AUPR, 0.5 for AUROC); green = the model itself."),
                    box("PR and ROC curves", fig("sg-curves", 300), tip="For the TF chosen on the left; colours as in the bars.")),
                box("Locus browser", fig("sg-locus-fig", 520),
                    tip="ChIP-seq signal (10 bp bins), peaks, each method's scored segments (height = score), significant "
                        "k-mers of the TF's model (click to select) and the sequence (windows <= 300 bp)."),
                box("Values", table("sg-tab", page_size=8), tip="The numbers behind the bar plot (AUPR when 'AUPR / Model' is chosen)."),
                about="Segmentation eval (chr1) and a genome browser for its tracks.")


def dcc_check(id, options, value):
    from dash import dcc
    return dcc.Checklist(id=id, options=options, value=value, className="check", inline=True, persistence=True,
                         persistence_type="session")


def dcc_num(id, v):
    from ..ui import num
    return num(id, v, min=1, max=11)


def parse_locus(s):
    m = re.match(r"^\s*(chr\w+)\s*[:\s]\s*([\d,]+)\s*[-\s]\s*([\d,]+)\s*$", s or "")
    if not m:
        return None
    a, b = int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))
    return (m.group(1), min(a, b), max(a, b)) if b != a else None


def fmt(c, a, b):
    return f"{c}:{a:,}-{b:,}"


def rects(start, end, y, y0=0.0):
    x = np.c_[start, start, end, end, np.full(len(start), np.nan)].ravel()
    yy = np.c_[np.full(len(start), y0), y, y, np.full(len(start), y0), np.full(len(start), np.nan)].ravel()
    return x, yy


def annotate(seq, start, model, kmin, kmax, pthr):
    """Greedy non-overlapping significant k-mers in the window (longest-first ties broken by score)."""
    v, p = S.value(model, "score"), np.asarray(S.arr(model, "pval"))
    kmax = min(kmax, S.meta(model)["maxk"])
    cand = []
    for k in range(kmin, kmax + 1):
        subs = [seq[i:i + k] for i in range(len(seq) - k + 1)]
        ok = [i for i, s in enumerate(subs) if set(s) <= set("ACGT")]
        if not ok:
            continue
        g = encode([subs[i] for i in ok])
        for i, gg in zip(ok, g):
            if p[gg] < pthr:
                cand.append((float(v[gg]), k, i, subs[i]))
    cand.sort(key=lambda t: -t[0])
    taken = np.zeros(len(seq), bool)
    out = []
    for sc, k, i, s in cand:
        if not taken[i:i + k].any():
            taken[i:i + k] = True
            out.append((start + i, start + i + k, sc, s))
    return out


def register(app):
    @app.callback(Output("sg-bars", "figure"), Output("sg-tab", "data"), Output("sg-tab", "columns"),
                  Input("sg-metric", "value"), Input("sg-methods", "value"))
    def bars(metric, methods):
        d = S.seg_summary()
        if d.empty:
            return empty("no eval results prepared"), [], []
        d = d.sort_values("aupr").groupby(["model", "method"]).tail(1)
        tfs = list(S.seg_meta()["tfs"])
        ref = d[d.method == "Model"].set_index("model")["aupr"]
        f = go.Figure()
        ms = [m for m in METHODS if m in (methods or [])]
        for m in ms:
            x = d[d.method == m].set_index("model").reindex(tfs)
            y = x["aupr"] / ref.reindex(tfs) if metric == "rel" else x[metric]
            f.add_bar(x=tfs, y=y, name=m, marker_color=METHOD_COLOR.get(m, TAB10[7]),
                      hovertemplate="%{x} " + m + ": %{y:.3g}<extra></extra>")
        base = d.drop_duplicates("model").set_index("model").reindex(tfs)
        rnd = base["peak_bp_frac"] / ref.reindex(tfs) if metric == "rel" else (
            base["peak_bp_frac"] if metric == "aupr" else base["auroc_random"])
        for i, (t, r) in enumerate(zip(tfs, rnd)):
            f.add_shape(type="line", x0=i - 0.42, x1=i + 0.42, y0=r, y1=r, line=dict(dash="dash", width=1.2, color="black"))
        key(f, "random", "black", dash="dash")
        if metric == "rel":
            f.add_hline(y=1, line_color="darkgreen", line_width=1)
            key(f, "Model", "darkgreen", dash="solid")
        f.update_layout(barmode="group", bargap=0.25, legend=dict(orientation="h", y=1.12, x=0),
                        yaxis_title={"aupr": "AUPR", "auc_roc_genome": "AUROC", "rel": "AUPR / Model AUPR"}[metric])
        t = d[d.method.isin(ms + ["Model"])].copy()
        t["rel"] = t["aupr"] / t.model.map(ref)
        t = t.pivot(index="model", columns="method", values="aupr" if metric == "rel" else metric).reindex(tfs)
        rows = [{"model": tf, **{m: (None if np.isnan(t.loc[tf, m]) else round(float(t.loc[tf, m]), 4)) for m in t.columns}} for tf in tfs]
        cols = [{"name": "TF", "id": "model"}] + [{"name": m, "id": m} for m in METHODS if m in t.columns]
        return f, rows, cols

    @app.callback(Output("sg-curves", "figure"), Input("sg-tf", "value"), Input("sg-methods", "value"))
    def curves(tf, methods):
        if not tf:
            return empty()
        f = make_subplots(rows=1, cols=2, subplot_titles=["PR", "ROC (genome)"], horizontal_spacing=0.14)
        for m in [m for m in METHODS if m in (methods or []) + ["Model"]]:
            p = S.DATA / "seg" / tf / f"{m}.curves.npz"
            if not p.exists():
                continue
            z = np.load(p)
            c = METHOD_COLOR.get(m, TAB10[7])
            f.add_scatter(x=z["pr"][:, 0], y=z["pr"][:, 1], mode="lines", name=m, line=dict(color=c, width=1.3), row=1, col=1)
            f.add_scatter(x=z["roc"][:, 0], y=z["roc"][:, 1], mode="lines", name=m, line=dict(color=c, width=1.3),
                          showlegend=False, row=1, col=2)
        f.update_xaxes(title_text="recall", row=1, col=1)
        f.update_yaxes(title_text="precision", row=1, col=1)
        f.update_xaxes(title_text="FPR", row=1, col=2)
        f.update_yaxes(title_text="TPR", row=1, col=2)
        f.update_annotations(font_size=12)
        f.update_layout(showlegend=False, title=None)
        return f

    @app.callback(Output("sg-locus", "value"),
                  Input("sg-tf", "value"), Input("sg-left", "n_clicks"), Input("sg-right", "n_clicks"),
                  Input("sg-in", "n_clicks"), Input("sg-out", "n_clicks"), Input("sg-peak", "n_clicks"),
                  State("sg-locus", "value"))
    def nav(tf, *_):
        cur = parse_locus(_[-1])
        t = ctx.triggered_id
        if t in ("sg-peak", "sg-tf") or cur is None:
            if t == "sg-tf" and cur is not None:
                return no_update
            pk = S.seg_track(tf, "peaks")
            i = int(np.argmax(pk["signal"])) if t != "sg-peak" else int(np.random.default_rng().integers(len(pk["start"])))
            c = (int(pk["start"][i]) + int(pk["end"][i])) // 2
            return fmt("chr1", c - 250, c + 250)
        c, a, b = cur
        w = b - a
        if t == "sg-left":
            a, b = a - w // 2, b - w // 2
        elif t == "sg-right":
            a, b = a + w // 2, b + w // 2
        elif t == "sg-in":
            a, b = a + w // 4, b - w // 4
        elif t == "sg-out":
            a, b = a - w // 2, b + w // 2
        return fmt(c, max(0, a), max(a + 10, b))

    @app.callback(Output("sg-locus-fig", "figure"),
                  Input("sg-locus", "value"), Input("sg-tf", "value"), Input("sg-methods", "value"), Input("sg-ann", "value"),
                  Input("sg-kmin", "value"), Input("sg-kmax", "value"), Input("pthr", "value"))
    def locus(loc, tf, methods, ann, kmin, kmax, pthr):
        L = parse_locus(loc)
        if not tf or L is None:
            return empty("locus as chr1:start-end")
        c, a, b = L
        if b - a > 200_000:
            return empty("window > 200 kb")
        tm = S.seg_meta()["tfs"][tf]
        ms = [m for m in METHODS if m in (methods or []) + ["Model"] and any(x["method"] == m for x in tm["methods"])]
        model_ok = tf in S.models()
        show_ann = bool(ann) and model_ok and b - a <= 5000
        show_seq = b - a <= 300
        tracks = (["ChIP-seq"] if tm.get("chip") else []) + ["peaks"] + ms + (["k-mers"] if show_ann else []) + (["seq"] if show_seq else [])
        hts = [1.2 if t in ("ChIP-seq", "k-mers") else (0.35 if t in ("peaks", "seq") else 1) for t in tracks]
        f = make_subplots(rows=len(tracks), cols=1, shared_xaxes=True, vertical_spacing=0.02, row_heights=hts)
        for i, t in enumerate(tracks, 1):
            if t == "ChIP-seq":
                nb = int(min(1000, b - a))
                y = S.chip(tf, a, b, nb)
                x = np.linspace(a, b, nb + 1)
                f.add_scatter(x=np.repeat(x, 2)[1:-1], y=np.repeat(y, 2), mode="lines", fill="tozeroy", line=dict(color="#555", width=0.8),
                              fillcolor="rgba(85,85,85,0.35)", hoverinfo="skip", row=i, col=1)
            elif t == "peaks":
                w = S.seg_window(tf, "peaks", a, b)
                x, y = rects(w["start"], w["end"], 1)
                f.add_scatter(x=x, y=y, mode="lines", fill="toself", line=dict(width=0), fillcolor="#333", hoverinfo="skip", row=i, col=1)
            elif t == "k-mers":
                seq = S.fetch(c, a, b)
                hits = annotate(seq, a, tf, int(kmin or 6), int(kmax or 11), pthr or S.PVAL)
                if hits:
                    st, en, sc, ss = map(np.array, zip(*hits))
                    x, y = rects(st, en, sc)
                    f.add_scatter(x=x, y=y, mode="lines", fill="toself", line=dict(width=0.5, color=TAB10[0]),
                                  fillcolor="rgba(31,119,180,0.5)", hoverinfo="skip", row=i, col=1)
                    f.add_scatter(x=(st + en) / 2, y=sc, mode="markers+text" if b - a <= 250 else "markers",
                                  text=ss, textposition="top center", textfont=dict(size=9, family="Menlo, monospace"),
                                  marker=dict(size=4, color=TAB10[0]), customdata=ss,
                                  hovertemplate="%{customdata}<br>score %{y:.2f}<extra></extra>", row=i, col=1)
            elif t == "seq":
                seq = S.fetch(c, a, b)
                for base, col in BASE_COLOR.items():
                    pos = [a + j + 0.5 for j, ch in enumerate(seq) if ch == base]
                    f.add_scatter(x=pos, y=[0.5] * len(pos), mode="text", text=[base] * len(pos),
                                  textfont=dict(color=col, size=10 if b - a <= 150 else 7, family="Menlo, monospace"),
                                  hoverinfo="skip", row=i, col=1)
            else:
                w = S.seg_window(tf, t, a, b)
                x, y = rects(w["start"], w["end"], w["score"])
                col = METHOD_COLOR.get(t, TAB10[7])
                f.add_scatter(x=x, y=y, mode="lines", fill="toself", line=dict(width=0.6, color=col), fillcolor=col,
                              opacity=0.7, hoverinfo="skip", row=i, col=1)
            f.update_yaxes(title_text=t, title_standoff=4, title_font_size=10, row=i, col=1, nticks=3,
                           showticklabels=t not in ("peaks", "seq"), tickfont_size=9)
        f.update_xaxes(range=[a, b], showticklabels=False)
        f.update_xaxes(showticklabels=True, tickformat=",d", row=len(tracks), col=1)
        f.update_layout(showlegend=False, margin=dict(l=80, r=12, t=20, b=30), title=dict(text=f"{tf}  {fmt(c, a, b)}  ({b - a:,} bp)"))
        return f

    @app.callback(Output("kmer", "value", allow_duplicate=True), Input("sg-locus-fig", "clickData"), prevent_initial_call=True)
    def pick(cd):
        if cd and cd["points"] and isinstance(cd["points"][0].get("customdata"), str):
            return cd["points"][0]["customdata"]
        return no_update

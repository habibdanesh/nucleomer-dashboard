"""Map: 2D/3D embedding of the top k-mers of one or more models (plus optional random decoys).

Distance between two k-mers = number of single-base insertions/deletions along a path that stays inside the
selected set (either Nucleomer edges significant in any selected model, or every one-indel pair). Unconnected
pairs get (max finite distance + 1). Embedded with UMAP on the precomputed distances, or a spring layout.
"""
import threading
from functools import lru_cache

import numpy as np
import plotly.graph_objects as go
from dash import Input, Output, State, html, no_update
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path

from .. import store as S
from ..kmers import OFF, encode, k_of, revcomp_gid
from ..theme import TAB10, empty
from ..ui import box, check, ctl, dd, fig, num, page, radio, rng, txt


def layout():
    ms = S.models()
    side = [
        ctl("Models", dd("em-models", ms, ["GATA2"] if "GATA2" in ms else ms[:1], multi=True),
            "k-mers are pooled from all chosen models"),
        ctl("k", rng("em-k", 4, 11, [6, 11])),
        ctl("Top per k per model", num("em-top", 100, min=5, max=3000, step=5),
            "the strongest significant k-mers of each length, from each model"),
        ctl("Random decoys", num("em-decoy", 0, min=0, max=5000, step=100),
            "random k-mers with the same length distribution, as a background (grey)"),
        check("em-rc", "merge reverse complements", True, "a k-mer and its reverse complement are one point"),
        ctl("Distance", radio("em-dist", [{"label": "graph edges", "value": "edges"},
                                          {"label": "any 1-indel", "value": "indel"}], "edges"),
            "distance = number of single-base insertions/deletions between two k-mers, through significant graph edges "
            "of the chosen models (graph edges) or through any k-mer pair one indel apart (any 1-indel)"),
        ctl("Embedding", radio("em-meth", ["UMAP", "spring"], "UMAP")),
        ctl("View", radio("em-dim", ["2D", "3D"], "2D")),
        ctl("UMAP neighbours / min dist", html.Div([num("em-nn", 15, min=2, max=200), num("em-md", 0.3, min=0, max=1, step=0.05)],
                                                   className="pair"), "UMAP n_neighbors and min_dist"),
        html.Button("Run", id="em-run", n_clicks=0, className="btn", title="recompute the map with the settings above"),
        html.Hr(),
        ctl("Colour", dd("em-col", [{"label": "best model", "value": "best"},
                                    {"label": "value, one panel per model", "value": "facet"},
                                    {"label": "# models significant", "value": "nsig"},
                                    {"label": "k", "value": "k"}, {"label": "pattern", "value": "pattern"},
                                    {"label": "decoy", "value": "decoy"}] +
                         [{"label": f"score in {m}", "value": "m:" + m} for m in ms], "best"),
            "changes colours only (no Run needed). best model: the model in which the k-mer scores highest"),
        ctl("Pattern", txt("em-pat", "GATA", placeholder="regex, either strand"), "used by Colour = pattern"),
        check("em-edges", "draw edges", False),
        html.Div(id="em-stat", className="small muted"),
    ]
    return page(side, box("Map of top k-mers", fig("em-fig", 720, cls="fig"),
                          tip="Each point is a k-mer; nearby points are few insertions/deletions apart. Point size = best "
                              "score over the chosen models; red ring = selected k-mer. Click a point to select it."),
                about="Top k-mers of one or more models, laid out by edit distance. Set the options, then Run.")


def pick_kmers(models, k0, k1, top, merge, pthr):
    sel = set()
    for m in models:
        v, p = S.value(m, "score"), np.asarray(S.arr(m, "pval"))
        for k in range(k0, min(k1, S.meta(m)["maxk"]) + 1):
            g = np.arange(OFF[k], OFF[k + 1], dtype=np.int64)
            g = g[p[g] < pthr]
            if merge:
                c = S.canonical(g)
                vv = np.maximum(v[g], v[revcomp_gid(g)])
                o = np.argsort(-vv, kind="stable")
                g = c[o][np.sort(np.unique(c[o], return_index=True)[1])]
            else:
                g = g[np.argsort(-v[g], kind="stable")]
            sel.update(g[:top].tolist())
    return np.array(sorted(sel), np.int64)


_LOCK = threading.Lock()   # numba (UMAP) is not safe to run from two request threads at once


def compute(*args):
    with _LOCK:
        return _compute(*args)


@lru_cache(16)
def _compute(models, k0, k1, top, decoys, merge, dist, meth, dim, nn, md, pthr):
    g = pick_kmers(models, k0, k1, top, merge, pthr)
    n_sel = len(g)
    if decoys:
        rng_ = np.random.default_rng(0)
        pool = np.arange(OFF[k0], OFF[k1 + 1], dtype=np.int64)
        # match the length distribution of the selected k-mers
        kk = k_of(g)
        d = np.concatenate([rng_.choice(np.arange(OFF[k], OFF[k + 1]), size=int(round(decoys * (kk == k).mean())), replace=False)
                            for k in range(k0, k1 + 1) if (kk == k).any()])
        d = S.canonical(d) if merge else d
        d = np.setdiff1d(np.unique(d), g)
        g = np.r_[g, d]
    is_decoy = np.r_[np.zeros(n_sel, bool), np.ones(len(g) - n_sel, bool)]
    n = len(g)
    idx = {int(x): i for i, x in enumerate(g)}
    canon = S.canonical if merge else (lambda a: np.asarray(a, np.int64))
    rows, cols = [], []
    if dist == "edges":
        for m in models:
            E = S.edges(m)
            for st in (g, revcomp_gid(g)) if merge else (g,):
                for i, x in enumerate(st):
                    o = slice(E["out_ptr"][x], E["out_ptr"][x + 1])
                    d = E["dst"][o][E["pval"][o] < pthr]
                    for y in canon(d):
                        j = idx.get(int(y))
                        if j is not None and j != i:
                            rows.append(i)
                            cols.append(j)
    else:
        seqs = S.decode_np(g)
        for i, s in enumerate(seqs):
            if len(s) < 2:
                continue
            dels = encode([s[:p] + s[p + 1:] for p in range(len(s))])
            for y in canon(dels):
                j = idx.get(int(y))
                if j is not None:
                    rows.append(i)
                    cols.append(j)
    A = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    A = ((A + A.T) > 0).astype(float)
    ei, ej = A.nonzero()
    up = ei < ej
    if meth == "UMAP":
        D = shortest_path(A, unweighted=True, directed=False)
        fin = D[np.isfinite(D)].max() if np.isfinite(D).any() else 1
        D[~np.isfinite(D)] = fin + 1
        import inspect
        import warnings
        import umap
        # umap-learn >= 0.6 changes the default layout; keep the 0.5.x one so maps match across installs
        kw = {"compatibility_layout": True} if "compatibility_layout" in inspect.signature(umap.UMAP).parameters else {}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)   # "precomputed metric: no inverse_transform", "n_jobs overridden"
            X = umap.UMAP(metric="precomputed", n_neighbors=min(nn, n - 1), min_dist=md, n_components=3 if dim == "3D" else 2,
                          random_state=0, **kw).fit_transform(D)
    else:
        import networkx as nx
        G = nx.from_scipy_sparse_array(A)
        pos = nx.spring_layout(G, dim=3 if dim == "3D" else 2, seed=0, iterations=200, k=1.5 / np.sqrt(n))
        X = np.array([pos[i] for i in range(n)])
    return g, is_decoy, X, ei[up], ej[up]


def facet(X, V, models, kind, size, seqs, hov, sc, dec, ei, ej, draw_edges):
    """Same layout, one panel per model, coloured by that model's value (shared colour scale)."""
    from plotly.subplots import make_subplots
    n = len(models)
    nc = min(n, 3)
    nr = int(np.ceil(n / nc))
    f = make_subplots(rows=nr, cols=nc, subplot_titles=list(models), horizontal_spacing=0.02, vertical_spacing=0.06)
    lo, hi = np.nanpercentile(V[~dec], [1, 99]) if (~dec).any() else (0, 1)
    for j, m in enumerate(models):
        r, c = j // nc + 1, j % nc + 1
        if draw_edges and len(ei):
            f.add_scattergl(x=np.c_[X[ei, 0], X[ej, 0], np.full(len(ei), np.nan)].ravel(),
                            y=np.c_[X[ei, 1], X[ej, 1], np.full(len(ei), np.nan)].ravel(), mode="lines",
                            line=dict(color="#ddd", width=0.5), hoverinfo="skip", showlegend=False, row=r, col=c)
        o = np.argsort(np.nan_to_num(V[:, j], nan=-np.inf))
        f.add_scattergl(x=X[o, 0], y=X[o, 1], mode="markers", showlegend=False, customdata=seqs[o],
                        hovertext=np.array(hov)[o], hoverinfo="text",
                        marker=dict(size=size[o] * 0.8, color=V[o, j], coloraxis="coloraxis", line=dict(width=0)),
                        row=r, col=c)
    f.update_xaxes(visible=False)
    f.update_yaxes(visible=False)
    for j in range(n):
        f.update_yaxes(scaleanchor=f"x{j + 1 if j else ''}", row=j // nc + 1, col=j % nc + 1)
    f.update_layout(coloraxis=dict(colorscale="Viridis", cmin=lo, cmax=hi,
                                   colorbar=dict(title=dict(text=S.VALUE_LABEL[kind], side="right"), thickness=10)))
    f.update_annotations(font_size=12)
    return f


def register(app):
    @app.callback(Output("em-fig", "figure"), Output("em-stat", "children"),
                  Input("em-run", "n_clicks"), Input("em-col", "value"), Input("em-pat", "value"), Input("em-edges", "value"),
                  Input("value", "value"),
                  State("em-models", "value"), State("em-k", "value"), State("em-top", "value"), State("em-decoy", "value"),
                  State("em-rc", "value"), State("em-dist", "value"), State("em-meth", "value"), State("em-dim", "value"),
                  State("em-nn", "value"), State("em-md", "value"), State("pthr", "value"), State("kmer", "value"))
    def update(_, col, pat, draw_edges, kind, models, krange, top, decoys, merge, dist, meth, dim, nn, md, pthr, cur):
        if not models:
            return empty("choose models"), ""
        models = tuple(models)
        try:
            g, dec, X, ei, ej = compute(models, krange[0], krange[1], int(top or 100), int(decoys or 0), bool(merge), dist,
                                        meth, dim, int(nn or 15), float(md if md is not None else 0.3), float(pthr or S.PVAL))
        except Exception as e:
            return empty(f"failed: {e}"), ""
        seqs = S.decode_np(g)
        kk = k_of(g)
        # value per model (stronger strand when merged)
        V = np.stack([np.maximum(S.value(m, kind)[g], S.value(m, kind)[revcomp_gid(g)]) if merge else S.value(m, kind)[g]
                      for m in models], 1)
        Pv = np.stack([np.minimum(np.asarray(S.arr(m, "pval"))[g], np.asarray(S.arr(m, "pval"))[revcomp_gid(g)]) if merge
                       else np.asarray(S.arr(m, "pval"))[g] for m in models], 1)
        sc = np.nanmax(np.stack([np.maximum(S.value(m, "score")[g], S.value(m, "score")[revcomp_gid(g)]) for m in models], 1), 1)
        size = np.where(dec, 3, 3 + 12 * np.clip(np.nan_to_num(sc), 0, 1))
        hov = ["<br>".join([s] + [f"{m} {V[i, j]:.3f}" for j, m in enumerate(models)]) for i, s in enumerate(seqs)]
        d3 = X.shape[1] == 3
        if col == "facet" and not d3:
            return facet(X, V, models, kind, size, seqs, hov, sc, dec, ei, ej, draw_edges), \
                f"{int((~dec).sum())} k-mers + {int(dec.sum())} decoys, {len(ei)} edges"
        Sc = go.Scatter3d if d3 else go.Scattergl
        f = go.Figure()
        xyz = lambda idx: dict(x=X[idx, 0], y=X[idx, 1], **({"z": X[idx, 2]} if d3 else {}))
        if draw_edges and len(ei):
            ex = np.c_[X[ei, 0], X[ej, 0], np.full(len(ei), np.nan)].ravel()
            ey = np.c_[X[ei, 1], X[ej, 1], np.full(len(ei), np.nan)].ravel()
            kw = dict(x=ex, y=ey, **({"z": np.c_[X[ei, 2], X[ej, 2], np.full(len(ei), np.nan)].ravel()} if d3 else {}))
            f.add_trace(Sc(**kw, mode="lines", line=dict(color="#ccc", width=0.6 if not d3 else 1), hoverinfo="skip", showlegend=False))
        o = np.argsort(np.nan_to_num(sc, nan=-1))       # strongest on top
        cat = None
        if col == "best":
            cat = np.array(models)[np.nanargmax(np.nan_to_num(V, nan=-np.inf), 1)]
            cat = np.where(dec, "decoy", cat)
            order = list(models) + ["decoy"]
        elif col == "decoy":
            cat, order = np.where(dec, "decoy", "selected"), ["selected", "decoy"]
        elif col == "pattern":
            hit = np.isin(g, S.search(g, pat)) if pat else np.zeros(len(g), bool)
            cat, order = np.where(hit, pat or "", "other"), [pat or "", "other"]
        if cat is not None:
            for i, c in enumerate(order):
                idx = o[cat[o] == c]
                if not len(idx):
                    continue
                color = "#cfcfcf" if c in ("decoy", "other") else TAB10[i % 10]
                f.add_trace(Sc(**xyz(idx), mode="markers", name=c,
                               marker=dict(size=size[idx] * (0.6 if d3 else 1), color=color, line=dict(width=0)),
                               customdata=seqs[idx], hovertext=np.array(hov)[idx], hoverinfo="text"))
        else:
            if col == "nsig":
                c, ct = (Pv < (pthr or S.PVAL)).sum(1), "# models"
            elif col == "k":
                c, ct = kk, "k"
            else:
                m = col[2:]
                v = S.value(m, kind)
                c, ct = (np.maximum(v[g], v[revcomp_gid(g)]) if merge else v[g]), f"{S.VALUE_LABEL[kind]} ({m})"
            f.add_trace(Sc(**xyz(o), mode="markers", showlegend=False,
                           marker=dict(size=size[o] * (0.6 if d3 else 1), color=c[o], colorscale="Viridis", showscale=True,
                                       colorbar=dict(title=dict(text=ct, side="right"), thickness=10), line=dict(width=0)),
                           customdata=seqs[o], hovertext=np.array(hov)[o], hoverinfo="text"))
        gc = S.parse_kmer(cur)
        if gc is not None:
            hit = np.where((g == gc) | (g == revcomp_gid([gc])[0]))[0]
            if len(hit):
                f.add_trace(Sc(**xyz(hit), mode="markers", showlegend=False, hoverinfo="skip",
                               marker=dict(size=14, color="rgba(0,0,0,0)", line=dict(width=2, color="#d62728"))))
        if d3:
            f.update_layout(scene=dict(xaxis_visible=False, yaxis_visible=False, zaxis_visible=False), margin=dict(l=0, r=0, t=0, b=0))
        else:
            f.update_layout(xaxis_visible=False, yaxis_visible=False, yaxis_scaleanchor="x")
        f.update_layout(legend=dict(x=1, xanchor="right", y=1, title_text={"best": "best model", "decoy": "", "pattern": "pattern"}.get(col, "")),
                        uirevision=str((models, krange, top, decoys, dist, meth, dim)))
        stat = f"{int((~dec).sum())} k-mers + {int(dec.sum())} decoys, {len(ei)} edges"
        return f, stat

    @app.callback(Output("kmer", "value", allow_duplicate=True), Input("em-fig", "clickData"), prevent_initial_call=True)
    def pick(c):
        if c and c["points"] and c["points"][0].get("customdata"):
            return c["points"][0]["customdata"]
        return no_update

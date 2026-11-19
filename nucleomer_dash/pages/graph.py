"""Graph: the Nucleomer graph around the selected k-mer (breadth-first, edges filtered by p-value)."""
import networkx as nx
import numpy as np
import plotly.graph_objects as go
from dash import Input, Output, no_update

from .. import store as S
from ..kmers import k_of
from ..theme import empty
from ..ui import box, check, ctl, dd, fig, num, page, radio, slider


def layout():
    ms = S.models()
    side = [
        ctl("Radius", slider("gr-rad", 1, 4, 2), "how many edges away from the selected k-mer to go"),
        ctl("Max nodes", num("gr-max", 300, min=10, max=3000, step=10),
            "if the neighbourhood is larger, the strongest k-mers (by score) are kept"),
        check("gr-sig", "significant edges only", True, "only edges with edge p below the p threshold in the top bar"),
        ctl("Layout", radio("gr-lay", [{"label": "k levels", "value": "lattice"}, {"label": "spring", "value": "spring"}], "lattice"),
            "k levels: one row per k-mer length (shorter at the bottom); spring: force-directed"),
        ctl("View", radio("gr-dim", ["2D", "3D"], "2D")),
        ctl("Colour", dd("gr-col", [{"label": "value", "value": "value"}, {"label": "-log10 p", "value": "nlp"},
                                     {"label": "k", "value": "k"}] +
                         [{"label": f"value in {m}", "value": "m:" + m} for m in ms], "value"),
            "'value in X' colours the same graph by another model's values"),
        check("gr-lab", "labels", False),
    ]
    return page(side, box("Nucleomer graph around the selected k-mer", fig("gr-fig", 680),
                          tip="Nodes are k-mers; edges join k-mers one insertion or deletion apart. Node size = score "
                              "(model's best k-mer = largest); red ring = selected k-mer. Click a node to move to it."),
                about="The model's Nucleomer graph, starting from the selected k-mer.")


def ego(m, g, radius, max_nodes, sig, pthr):
    E = S.edges(m)
    v = S.value(m, "score")

    def nbrs(x):
        o = slice(E["out_ptr"][x], E["out_ptr"][x + 1])
        i = E["in_order"][E["in_ptr"][x]:E["in_ptr"][x + 1]]
        a, pa = E["dst"][o], E["pval"][o]
        b, pb = E["src"][i], E["pval"][i]
        n = np.r_[a, b]
        return n[np.r_[pa, pb] < pthr] if sig else n

    nodes, frontier = {g}, [g]
    for _ in range(radius):
        nxt = set()
        for x in frontier:
            nxt.update(int(y) for y in nbrs(x))
        nxt -= nodes
        if len(nodes) + len(nxt) > max_nodes:
            nxt = sorted(nxt, key=lambda y: -v[y])[:max_nodes - len(nodes)]
            nodes.update(nxt)
            break
        nodes.update(nxt)
        frontier = nxt
    nodes = np.array(sorted(nodes), np.int64)
    G = nx.Graph()
    G.add_nodes_from(nodes.tolist())
    for x in nodes:
        o = slice(E["out_ptr"][x], E["out_ptr"][x + 1])
        d, p = E["dst"][o], E["pval"][o]
        keep = np.isin(d, nodes) & ((p < pthr) if sig else True)
        G.add_edges_from((int(x), int(y), {"p": float(pp)}) for y, pp in zip(d[keep], p[keep]))
    return G


def register(app):
    @app.callback(Output("gr-fig", "figure"),
                  Input("kmer", "value"), Input("model", "value"), Input("value", "value"), Input("pthr", "value"),
                  Input("gr-rad", "value"), Input("gr-max", "value"), Input("gr-sig", "value"), Input("gr-lay", "value"),
                  Input("gr-dim", "value"), Input("gr-col", "value"), Input("gr-lab", "value"))
    def update(s, m, kind, pthr, rad, mx, sig, lay, dim, col, lab):
        g = S.parse_kmer(s)
        if g is None or not m:
            return empty("enter a k-mer above")
        if len(s.strip()) > S.meta(m)["maxk"]:
            return empty(f"{m} was run up to k = {S.meta(m)['maxk']}")
        pthr = pthr or S.PVAL
        G = ego(m, g, rad or 2, mx or 300, bool(sig), pthr)
        nodes = np.array(list(G.nodes), np.int64)
        if len(nodes) == 1:
            return empty("no edges at this threshold")
        d3 = dim == "3D"
        pos = nx.spring_layout(G, dim=3 if d3 else 2, seed=0, iterations=100)
        P = np.array([pos[int(x)] for x in nodes])
        kk = k_of(nodes)
        if lay == "lattice":                       # one level per k-mer length
            P[:, -1] = kk
            if not d3:
                ix = {int(x): i for i, x in enumerate(nodes)}
                nb = [[ix[int(y)] for y in G.neighbors(int(x))] for x in nodes]
                for _ in range(4):                 # barycentre sweeps, then spread each level evenly
                    P[:, 0] = [np.mean(P[n, 0]) if n else P[i, 0] for i, n in enumerate(nb)]
                    for k_ in np.unique(kk):
                        i = np.where(kk == k_)[0]
                        P[i[np.argsort(P[i, 0], kind="stable")], 0] = np.linspace(-1, 1, len(i) + 2)[1:-1] if len(i) > 1 else 0
        idx = {int(x): i for i, x in enumerate(nodes)}
        seqs = S.decode_np(nodes)
        v = S.value(m, kind)[nodes]
        if col == "value":
            c, ct = v, S.VALUE_LABEL[kind]
        elif col == "nlp":
            c, ct = S.value(m, "nlp")[nodes], "-log10 p"
        elif col == "k":
            c, ct = kk, "k"
        else:
            mm = col[2:]
            c, ct = S.value(mm, kind)[nodes], f"{S.VALUE_LABEL[kind]} ({mm})"
        sc = np.clip(S.value(m, "score")[nodes], 0, 1)
        size = 5 + 14 * np.nan_to_num(sc)
        is_sel = nodes == g
        ex, ey, ez = [], [], []
        for a, b in G.edges:
            pa, pb = P[idx[a]], P[idx[b]]
            ex += [pa[0], pb[0], None]
            ey += [pa[1], pb[1], None]
            if d3:
                ez += [pa[2], pb[2], None]
        hov = [f"{s_}<br>k {k_}<br>{S.VALUE_LABEL[kind]} {v_:.3f}" for s_, k_, v_ in zip(seqs, kk, v)]
        mk = dict(size=size, color=c, colorscale="Viridis", showscale=True, opacity=0.9,
                  colorbar=dict(title=dict(text=ct, side="right"), thickness=10),
                  line=dict(width=np.where(is_sel, 3, 0.5), color=np.where(is_sel, "#d62728", "white")))
        mode = "markers+text" if lab else "markers"
        if d3:
            f = go.Figure(go.Scatter3d(x=ex, y=ey, z=ez, mode="lines", line=dict(color="#bbb", width=1), hoverinfo="skip"))
            f.add_scatter3d(x=P[:, 0], y=P[:, 1], z=P[:, 2], mode=mode, text=seqs if lab else None,
                            marker={**mk, "size": size * 0.6, "line": dict(width=0)}, customdata=seqs,
                            hovertext=hov, hoverinfo="text", textfont_size=9)
            f.update_layout(showlegend=False, scene=dict(xaxis_visible=False, yaxis_visible=False,
                                                         zaxis=dict(title="k" if lay == "lattice" else "", visible=lay == "lattice")),
                            margin=dict(l=0, r=0, t=0, b=0))
        else:
            f = go.Figure(go.Scattergl(x=ex, y=ey, mode="lines", line=dict(color="#ccc", width=0.7), hoverinfo="skip"))
            f.add_scattergl(x=P[:, 0], y=P[:, 1], mode=mode, text=seqs if lab else None, textposition="top center",
                            textfont_size=9, marker=mk, customdata=seqs, hovertext=hov, hoverinfo="text")
            f.update_layout(showlegend=False, xaxis_visible=False,
                            yaxis=dict(visible=lay == "lattice", title="k", dtick=1))
        f.update_layout(title=f"{len(nodes)} k-mers, {G.number_of_edges()} edges   (size = score, red ring = {s.strip().upper()})")
        return f

    @app.callback(Output("kmer", "value", allow_duplicate=True), Input("gr-fig", "clickData"), prevent_initial_call=True)
    def pick(c):
        if c and c["points"] and c["points"][0].get("customdata"):
            return c["points"][0]["customdata"]
        return no_update

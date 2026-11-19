"""Plot defaults: plain white axes, matplotlib's tab10 colours and viridis, small margins."""
import plotly.graph_objects as go
import plotly.io as pio

TAB10 = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
FAMILY_COLOR = {"BPNet": TAB10[0], "Cherimoya": TAB10[1], "ProCapNet": TAB10[2]}
METHOD_COLOR = {"Nucleomer": TAB10[0], "FIMO": TAB10[1], "Model": TAB10[2], "DeepLIFT": TAB10[3], "Gradient": TAB10[5],
                "ISM": TAB10[4]}
BASE_COLOR = {"A": "#2ca02c", "C": "#1f77b4", "G": "#ff7f0e", "T": "#d62728"}
FONT = "Helvetica Neue, Helvetica, Arial, sans-serif"

pio.templates["nm"] = go.layout.Template(pio.templates["simple_white"])
t = pio.templates["nm"].layout
t.font = dict(family=FONT, size=12, color="#222")
t.colorway = TAB10
t.margin = dict(l=56, r=12, t=28, b=44)
t.colorscale.sequential = "Viridis"
t.coloraxis.colorbar = dict(thickness=10, outlinewidth=0, len=0.8)
t.hoverlabel = dict(font_family=FONT, font_size=11, bgcolor="white")
t.legend = dict(font_size=11, itemsizing="constant")
t.xaxis.update(ticks="outside", showgrid=False, zeroline=False)
t.yaxis.update(ticks="outside", showgrid=False, zeroline=False)
t.title = dict(font_size=12, x=0.0, xanchor="left", xref="paper")
pio.templates.default = "nm"

GRAPH_CONFIG = {"displaylogo": False, "toImageButtonOptions": {"format": "svg", "scale": 1},
                "modeBarButtonsToRemove": ["select2d", "lasso2d", "autoScale2d", "sendChartToCloud"], "showSendToCloud": False}
GRAPH_CONFIG_SELECT = {**GRAPH_CONFIG, "modeBarButtonsToRemove": ["autoScale2d", "sendChartToCloud"]}


def key(f, name, color, symbol="square", dash=None, row=None, col=None):
    """Legend-only entry (for things that have no legend of their own: heatmaps, shapes, marker colours)."""
    kw = dict(row=row, col=col) if row else {}
    if dash:
        f.add_scatter(x=[None], y=[None], mode="lines", name=name, line=dict(color=color, dash=dash, width=1.2), **kw)
    else:
        f.add_scatter(x=[None], y=[None], mode="markers", name=name, marker=dict(color=color, symbol=symbol, size=9,
                      line=dict(width=0) if symbol != "square-open" else dict(width=1.5, color=color)), **kw)
    return f


def legend_top(f, **kw):
    """Horizontal legend above the plot."""
    f.update_layout(showlegend=True, legend=dict(orientation="h", x=0, xanchor="left", y=1.0, yanchor="bottom",
                                                 bgcolor="rgba(255,255,255,0)", **kw), margin_t=34)
    return f


def empty(msg=""):
    f = go.Figure()
    f.update_layout(xaxis_visible=False, yaxis_visible=False,
                    annotations=[dict(text=msg, showarrow=False, font_color="#999")] if msg else [])
    return f

"""Small layout helpers so every page looks the same."""
from dash import dcc, html, dash_table

from .theme import GRAPH_CONFIG, empty

P = dict(persistence=True, persistence_type="session")


def hint(text):
    """Small '?' that shows text on hover."""
    return html.Span("?", className="hint", **{"data-tip": text}) if text else None


def ctl(label, comp, tip=None):
    return html.Div([html.Label([label, hint(tip)]), comp], className="ctl")


def dd(id, options, value, multi=False, **kw):
    opts = [o if isinstance(o, dict) else {"label": str(o), "value": o} for o in options]
    return dcc.Dropdown(id=id, options=opts, value=value, multi=multi, clearable=multi, searchable=True, **P, **kw)


def radio(id, options, value):
    opts = [o if isinstance(o, dict) else {"label": str(o), "value": o} for o in options]
    return dcc.RadioItems(id=id, options=opts, value=value, inline=True, className="radio", **P)


def check(id, label, on=False, tip=None):
    return dcc.Checklist(id=id, options=[{"label": html.Span([label, hint(tip)]), "value": 1}], value=[1] if on else [],
                         className="check", **P)


def num(id, value, **kw):
    return dcc.Input(id=id, type="number", value=value, debounce=True, className="num", **P, **kw)


def txt(id, value="", placeholder="", **kw):
    return dcc.Input(id=id, type="text", value=value, placeholder=placeholder, debounce=True, className="txt", **P, **kw)


def rng(id, lo, hi, value, step=1, marks=None):
    return dcc.RangeSlider(id=id, min=lo, max=hi, step=step, value=value, allowCross=False,
                           marks=marks or {i: str(i) for i in range(lo, hi + 1)}, className="slider", **P)


def slider(id, lo, hi, value, step=1, marks=None):
    return dcc.Slider(id=id, min=lo, max=hi, step=step, value=value,
                      marks=marks or {i: str(i) for i in range(lo, hi + 1)}, className="slider", **P)


def fig(id, h=360, config=None, cls="fig"):
    return dcc.Loading(dcc.Graph(id=id, figure=empty(), config=config or GRAPH_CONFIG, style={"height": f"{h}px"}),
                       type="dot", color="#999", delay_show=300, className=cls)


def page(side, *main, about=None):
    side = ([html.Div(about, className="about")] if about else []) + list(side)
    return html.Div([html.Div(side, className="side"), html.Div(list(main), className="main")], className="page")


def row(*children, cls="row"):
    return html.Div(list(children), className=cls)


def box(title, *children, tip=None, cls="box"):
    head = [html.Div([title, hint(tip)], className="box-title")] if title else []
    return html.Div(head + list(children), className=cls)


def table(id, page_size=15, tips=None, **kw):
    """tips: {column id: header tooltip}."""
    return dash_table.DataTable(
        id=id, page_size=page_size, tooltip_header=tips or {}, tooltip_delay=200, tooltip_duration=None, sort_action="native", filter_action="none", cell_selectable=True,
        style_as_list_view=True,
        style_table={"overflowX": "auto"},
        style_cell={"fontFamily": "Menlo, monospace", "fontSize": 11, "padding": "2px 8px", "textAlign": "right",
                    "border": "none", "minWidth": 56, "maxWidth": 180, "overflow": "hidden", "textOverflow": "ellipsis"},
        style_header={"fontFamily": "inherit", "fontWeight": 600, "borderBottom": "1px solid #ccc", "background": "white"},
        style_cell_conditional=[{"if": {"column_id": c}, "textAlign": "left"} for c in ("seq", "rc", "kmer", "nbr", "core", "model", "method")],
        style_data_conditional=[{"if": {"state": "active"}, "backgroundColor": "#eef3fa", "border": "none"}],
        **kw)

#!/usr/bin/env python3
"""Nucleomer dashboard.  Run:  ./run.sh   (or  python app.py [--port 8050])  then open http://127.0.0.1:8050"""
import argparse
import json
import sys
from pathlib import Path

from nucleomer_dash import store as S

ROOT = Path(__file__).resolve().parent
if not (S.DATA / "manifest.json").exists():
    sys.exit(f"No data in {S.DATA}. Run ./setup.sh (or: python get_data.py).")
DATA_VERSION = S.manifest().get("data_version", "?")
WANT = json.loads((ROOT / "data.json").read_text()).get("version")
if WANT and WANT != DATA_VERSION:
    print(f"Note: data {DATA_VERSION} is installed but this version of the code expects {WANT}. "
          f"Run: python get_data.py", file=sys.stderr)

from dash import Dash, Input, Output, dcc, html  # noqa: E402

from nucleomer_dash import __version__  # noqa: E402
from nucleomer_dash import theme  # noqa: E402,F401  (registers the plot template)
from nucleomer_dash.pages import compare, embed, explore, graph, help, kmer, seg, stats  # noqa: E402
from nucleomer_dash.ui import hint  # noqa: E402

PAGES = {"K-mers": explore, "K-mer": kmer, "Graph": graph, "Map": embed, "Compare": compare, "Stats": stats,
         "Segmentation": seg, "Help": help}
P = dict(persistence=True, persistence_type="session")

app = Dash(__name__, title="Nucleomer", suppress_callback_exceptions=True, update_title=None)
ms = S.models()
app.layout = html.Div([
    html.Div([
        html.Div([html.Span("Nucleomer", className="brand"),
                  html.Span(f"v{__version__}", className="ver", title=f"dashboard v{__version__}, data {DATA_VERSION}")],
                 className="brand-box"),
        dcc.Tabs(id="tab", value="K-mers", className="tabs", **P,
                 children=[dcc.Tab(label=t, value=t, className="tab", selected_className="tab--sel") for t in PAGES]),
        html.Div([
            html.Label(["Model", hint("sequence-to-activity model whose Nucleomer results are shown")]),
            dcc.Dropdown(id="model", options=ms, value="GATA2" if "GATA2" in ms else ms[0], clearable=False,
                         className="dd-model", **P),
            html.Label(["Value", hint("prediction: mean model output over backgrounds with the k-mer inserted. "
                                      "lift: prediction minus the median of all k-mers of the same length. "
                                      "score: lift / the model's best lift (1 = best k-mer).")]),
            dcc.Dropdown(id="value", options=[{"label": v, "value": k} for k, v in S.VALUE_LABEL.items() if k != "nlp"],
                         value="score", clearable=False, className="dd-value", **P),
            html.Label(["p <", hint("significance threshold for k-mers (node p) and graph edges (edge p), used in every tab")]),
            dcc.Input(id="pthr", type="number", value=S.PVAL, min=0, max=1, step=0.001, debounce=True, className="num", **P),
            html.Label(["K-mer", hint("the selected k-mer (type one, or click a k-mer in any plot or table)")]),
            dcc.Input(id="kmer", type="text", value="AGATAAG", debounce=True, className="txt kmer-in", spellCheck=False, **P),
        ], className="globals"),
    ], className="top"),
    html.Div(id="content"),
])


@app.callback(Output("content", "children"), Input("tab", "value"))
def render(t):
    return PAGES.get(t, explore).layout()


for p in PAGES.values():
    p.register(app)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=8050)
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--debug", action="store_true")
    a = ap.parse_args()
    app.run(host=a.host, port=a.port, debug=a.debug)

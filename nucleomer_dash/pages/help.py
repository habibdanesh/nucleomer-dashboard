"""Help: renders docs/DASHBOARD.md."""
from pathlib import Path

from dash import dcc, html

DOC = Path(__file__).resolve().parents[2] / "docs" / "DASHBOARD.md"


def layout():
    text = DOC.read_text() if DOC.exists() else "docs/DASHBOARD.md not found."
    return html.Div(dcc.Markdown(text, className="help", link_target="_blank"), style={"padding": "10px"})


def register(app):
    pass

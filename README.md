# Nucleomer dashboard

An interactive browser app for exploring Nucleomer results: all k-mers of 13 models with their marginalisation values
and p-values, the k-mer graphs, maps of top k-mers, model-vs-model comparisons and the segmentation eval.

It runs locally: a small Python server plus your web browser. Nothing is uploaded anywhere.

## Quick start (macOS / Linux)

Requirements: Python 3.10 or newer, about 6 GB of free disk space (5 GB after setup) and 4 GB of RAM.
```bash
git clone https://github.com/habibdanesh/nucleomer-dashboard.git
cd nucleomer-dashboard
./setup.sh      # once: creates .venv, installs packages, downloads the data (~1.5 GB)
./run.sh        # starts the dashboard and opens http://127.0.0.1:8050
```

Stop it with Ctrl+C. Next time only `./run.sh` is needed.

**Update to a new version:**

```bash
git pull
./setup.sh      # installs new packages and downloads new data only if needed
```

## Using it

- The **top bar** applies to every tab: the model, the value to plot (score, lift or prediction), the p-value
  threshold, and the selected k-mer.
- **Click** any k-mer (a point, a table row, a graph node) to select it everywhere; the K-mer tab shows it in detail.
- Hover over a **?** to see what a plot, column or control means. The **Help** tab has the full reference.
- The camera icon above each plot downloads it as SVG.

| Tab | What it shows |
|---|---|
| K-mers | all k-mers of one model: filters, density plots, table |
| K-mer | one k-mer: value in every model, per background, single-base edits, graph neighbours |
| Graph | the Nucleomer graph around the selected k-mer (2D/3D) |
| Map | UMAP / spring layout of the top k-mers of one or more models, with random decoys |
| Compare | model A vs B, correlation of all models, top k-mers of A looked up in B |
| Stats | distributions per model, summary table |
| Segmentation | eval accuracy, PR/ROC curves, chr1 locus browser |
| Help | reference: every view and how it is computed |

## Troubleshooting

- **The download fails or is slow.** Open the data link from `data.json` in a browser, download the `.tar.gz`, and
  run `./setup.sh ~/Downloads/<file>.tar.gz`.
- **"Python 3.10 or newer was not found".** Install Python (`brew install python@3.12`, or from python.org) and run
  `./setup.sh` again.
- **Port 8050 is in use.** `./run.sh --port 8051`.
- **Windows.** Not tested. In a terminal: `python -m venv .venv`, `.venv\Scripts\pip install -r requirements.txt`,
  `.venv\Scripts\python get_data.py`, `.venv\Scripts\python app.py`, then open http://127.0.0.1:8050.

## Repository

```
app.py, nucleomer_dash/   the app (one module per tab in nucleomer_dash/pages/)
setup.sh, run.sh          install and start
get_data.py, data.json    download + checksum of the data bundle
docs/DASHBOARD.md         reference (also the Help tab)
build/                    maintainers: raw Nucleomer results -> data bundle
CHANGELOG.md
```

## Maintainers: releasing a new version

```bash
PY=/path/to/python-with-pysam-and-pyBigWig    # pip install -r requirements.txt -r build/requirements.txt
$PY build/prepare.py --force                  # Nucleomer results (paths in build/sources.json) -> data/models
$PY build/prepare_seg.py                      # eval results -> data/seg, data/genome
$PY build/bundle.py                           # -> dist/nucleomer-dashboard-data-<date>.tar.gz, updates data.json
```

Upload the `.tar.gz` (Google Drive: share as "Anyone with the link"), paste the share link into `data.json` → `url`,
bump `nucleomer_dash/__init__.py`, add a `CHANGELOG.md` entry, commit, tag `vX.Y.Z` and push. Code-only changes
need no new bundle.

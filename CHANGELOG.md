# Changelog

Versions are git tags (`vX.Y.Z`). The data bundle has its own version (a date), listed in `data.json`.

## 0.2.1 (2026-11-18)

- `data.json` now holds the Google Drive link to data bundle 2026.11.18, so `./setup.sh` downloads it automatically (tested from a fresh clone: 1.3 GB in about 40 s, checksum OK).
- README: clone command.

## 0.2.0 (2026-11-18)

- **Standalone repository.** Install with `./setup.sh` (Python venv, packages, data download with checksum); start with `./run.sh`.
- **Self-contained data bundle.** The app reads only `data/`, so the original result files, genome and bigWigs are no longer needed. On disk it went from 6 GB to 3.4 GB; the download is a single `.tar.gz`.
  - p-values are stored as integers (p × 1000, exact).
  - Edges keep only the out-degree, target, p and Δ prediction; the other columns are rebuilt on load.
  - Per-background predictions are kept for k ≤ 6 and the top 25,000 k-mers per model (plus reverse complements and the median k-mer of each length).
  - The bundle includes chr1 sequence and ChIP-seq signal (10 bp bins).
  - The unused prediction SD was dropped.
- **Self-descriptive.** Every plot has a title and a **?** tooltip. Controls and table columns have tooltips, and each tab has a one-line description. Legends now name every colour, density and line. A **Help** tab renders `docs/DASHBOARD.md`.
- **UMAP.** The Map uses the same UMAP layout on umap-learn 0.5 and 0.6+ (`compatibility_layout`), so maps match across installs.
- **Docs.** README is now a quick start; the reference moved to `docs/DASHBOARD.md`; the build scripts moved to `build/`.

## 0.1.0 (2026-11-18)

- First version: 13 models; 7 tabs (K-mers, K-mer, Graph, Map, Compare, Stats, Segmentation); data converters; eval v7.

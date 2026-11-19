# Nucleomer dashboard — reference

What every view shows and how it is computed. (This page is also the **Help** tab of the dashboard.)

## Layout

- **Top bar (applies to every tab).**
  - **Model:** the primary model.
  - **Value:** the per-k-mer quantity that is plotted.
  - **p <:** the significance threshold. The default 0.01 is Nucleomer's own node and edge threshold.
  - **K-mer:** the selected k-mer.
- **Selecting a k-mer.** Clicking a k-mer anywhere (point, table row, graph node, locus label) selects it. Every tab then follows that selection.
- **Remembered settings.** Controls keep their values for the browser session.
- **Exporting figures.** Every plot has a camera button, which downloads an SVG for drafting figures. Plotly's "Share chart" button, which uploads to Plotly Cloud, is removed, so nothing leaves the machine.
- **Help in the app.** Hover over a **?** for what a plot, column or control shows. The grey line at the top of each tab's left panel says what the tab is for.

### Values (top-bar "Value")

| value | definition |
|---|---|
| prediction | mean model output over the backgrounds with the k-mer inserted at the centre (= Nucleomer node `pred`; model units, e.g. log counts) |
| lift | prediction − median prediction of all k-mers of the same length |
| score | lift / max lift over all significant (p < 0.01) k-mers of the model. 1 = the model's best k-mer. It is comparable across models |

- **p-value:** this is Nucleomer's node p-value, a bootstrap over k-mers with the same dinucleotide composition (1,000 samples). Its resolution is therefore 1e-3, so p = 0 is shown as "<0.001" and as −log10 p = 3.
- **Reverse complements:** they are *not* merged in the stored data, because the models are not strand-symmetric. Where a view offers "merge reverse complements", each k-mer/RC pair becomes one item. That item is the canonical (lexicographically smaller) sequence, takes the larger of the two values, and is significant if either strand is.

### Tabs

**K-mers** — every k-mer (k = 1..maxk) of the primary model.
- *Left plot:* density of all k-mers on a value × k grid. Each k row is scaled to its own maximum (log counts), so short k stay visible.
- *Right plot:* density of value × −log10 p.
- *Blue points and the table:* the k-mers passing the filters. Filters are the k range, significant only, regex pattern (matched on either strand), and RC merge. The top N are shown, sorted by value.
- *Table columns:* the `core` column is the k-mer's core from Nucleomer's `find_cores`, which walks the deletion DAG while keeping ≥ 80 % of the prediction.

**K-mer** — one k-mer in depth.
- *Info panel:* sequence and RC, value, p, the core and the core's alignment in the k-mer. TOMTOM hits are shown if the k-mer is one of the model's top cores.
- *Across models:* the value of the k-mer (dark) and its RC (light) in every model.
  - Bars are coloured by model family (BPNet blue, Cherimoya orange, ProCapNet green).
  - A black outline means p < threshold.
- *Per background:*
  - Each point is one background sequence. The value plotted is (prediction with the k-mer inserted) − (prediction of the same background without it).
  - Groups: the k-mer, its RC, and the median k-mer of the same length as a reference.
  - Per-background predictions are stored for all k ≤ 6, the model's top 25,000 significant k-mers by lift and their reverse complements, and the median k-mer of each length. For other k-mers the panel says so.
- *Single-base edits:* the change in value when one base is substituted, deleted or inserted. Every edit is looked up in the full marginalisation, which is in effect an in-silico mutagenesis of the k-mer itself.
  - Substitutions are shown as a 4 × k grid.
  - Deletions are shown as 1 × k.
  - Insertions are shown as 4 × (k+1), with position = "before base i".
- *Graph neighbours:* every Nucleomer edge touching the k-mer (out = from this k-mer, in = into it). Columns are the edge p, the change in value, and the neighbour's node p. Clicking a row moves to that neighbour.

**Graph** — the Nucleomer graph around the selected k-mer.
- *How nodes are collected:* breadth-first search up to *radius* edges, using only edges with edge p < threshold (or all candidate edges if unticked). Each extra hop is a single-base insertion or deletion.
- *Cap:* when the next ring would exceed *max nodes*, the strongest k-mers (by score) are kept.
- *Layouts:*
  - **k levels** puts one row per k-mer length. Order within a row comes from 4 barycentre sweeps (each node is moved to the mean x of its neighbours), which reduces crossings.
  - **spring** is a networkx spring layout.
- *Display:* 2D or 3D. Colour can be value, −log10 p, k, or the value in another model. Size is the score. The selected k-mer has a red outline.

**Map** — an embedding of the top k-mers of one or more models.
- *Selection:* per model and per k in range, take the top N significant k-mers by score, with RC merged if ticked. The final set is the union over models.
- *Random decoys:* uniformly sampled k-mers of the same lengths, in the same length proportions as the selected set. This gives a background spectrum, the same idea as random background k-mers in KMAP.
- *Distance:* shortest-path length (in single-base indels) inside the selected set. Unconnected pairs get max + 1.
  - **graph edges** uses Nucleomer edges with p < threshold in *any* selected model.
  - **any 1-indel** connects every pair that differs by one insertion or deletion, whether or not it is significant, so decoys join the map. It is a pure sequence distance.
- *Embedding:* UMAP (`metric="precomputed"`; neighbours and min dist adjustable; seed 0) or a spring layout. Either can be 2D or 3D.
- *Colouring options:*
  - best model (argmax value; decoys grey);
  - **value, one panel per model** (same layout, shared colour scale);
  - # models significant;
  - k;
  - regex pattern;
  - decoy flag;
  - value in any model.
- *Display:* size is the max score over the selected models. The selected k-mer is circled.
- *Run button:* the layout is recomputed only on **Run**. Colouring changes are instant.

**Compare** — models against each other at one k (RC merged by default). Every value comes from the *complete* marginalisation (all 4^k k-mers), not the pruned graph, so a k-mer that is absent from model B's graph still has a value and a p in B.
- *A vs B:* every k-mer, coloured by where it is significant (A only, B only, both). The title gives the Spearman correlation over all k-mers. For k ≥ 8 the non-significant bulk is drawn as a density.
- *All models:* Spearman correlation of the per-k-mer value for every pair of models. Models are ordered by average-linkage clustering on 1 − r.
- *Top k-mers of A, looked up in B:* A's top significant k-mers, each with its value, p and percentile in B. For example, a TF model's top k-mers looked up in the accessibility (ATAC) model.
- *Top k-mers × models:* the union of each model's top N significant k-mers, as a heatmap of value per model. Rows are clustered (average linkage, cosine). Clicking a cell selects that k-mer.

**Stats** — the primary model, plus one table over all models.
- *Panels:*
  - fraction of k-mers with p < threshold, per k;
  - node p-value distribution per k (k ≥ 4; darker = longer);
  - value distribution per k (box = quartiles, whiskers = 0.1–99.9 %);
  - edge p-value distribution, split into insertion and deletion edges;
  - background predictions without insertion (`pred.before`).
- *Table:* maxk, number of backgrounds, number of significant k-mers, graph size, and top cores for every model.

**Segmentation** — eval v7 (chr1, top 1M hits).
- *Bars:* AUPR, AUROC (genome), or AUPR / Model AUPR for the ticked methods. The dashed line is random: π = `peak_bp_frac` for AUPR (also divided by the Model AUPR in the relative view), and `auroc_random` for AUROC. The solid line at 1 is the Model.
- *Curves:* PR and genome-wide ROC for one TF. Curves are thinned to ≤ 3,000 points, evenly spaced along the curve.
- *Locus browser:* chr1 only, windows up to 200 kb. Navigate by typing a locus, with the ‹ − + › buttons, or with "peak" (a random ChIP-seq peak). The default is the strongest peak of the TF. Tracks:
  - ChIP-seq signal from the ENCODE bigWig, where one exists (GATA2, ELF4, JUNB, KLF1, SP1), stored as 10 bp bin means;
  - peaks;
  - each method's scored segments, exactly as used in the eval;
  - **k-mers** (windows ≤ 5 kb): a live annotation computed from the Nucleomer results of the TF's model. Every substring of length kmin..kmax with p < threshold is a candidate, and non-overlapping k-mers are taken greedily by score. Labels appear at ≤ 250 bp; click one to select it;
  - the sequence, at ≤ 300 bp.

## Data

The app reads only `data/` (the data bundle, downloaded by `./setup.sh`). The bundle is made from the raw Nucleomer
results by `build/prepare.py` and `build/prepare_seg.py`; each script's docstring gives the exact format of every file.

**Indexing.** Every k-mer of length 1..11 has a global id, `gid = OFF[k] + base-4 index` (A=0, C=1, G=2, T=3, first base
most significant). This is the row order of Nucleomer's `kmers_k*.fa` and `pred.after.k*.npy`. All per-model arrays
have the same N = 5,592,404 rows (`nucleomer_dash/kmers.py`), so no sequences are stored.

`data/models/<model>/`:

| file | content |
|---|---|
| `pred.npy` | mean over backgrounds of `marginalization/pred.after.k*.npy` (float32; NaN for k > maxk) |
| `pval.npy` | node p-value × 1000 (`graph/node_pvals_k*.npy`; exact, since p comes from 1000 bootstrap samples) |
| `core.npy` | gid of the core from `cores/all_cores.maxk*.json` (−1 if not in the graph) |
| `edges.npz` | every graph edge (`graph_edges.maxk*.csv`): out-degree per k-mer, target, edge p × 1000, Δ prediction. Insertion or deletion follows from the lengths |
| `bg.npz` | per-background predictions for a subset of k-mers (see *Per background*) |
| `meta.json` | family, assay, maxk, number of backgrounds, `pred.before`, median k-mer per k, graph size, top cores and their TOMTOM hits |
| `tomtom/` | TOMTOM outputs of the top cores |

`data/seg/` (eval v7, `08-expanded-eval-of-v5/results-chr1-v7`) and `data/genome/`:

| file | content |
|---|---|
| `seg/summary.tsv` | the eval summary |
| `seg/<TF>/<method>.npz` | scored segments on chr1 (best run per method) |
| `seg/<TF>/<method>.curves.npz` | PR and ROC, thinned to ≤ 3,000 points |
| `seg/<TF>/peaks.npz` | chr1 ChIP-seq peaks |
| `seg/<TF>/chip.npy` | ChIP-seq bigWig on chr1, 10 bp bin means |
| `genome/chr1.npy` | hg38 chr1 sequence |

**Models (13):**
- BPNet: GATA2, JUNB, ELF4, SOX6, MYC, CTCF, ATAC.
- Cherimoya (old checkpoints, as in eval v7): KLF1, TEAD4, ZBTB33, NRF1, SP1.
- ProCapNet: PRO-cap.
- maxk is 11 except ATAC (10) and PRO-cap (8): only those k were marginalised.
- Malinois is not included: its results folder has only k-mer FASTAs and backgrounds, no predictions.

**Checks.** Converted values match Nucleomer's own `graph_nodes.csv`: prediction to ≤ 0.001 (fp16 accumulation in the
original), p to 3e-8, edge p and Δ prediction to the CSV's 4 decimals (SOX6, 5,000 random nodes and edges). Edge types in the CSV match the k-mer lengths for every edge (asserted while building).

**Memory.** Node arrays are memory-mapped (~60 MB per model); a model's edges take 100–500 MB once loaded, and up to 4
models are kept. A typical session uses 1–3 GB of RAM.

## Related work and what it suggests

- **Global importance analysis.** GIA ([Koo et al. 2021, PLoS Comput Biol](https://repository.cshl.edu/id/eprint/40066/)) inserts patterns into background sequences to measure their population-level effect. That is the same principle as Nucleomer's marginalisation. GIA was then extended to the number of motifs, their spacing, and context such as GC content.
- **tangermeme.** [tangermeme](https://www.biorxiv.org/content/10.1101/2025.08.08.669296.full.pdf) is a toolkit for "everything but the model". It offers `marginalize`, `ablate` and `space` (motif-pair spacing) experiments, and its [docs](https://tangermeme.readthedocs.io/en/latest/tutorials/Tutorial_B1_Marginalization.html) stress that the before/after difference is the quantity of interest. That is why the *Per background* panel shows differences.
- **KMAP.** [KMAP (Genome Research 2025)](https://genome.cshlp.org/content/35/5/1234.full.pdf) is a UMAP-like 2D embedding built specifically for k-mers. It shows motif clusters against a cloud of random k-mers and notes that plain UMAP on k-mers can collapse. The Map tab's decoys and the "any 1-indel" distance follow that idea.
- **ChromBPNet / TF-MoDISco reports.** These match learned motifs to databases with TOMTOM and use *marginal footprints* to verify unclear motifs ([example report](https://mitra.stanford.edu/kundaje/oak/akundaje/anusri/chrombpnet_data2/chrombpnet_model/evaluation/k562_overall_report.html)). They are static HTML, and none of these tools offers an interactive k-mer-level explorer.

## Roadmap (not yet built)

**Needs new model runs:**
- motif-pair spacing and orientation (GIA / tangermeme `space`), e.g. GATA–E-box (TAL1) at 7–8 bp, which maxk = 11 cannot capture;
- Nucleomer for the new Cherimoya fold models.

**From existing data:**
- k-mer clusters → PWM / logo (align the members of a Map cluster, or a core and its supersets);
- JASPAR best-match annotation for colouring (needs the JASPAR 2024 MEME file locally);
- a "figure mode" that re-draws the current view in matplotlib with the paper style;
- saved/shareable views (URL state).

**Release:** a public host (e.g. a lab server or Hugging Face Space) and an archived data bundle (Zenodo) with the paper.


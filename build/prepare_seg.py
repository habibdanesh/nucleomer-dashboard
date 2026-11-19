#!/usr/bin/env python3
"""Convert one segmentation-eval results folder (08-expanded-eval-of-v5/results-chr1-v*) for the dashboard.

Writes data/seg/:
  summary.tsv                 copy of the eval summary
  meta.json                   TFs, methods (best run of each), source file names, region
  <TF>/<method>.npz           segments on chr1, sorted by start: dstart uint32 (start - previous start),
                              length uint16, score float16
  <TF>/<method>.curves.npz    PR (recall, precision) and genome ROC (fpr, tpr), <= 3000 points each
  <TF>/peaks.npz              ChIP-seq peaks on chr1: start, end (int32), signal (float32)
  <TF>/chip.npy               ChIP-seq bigWig on chr1, mean over 10 bp bins (float16), if a bigWig exists
and data/genome/chr1.npy      chr1 sequence as uint8 ASCII (upper case)

Usage:  python build/prepare_seg.py [--root DIR]
"""
import argparse
import glob
import json
import shutil

import numpy as np
import pandas as pd

from prepare import DATA, load_sources

CHROM = "chr1"
BIN = 10


def thin(df, x, y, n=3000):
    """Keep <= n points, evenly spaced along the curve's arc length (keeps the shape, drops duplicates)."""
    a = df[[x, y]].to_numpy(float)
    a = a[np.isfinite(a).all(1)]
    if len(a) <= n:
        return a
    d = np.r_[0, np.cumsum(np.hypot(*np.diff(a / (np.ptp(a, 0) + 1e-12), axis=0).T))]
    keep = np.unique(np.searchsorted(d, np.linspace(0, d[-1], n)).clip(0, len(a) - 1))
    return a[keep]


def save_segments(path, s):
    s = s.sort_values("start")
    st, ln, sc = s.start.to_numpy(np.int64), (s.end - s.start).to_numpy(np.int64), s.score.to_numpy(np.float64)
    assert ln.min() > 0 and ln.max() < 65536 and np.abs(sc).max() < 65000
    np.savez(path, dstart=np.diff(st, prepend=0).astype(np.uint32), length=ln.astype(np.uint16), score=sc.astype(np.float16))


def chip_bins(bw_path, length):
    import pyBigWig
    bw = pyBigWig.open(bw_path)
    v = np.nan_to_num(bw.values(CHROM, 0, length, numpy=True))
    n = length // BIN
    return v[:n * BIN].reshape(n, BIN).mean(1).astype(np.float16)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root")
    a = ap.parse_args()
    reg, base = load_sources(a.root)
    seg = reg["segmentation"]
    res, cfg_path = base / seg["results"], base / seg["config"]
    cfg = json.loads(cfg_path.read_text())
    out = DATA / "seg"
    out.mkdir(parents=True, exist_ok=True)

    import pysam
    fa = pysam.FastaFile(str((base / seg["genome"]).resolve()))
    L = fa.get_reference_length(CHROM)
    (DATA / "genome").mkdir(exist_ok=True)
    np.save(DATA / "genome" / f"{CHROM}.npy", np.frombuffer(fa.fetch(CHROM).upper().encode(), np.uint8))

    shutil.copy(res / "summary.tsv", out / "summary.tsv")
    summ = pd.read_csv(res / "summary.tsv", sep="\t")
    meta = {"source": str(res.relative_to(base)), "genome": "hg38", "chrom": CHROM, "chrom_length": L, "chip_bin": BIN,
            "tfs": {}}
    for m in cfg["models"]:
        tf = m["name"]
        d = out / tf
        d.mkdir(exist_ok=True)
        bed = (cfg_path.parent / m["chipseq_bed"]).resolve()
        p = pd.read_csv(bed, sep="\t", header=None, usecols=[0, 1, 2, 6])
        p = p[p[0] == CHROM].sort_values(1)
        np.savez(d / "peaks.npz", start=p[1].to_numpy(np.int32), end=p[2].to_numpy(np.int32), signal=p[6].to_numpy(np.float32))
        bw = sorted(glob.glob(str(bed.parent / f"K562.{tf}.*.bigWig")))
        if bw:
            np.save(d / "chip.npy", chip_bins(bw[0], L))
        methods = []
        for mdir in sorted((res / "segments" / tf).iterdir()):
            r = summ[(summ.model == tf) & (summ.method == mdir.name)]
            if r.empty:
                continue
            run = r.sort_values("aupr").iloc[-1]["run"]           # best run if a method was swept
            s = pd.read_csv(mdir / f"{run}.scored.tsv", sep="\t", usecols=["start", "end", "score"])
            save_segments(d / f"{mdir.name}.npz", s)
            sc = res / "scores" / tf / mdir.name / str(run)
            np.savez(d / f"{mdir.name}.curves.npz", pr=thin(pd.read_csv(sc / "pr.tsv", sep="\t"), "recall", "precision"),
                     roc=thin(pd.read_csv(sc / "roc_genome.tsv", sep="\t"), "fpr", "tpr"))
            methods.append({"method": mdir.name, "run": str(run), "n": int(len(s))})
            print(f"{tf} {mdir.name} {run}: {len(s):,} segments", flush=True)
        meta["tfs"][tf] = {"family": m["family"], "peaks": bed.name, "chip": bool(bw), "chip_file": bw[0].split("/")[-1] if bw else None,
                           "methods": methods, "n_peaks": int(len(p))}
    regions = pd.read_csv(cfg_path.parent / cfg["target_bed"], sep=r"\s+", header=None)
    meta["region"] = [str(regions.iloc[0, 0]), int(regions[1].min()), int(regions[2].max())]
    (out / "meta.json").write_text(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()

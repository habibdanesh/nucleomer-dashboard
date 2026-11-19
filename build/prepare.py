#!/usr/bin/env python3
"""Convert Nucleomer results (one folder per model) into the compact arrays the dashboard reads.

For every model in build/sources.json, writes data/models/<name>/:
  pred.npy    float32[N]   mean prediction over backgrounds (= graph node 'pred'); NaN for k > maxk
  pval.npy    uint16[N]    node p-value x 1000 (bootstrap with 1000 samples, so exact); 65535 = not computed
  core.npy    int32[N]     gid of the k-mer's core (cores/all_cores.*.json); -1 if not in the graph
  edges.npz   out_deg uint8[N]  number of graph edges leaving each k-mer (edges are sorted by source, then target)
              dst int32[E]      target k-mer;  insertion if len(target) > len(source), else deletion
              pval uint16[E]    edge p-value x 1000
              dpred float16[E]  pred(target) - pred(source)
  bg.npz      gid int32[n], after float16[n, B]: per-background predictions for a subset of k-mers
              (all k <= 6, the top BG_TOP k-mers by lift and their reverse complements, the median k-mer of each k)
  tomtom/     TOMTOM tables of the top cores
  meta.json   model info, background predictions without insertion, median k-mer per k, top cores, TOMTOM hits
N = 5,592,404 = all k-mers of length 1..11, indexed as in nucleomer_dash/kmers.py.

Usage:  python build/prepare.py [--models GATA2 ATAC ...] [--root DIR] [--force]
"""
import argparse
import json
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from nucleomer_dash.kmers import KMAX, N, OFF, encode, revcomp_gid  # noqa: E402

DATA = Path(os.environ.get("NUCLEOMER_DATA", ROOT / "data"))
P_NA = 65535
BG_TOP = 25_000


def load_sources(root=None):
    reg = json.loads((ROOT / "build" / "sources.json").read_text())
    base = Path(root or os.environ.get("NUCLEOMER_RESULTS") or (ROOT / reg["root"])).resolve()
    for m in reg["models"]:
        m["src"] = base / m["path"]
    return reg, base


def milli(p):
    """p-values that are multiples of 1/1000 -> uint16 (P_NA for NaN)."""
    q = np.round(np.asarray(p, np.float64) * 1000)
    ok = np.isfinite(q)
    assert np.abs(q[ok] - np.asarray(p, np.float64)[ok] * 1000).max(initial=0) < 1e-3, "p-values are not multiples of 1e-3"
    return np.where(ok, q, P_NA).astype(np.uint16)


def detect_maxk(src):
    ks = [k for k in range(1, KMAX + 1)
          if (src / f"marginalization/pred.after.k{k}.npy").exists() and (src / f"graph/node_pvals_k{k}.npy").exists()]
    assert ks == list(range(1, max(ks) + 1)), f"gaps in k: {ks}"
    return max(ks)


def read_tomtom(path, top=3):
    try:
        t = pd.read_csv(path, sep="\t", comment="#")
    except Exception:
        return []
    t.columns = [c.strip() for c in t.columns]
    t = t.dropna(subset=["Target Name"]).head(top)
    return [{"target": str(r["Target Name"]).strip(), "p": float(r["p-value"])} for _, r in t.iterrows()]


def prepare(m, base, force=False):
    out = DATA / "models" / m["name"]
    if (out / "meta.json").exists() and not force:
        print(f"[{m['name']}] exists, skip")
        return
    out.mkdir(parents=True, exist_ok=True)
    src, t0 = m["src"], time.time()
    maxk = detect_maxk(src)

    # nodes
    pred = np.full(N, np.nan, np.float32)
    pval = np.full(N, np.nan, np.float32)
    for k in range(1, maxk + 1):
        P = np.load(src / f"marginalization/pred.after.k{k}.npy", mmap_mode="r")
        n_bg = P.shape[1]
        for s in range(0, P.shape[0], 1 << 18):
            pred[OFF[k] + s:OFF[k] + s + min(1 << 18, P.shape[0] - s)] = np.asarray(P[s:s + (1 << 18)], np.float32).mean(1)
        pval[OFF[k]:OFF[k + 1]] = np.load(src / f"graph/node_pvals_k{k}.npy")
    np.save(out / "pred.npy", pred)
    np.save(out / "pval.npy", milli(pval))
    print(f"[{m['name']}] maxk={maxk} B={n_bg} nodes {time.time() - t0:.0f}s", flush=True)

    # per-background predictions for a subset (the full matrices are ~1 GB per model)
    med = {k: int(OFF[k] + np.argsort(pred[OFF[k]:OFF[k + 1]])[(4 ** k) // 2]) for k in range(1, maxk + 1)}
    lift = pred.copy()
    for k in range(1, maxk + 1):
        lift[OFF[k]:OFF[k + 1]] -= np.median(pred[OFF[k]:OFF[k + 1]])
    sig = np.where(pval < 0.01)[0]
    top = sig[np.argsort(-lift[sig])[:BG_TOP]]
    sel = np.unique(np.r_[np.arange(OFF[min(6, maxk) + 1]), top, revcomp_gid(top), list(med.values())]).astype(np.int64)
    after = np.empty((len(sel), n_bg), np.float16)
    for k in range(1, maxk + 1):
        i = np.where((sel >= OFF[k]) & (sel < OFF[k + 1]))[0]
        if len(i):
            P = np.load(src / f"marginalization/pred.after.k{k}.npy", mmap_mode="r")
            after[i] = P[sel[i] - OFF[k]]
    np.savez(out / "bg.npz", gid=sel.astype(np.int32), after=after)

    # edges
    e = pd.read_csv(src / f"graph/graph_edges.maxk{maxk}.csv", usecols=["source", "target", "pval", "dPred", "edge_type"],
                    dtype={"source": str, "target": str, "pval": np.float32, "dPred": np.float32, "edge_type": str})
    s_, d_ = encode(e.source.to_numpy()), encode(e.target.to_numpy())
    assert np.array_equal(e.edge_type.to_numpy() == "del", d_ < s_), "edge_type does not match k-mer lengths"
    o = np.lexsort((d_, s_))
    deg = np.bincount(s_, minlength=N)
    assert deg.max() < 256
    np.savez(out / "edges.npz", out_deg=deg.astype(np.uint8), dst=d_[o].astype(np.int32),
             pval=milli(e.pval.to_numpy()[o]), dpred=e.dPred.to_numpy()[o].astype(np.float16))
    n_edges, n_nodes = len(e), int(len(np.unique(np.r_[s_, d_])))
    del e, s_, d_
    print(f"[{m['name']}] {n_edges:,} edges {time.time() - t0:.0f}s", flush=True)

    # cores, TOMTOM
    core = np.full(N, -1, np.int32)
    cj = src / f"cores/all_cores.maxk{maxk}.json"
    if cj.exists():
        c = json.loads(cj.read_text())
        core[encode(list(c.keys()))] = encode([v[0] for v in c.values()])
    np.save(out / "core.npy", core)
    tj = sorted((src / "cores").glob(f"top_*_cores.maxk{maxk}.json"), key=lambda p: -int(p.name.split("_")[1]))
    top_cores = json.loads(tj[0].read_text()) if tj else {}
    shutil.rmtree(out / "tomtom", ignore_errors=True)
    (out / "tomtom").mkdir()
    tomtom = {}
    for f in (src / "tomtom").glob("tomtom.*.txt"):
        shutil.copy(f, out / "tomtom" / f.name)
        tomtom[f.name.split(".")[1]] = read_tomtom(f)

    meta = dict(name=m["name"], family=m["family"], assay=m["assay"], maxk=maxk, n_backgrounds=n_bg,
                source=str(src.relative_to(base)) if src.is_relative_to(base) else src.name,
                pred_before=np.load(src / "marginalization/pred.before.npy").astype(float).ravel().tolist(),
                median_kmer={str(k): v for k, v in med.items()}, n_bg_kmers=int(len(sel)),
                n_edges=n_edges, n_graph_nodes=n_nodes,
                top_cores={k: [v[0], float(v[1])] for k, v in top_cores.items()}, tomtom=tomtom,
                prepared=time.strftime("%Y-%m-%d %H:%M"))
    (out / "meta.json").write_text(json.dumps(meta, indent=1))
    print(f"[{m['name']}] done {time.time() - t0:.0f}s", flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--models", nargs="*")
    p.add_argument("--root", help="folder the paths in sources.json are relative to")
    p.add_argument("--force", action="store_true")
    a = p.parse_args()
    reg, base = load_sources(a.root)
    for m in reg["models"]:
        if not a.models or m["name"] in a.models:
            prepare(m, base, a.force)
    (DATA / "models" / "order.json").write_text(json.dumps([m["name"] for m in reg["models"]]))


if __name__ == "__main__":
    main()

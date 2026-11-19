"""Read-only access to the data bundle (data/models, data/seg, data/genome). Arrays are memory-mapped or cached.
The layout of every file is documented in build/prepare.py and build/prepare_seg.py."""
import json
import os
import re
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from .kmers import KMAX, OFF, ALPH, k_of, revcomp_gid, encode

DATA = Path(os.environ.get("NUCLEOMER_DATA", Path(__file__).resolve().parent.parent / "data"))
P_NA = 65535           # missing value in the uint16 p-value arrays (p x 1000)
PVAL = 0.01            # Nucleomer's node and edge p-value threshold
P_FLOOR = 1e-3         # p-values are bootstrap estimates with 1000 samples; p = 0 is shown as < 1e-3


def OFF_K(k):
    """Slice of the global index holding all k-mers of length k."""
    return slice(int(OFF[k]), int(OFF[k + 1]))


def manifest():
    p = DATA / "manifest.json"
    return json.loads(p.read_text()) if p.exists() else {}


@lru_cache(None)
def models():
    """Models in the bundle, in the order of build/sources.json."""
    p = DATA / "models" / "order.json"
    order = json.loads(p.read_text()) if p.exists() else sorted(x.name for x in (DATA / "models").glob("*"))
    return [m for m in order if (DATA / "models" / m / "meta.json").exists()]


@lru_cache(None)
def meta(m):
    return json.loads((DATA / "models" / m / "meta.json").read_text())


def unmilli(q):
    """uint16 p x 1000 -> float32 p (NaN where missing)."""
    q = np.asarray(q)
    return np.where(q == P_NA, np.nan, q / 1000).astype(np.float32)


@lru_cache(None)
def arr(m, name):
    a = np.load(DATA / "models" / m / f"{name}.npy", mmap_mode="r")
    return unmilli(a) if name == "pval" else a


@lru_cache(None)
def value(m, kind="score"):
    """Per-k-mer value used throughout the app.
    pred  = mean prediction with the k-mer inserted (model units)
    lift  = pred - median pred of all k-mers of the same length
    score = lift / max lift over significant k-mers of any length (1 = the model's best k-mer)
    nlp   = -log10 p (p = 0 -> 3)"""
    if kind == "nlp":
        return -np.log10(np.maximum(np.asarray(arr(m, "pval")), P_FLOOR)).astype(np.float32)
    pred = np.asarray(arr(m, "pred"))
    if kind == "pred":
        return pred
    med = np.full(len(pred), np.nan, np.float32)
    for k in range(1, KMAX + 1):
        med[OFF[k]:OFF[k + 1]] = np.nanmedian(pred[OFF[k]:OFF[k + 1]]) if np.isfinite(pred[OFF[k]]) else np.nan
    lift = pred - med
    if kind == "lift":
        return lift
    sig = np.asarray(arr(m, "pval")) < PVAL
    return (lift / np.nanmax(lift[sig])).astype(np.float32)


VALUE_LABEL = {"score": "score", "lift": "lift", "pred": "prediction", "nlp": "-log10 p"}


@lru_cache(4)
def edges(m):
    """Graph edges with CSR indices both ways: out_ptr over src (edges are sorted by src), in_order/in_ptr over dst."""
    with np.load(DATA / "models" / m / "edges.npz") as z:
        deg, dst, pv, dpred = z["out_deg"], z["dst"], z["pval"], z["dpred"]
    n = len(deg)
    out_ptr = np.r_[0, np.cumsum(deg, dtype=np.int64)].astype(np.int32)
    src = np.repeat(np.arange(n, dtype=np.int32), deg)
    in_order = np.argsort(dst, kind="stable").astype(np.int32)
    in_ptr = np.searchsorted(dst[in_order], np.arange(n + 1)).astype(np.int32)
    return dict(src=src, dst=dst, pval=unmilli(pv), dpred=dpred.astype(np.float32), is_del=k_of(dst) < k_of(src),
                out_ptr=out_ptr, in_order=in_order, in_ptr=in_ptr)


def neighbours(m, g):
    """All graph edges touching k-mer g (both directions) as a DataFrame."""
    E = edges(m)
    o = slice(E["out_ptr"][g], E["out_ptr"][g + 1])
    i = E["in_order"][E["in_ptr"][g]:E["in_ptr"][g + 1]]
    df = pd.concat([
        pd.DataFrame({"dir": "out", "nbr": E["dst"][o], "pval": E["pval"][o], "dpred": E["dpred"][o], "is_del": E["is_del"][o]}),
        pd.DataFrame({"dir": "in", "nbr": E["src"][i], "pval": E["pval"][i], "dpred": E["dpred"][i], "is_del": E["is_del"][i]}),
    ], ignore_index=True)
    return df


@lru_cache(4)
def _bg(m):
    with np.load(DATA / "models" / m / "bg.npz") as z:
        return z["gid"], z["after"]


def bg_preds(m, g):
    """Per-background predictions with k-mer g inserted (None if g is not in the bundled subset), and without."""
    before = np.asarray(meta(m)["pred_before"], np.float32)
    gid, after = _bg(m)
    i = np.searchsorted(gid, g)
    return (after[i].astype(np.float32) if i < len(gid) and gid[i] == g else None), before


def median_kmer(m, k):
    """Global id of the k-mer with the median prediction among all k-mers of length k."""
    return int(meta(m)["median_kmer"][str(k)])


def decode_np(gids):
    """Global ids -> numpy unicode array (vectorised per length)."""
    g = np.asarray(gids, np.int64)
    out = np.empty(len(g), dtype=f"U{KMAX}")
    k = k_of(g)
    for kk in np.unique(k):
        sel = np.where(k == kk)[0]
        idx = g[sel] - OFF[kk]
        dig = np.empty((len(sel), kk), np.uint8)
        for j in range(kk - 1, -1, -1):
            dig[:, j] = ALPH[idx % 4]
            idx = idx // 4
        out[sel] = dig.view(f"S{kk}").ravel().astype(f"U{KMAX}")
    return out


def canonical(gids):
    g = np.asarray(gids, np.int64)
    return np.minimum(g, revcomp_gid(g))


def search(gids, pattern):
    """Subset of gids whose sequence (either strand) matches a regex / IUPAC-free pattern."""
    if not pattern:
        return gids
    pat = re.compile(pattern.upper())
    s = pd.Series(decode_np(gids))
    r = pd.Series(decode_np(revcomp_gid(gids)))
    hit = s.str.contains(pat) | r.str.contains(pat)
    return np.asarray(gids)[hit.to_numpy()]


def parse_kmer(s):
    s = (s or "").strip().upper()
    if not s or len(s) > KMAX or set(s) - set("ACGT"):
        return None
    return encode(s)


# ---------- segmentation ----------
@lru_cache(None)
def seg_meta():
    p = DATA / "seg" / "meta.json"
    return json.loads(p.read_text()) if p.exists() else {"tfs": {}}


@lru_cache(None)
def seg_summary():
    p = DATA / "seg" / "summary.tsv"
    return pd.read_csv(p, sep="\t") if p.exists() else pd.DataFrame()


@lru_cache(24)
def seg_track(tf, name):
    """Segments (or peaks) of one track as start, end (int64), score/signal."""
    with np.load(DATA / "seg" / tf / f"{name}.npz") as z:
        t = {k: z[k] for k in z.files}
    if "dstart" in t:
        st = np.cumsum(t.pop("dstart"), dtype=np.int64)
        t = {"start": st, "end": st + t.pop("length"), "score": t.pop("score").astype(np.float32)}
    return t


def seg_window(tf, name, start, end):
    t = seg_track(tf, name)
    # segments are sorted by start and short (<= ~1 kb), so a 1 kb look-back is enough for overlaps
    i0 = np.searchsorted(t["start"], start - 1000)
    i1 = np.searchsorted(t["start"], end)
    keep = slice(i0, i1)
    s, e = t["start"][keep], t["end"][keep]
    m = e > start
    return {k: v[keep][m] for k, v in t.items()}


@lru_cache(None)
def genome():
    return np.load(DATA / "genome" / f"{seg_meta()['chrom']}.npy", mmap_mode="r")


def fetch(chrom, start, end):
    """Upper-case sequence of [start, end) on the bundled chromosome."""
    assert chrom == seg_meta()["chrom"], f"only {seg_meta()['chrom']} is bundled"
    g = genome()
    return bytes(g[max(0, start):min(end, len(g))]).decode()


@lru_cache(16)
def _chip(tf):
    p = DATA / "seg" / tf / "chip.npy"
    return np.load(p, mmap_mode="r") if p.exists() else None


def chip(tf, start, end, nbins):
    """Mean ChIP-seq signal in nbins equal bins over [start, end) (from the 10 bp bins in the bundle)."""
    c, b = _chip(tf), seg_meta()["chip_bin"]
    if c is None:
        return None
    edges_ = np.linspace(start, end, nbins + 1)
    i = np.clip((edges_ // b).astype(np.int64), 0, len(c))
    cs = np.r_[0, np.cumsum(np.asarray(c[i[0]:i[-1] + 1], np.float64))]
    j = np.maximum(i[1:], i[:-1] + 1) - i[0]
    return (cs[np.minimum(j, len(cs) - 1)] - cs[i[:-1] - i[0]]) / (j - (i[:-1] - i[0]))


def fmt_p(p):
    """p-value as text; Nucleomer's bootstrap (1000 samples) cannot resolve p below 1e-3."""
    p = float(p)
    return "" if not np.isfinite(p) else ("<0.001" if p < P_FLOOR else f"{p:.3g}")

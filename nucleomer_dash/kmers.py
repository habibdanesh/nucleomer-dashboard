"""K-mer <-> integer codec shared by prepare.py and the app.

Every k-mer of length 1..KMAX has a global id  gid = OFF[k] + idx,  where idx is the base-4 number
of the k-mer (A=0, C=1, G=2, T=3, first base most significant). This is the row order of
Nucleomer's kmers_k{k}.fa / pred.after.k{k}.npy (itertools.product over "ACGT"), so per-model arrays
are plain concatenations over k and no sequence strings need to be stored.
"""
import numpy as np

KMAX = 11
OFF = np.concatenate([[0, 0], np.cumsum([4 ** k for k in range(1, KMAX + 1)])]).astype(np.int64)
N = int(OFF[KMAX + 1])                      # 5,592,404 k-mers of length 1..11
LUT = np.full(256, 255, np.uint8)
for i, c in enumerate(b"ACGT"):
    LUT[c] = i
ALPH = np.frombuffer(b"ACGT", np.uint8)
COMP = str.maketrans("ACGT", "TGCA")


def k_of(gid):
    """Length of each global id."""
    return np.searchsorted(OFF, np.asarray(gid), side="right") - 1


def encode(seqs):
    """Sequences (str or iterable of str, upper-case ACGT, len <= KMAX) -> global ids (int64)."""
    scalar = isinstance(seqs, str)
    a = np.asarray([seqs] if scalar else list(seqs), dtype=f"S{KMAX}")
    if a.size == 0:
        return np.zeros(0, np.int64)
    L = np.char.str_len(a)
    code = LUT[a.view(np.uint8).reshape(len(a), KMAX)].astype(np.int64)
    if (code[np.arange(KMAX)[None, :] < L[:, None]] == 255).any():
        raise ValueError("non-ACGT character")
    idx = np.zeros(len(a), np.int64)
    for j in range(KMAX):
        m = j < L
        idx[m] = idx[m] * 4 + code[m, j]
    out = OFF[L] + idx
    return int(out[0]) if scalar else out


def decode(gids):
    """Global ids -> list of sequences."""
    g = np.atleast_1d(np.asarray(gids, np.int64))
    k = k_of(g)
    idx = g - OFF[k]
    out = np.empty((len(g), KMAX), np.uint8)
    for j in range(KMAX - 1, -1, -1):
        out[:, j] = ALPH[idx % 4]
        idx //= 4
    # digits were written right-aligned; take the last k characters
    return [bytes(r[KMAX - kk:]).decode() for r, kk in zip(out, k)]


def revcomp_gid(gids):
    """Global id of the reverse complement (vectorised)."""
    g = np.atleast_1d(np.asarray(gids, np.int64))
    k = k_of(g)
    idx = g - OFF[k]
    rc = np.zeros_like(idx)
    for _ in range(KMAX):
        m = k > 0
        rc[m] = rc[m] * 4 + (3 - idx[m] % 4)
        idx[m] //= 4
        k = k - m
    return OFF[k_of(g)] + rc


def rc(s):
    return s.translate(COMP)[::-1]


def gids_of_k(k):
    return np.arange(OFF[k], OFF[k + 1], dtype=np.int64)

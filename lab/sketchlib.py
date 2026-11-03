"""Scaled MinHash sketching (same as build_index.py) for query sequences."""
import numpy as np

K, SCALE = 21, 50
THRESH = np.uint64((1 << 64) // SCALE)
LUT = np.full(256, 4, np.uint8)
for v, ch in enumerate("ACGT"):
    LUT[ord(ch)] = v
    LUT[ord(ch.lower())] = v
LUT[ord("U")] = 3
LUT[ord("u")] = 3
C1, C2, C3 = np.uint64(0x9E3779B97F4A7C15), np.uint64(0xBF58476D1CE4E5B9), np.uint64(0x94D049BB133111EB)


def splitmix(x):
    with np.errstate(over="ignore"):
        z = x + C1
        z = (z ^ (z >> np.uint64(30))) * C2
        z = (z ^ (z >> np.uint64(27))) * C3
        return z ^ (z >> np.uint64(31))


def sketch_codes(a):
    L = len(a)
    if L < K:
        return np.empty(0, np.uint32)
    n = L - K + 1
    cs = np.concatenate(([0], np.cumsum((a > 3).astype(np.int32))))
    ok = (cs[K:] - cs[:-K]) == 0
    a64 = a.astype(np.uint64)
    fwd = np.zeros(n, np.uint64)
    rc = np.zeros(n, np.uint64)
    two = np.uint64(2)
    three = np.uint64(3)
    for j in range(K):
        seg = a64[j:j + n]
        fwd = (fwd << two) | (seg & three)
        rc |= ((three - (seg & three)) << np.uint64(2 * j))
    canon = np.minimum(fwd, rc)[ok]
    h = splitmix(canon)
    sel = h[h < THRESH]
    return np.unique((sel >> np.uint64(32)).astype(np.uint32))


def sketch(seq):
    return sketch_codes(LUT[np.frombuffer(seq, np.uint8)])



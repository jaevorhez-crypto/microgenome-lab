"""Evolution: distance trees from the stored sketches. Mash-style distance D = -ln(2j / (1 + j)) / k from the Jaccard index j of two scaled sketches, then neighbour-joining.
Reliable for close relatives (same genus / family, roughly above 80-85 % identity); distant sequences saturate at distance 1.0. 16S sketches hold only about 30 hashes, so 16S trees are coarse."""
import json
import math
from pathlib import Path

import numpy as np

from gpu import APP
DB = APP / "data"
K = 21
_names = _lineage = _off = _h = None


def load():
    global _names, _lineage, _off, _h
    if _names is None:
        _names = json.load((DB / "index" / "names.json").open())
        _lineage = json.load((DB / "index" / "lineage.json").open())
        _off = np.fromfile(DB / "index" / "offsets.bin", dtype="<u4").astype(np.int64)
        _h = np.memmap(DB / "index" / "hashes.bin", dtype="<u4", mode="r")


def sketch_of(i):
    load()
    return np.asarray(_h[_off[i]:_off[i + 1]])


def members(term, kind="", limit=30, minhash=40):
    """indexed sequences whose lineage or title contains `term`, thinned evenly to `limit`"""
    load()
    t = term.strip().lower()
    out = []
    for i, n in enumerate(_names):
        if kind and n["kind"] != kind:
            continue
        lin = _lineage.get(n["id"], {}).get("lineage", [])
        if (t and (any(t in x.lower() for x in lin) or t in n["title"].lower())) and n["hashes"] >= minhash:
            out.append(i)
    if len(out) > limit:
        step = len(out) / limit
        out = [out[int(j * step)] for j in range(limit)]
    return out


def label(i):
    load()
    n = _names[i]
    lin = _lineage.get(n["id"], {})
    nm = lin.get("name") or " ".join(n["title"].split()[:3])
    return "%s (%s)" % (nm, n["id"])


def distance(a, b):
    inter = len(np.intersect1d(a, b, assume_unique=True))
    union = len(a) + len(b) - inter
    j = inter / union if union else 0.0
    if j <= 0:
        return 1.0
    return min(1.0, -math.log(2 * j / (1 + j)) / K)


def matrix(sketches):
    n = len(sketches)
    D = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            D[i][j] = D[j][i] = distance(sketches[i], sketches[j])
    return D


def nj(D, names):
    """neighbour-joining; returns a nested tree {name, children:[{..., length}]} (unrooted tree drawn from its last join)"""
    n = len(D)
    if n < 2:
        return {"name": names[0], "children": []}
    nodes = [{"name": nm, "children": []} for nm in names]
    D = [row[:] for row in D]
    idx = list(range(n))
    while len(idx) > 3:
        m = len(idx)
        r = {i: sum(D[i][j] for j in idx) for i in idx}
        best, bi, bj = None, None, None
        for ai in range(m):
            for bj_ in range(ai + 1, m):
                i, j = idx[ai], idx[bj_]
                q = (m - 2) * D[i][j] - r[i] - r[j]
                if best is None or q < best:
                    best, bi, bj = q, i, j
        dij = D[bi][bj]
        li = max(0.0, 0.5 * dij + (r[bi] - r[bj]) / (2 * (m - 2)))
        lj = max(0.0, dij - li)
        new = {"name": "", "children": [dict(nodes[bi], length=li), dict(nodes[bj], length=lj)]}
        u = len(D)
        for row in D:
            row.append(0.0)
        D.append([0.0] * (u + 1))
        nodes.append(new)
        for k in idx:
            if k not in (bi, bj):
                d = 0.5 * (D[bi][k] + D[bj][k] - dij)
                D[u][k] = D[k][u] = d
        idx = [k for k in idx if k not in (bi, bj)] + [u]
    if len(idx) == 3:
        a, b, c = idx
        la = max(0.0, 0.5 * (D[a][b] + D[a][c] - D[b][c]))
        lb = max(0.0, 0.5 * (D[a][b] + D[b][c] - D[a][c]))
        lc = max(0.0, 0.5 * (D[a][c] + D[b][c] - D[a][b]))
        return {"name": "", "children": [dict(nodes[a], length=la), dict(nodes[b], length=lb), dict(nodes[c], length=lc)]}
    a, b = idx
    return {"name": "", "children": [dict(nodes[a], length=D[a][b] / 2), dict(nodes[b], length=D[a][b] / 2)]}


def tree_for(indices, extra=None):
    """tree of the given indexed sequences (+ optionally one extra query sketch with a label)"""
    sks = [sketch_of(i) for i in indices]
    names = [label(i) for i in indices]
    if extra is not None:
        sks.append(extra[1])
        names.append(extra[0])
    D = matrix(sks)
    return nj(D, names), D, names

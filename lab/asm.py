"""Genome assembly from reads: a de Bruijn graph of canonical k-mers, low-coverage k-mers (sequencing errors) dropped, then unitigs (maximal unbranched paths).
Includes a read simulator so the result can be checked against a known genome. Pure Python, meant for small genomes (viruses, plasmids; up to a few hundred kb of reads per second)."""
import random

COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


def simulate_reads(genome, read_len=150, coverage=30, error=0.003, seed=1):
    rng = random.Random(seed)
    L = len(genome)
    n = int(L * coverage / read_len)
    reads = []
    for _ in range(n):
        p = rng.randrange(0, max(1, L - read_len + 1))
        r = genome[p:p + read_len]
        if rng.random() < 0.5:
            r = rc(r)
        if error:
            r = "".join((rng.choice([x for x in "ACGT" if x != c]) if rng.random() < error else c) for c in r)
        reads.append(r)
    return reads


def parse_reads(text, limit=2_000_000):
    """FASTA or FASTQ text -> list of upper-case ACGT strings (reads with other letters are split at them)"""
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    reads = []
    if lines and lines[0].startswith("@"):
        reads = [lines[i] for i in range(1, len(lines), 4)]
    else:
        cur = []
        for l in lines:
            if l.startswith(">"):
                if cur:
                    reads.append("".join(cur))
                cur = []
            else:
                cur.append(l)
        if cur:
            reads.append("".join(cur))
    out = []
    total = 0
    for r in reads:
        r = r.upper().replace("U", "T")
        for piece in "".join(c if c in "ACGT" else " " for c in r).split():
            out.append(piece)
            total += len(piece)
        if total > limit:
            break
    return out


def assemble(reads, k=31, min_cov=None):
    # count canonical k-mers
    cnt = {}
    for r in reads:
        for i in range(len(r) - k + 1):
            km = r[i:i + k]
            c = rc(km)
            if c < km:
                km = c
            cnt[km] = cnt.get(km, 0) + 1
    total_kmers = len(cnt)
    if min_cov is None:
        good = sorted(v for v in cnt.values() if v >= 3)          # error k-mers (count 1-2) are numerous: estimate the real coverage from the well-supported ones
        med = good[len(good) // 2] if good else 1
        min_cov = max(2, int(med * 0.15)) if med >= 8 else 1
    S = set()
    for km, c in cnt.items():
        if c >= min_cov:
            S.add(km)
            S.add(rc(km))
    kept = len(S) // 2
    seen = set()

    def nxt(km):
        t = km[1:]
        return [t + b for b in "ACGT" if (t + b) in S]

    def prv(km):
        t = km[:-1]
        return [b + t for b in "ACGT" if (b + t) in S]

    contigs = []
    for km in S:
        if km in seen:
            continue
        p = prv(km)
        # a unitig starts where the path is not a simple continuation
        if len(p) == 1 and len(nxt(p[0])) == 1:
            continue
        path = km
        seen.add(km)
        seen.add(rc(km))
        cur = km
        while True:
            n = nxt(cur)
            if len(n) != 1:
                break
            m = n[0]
            if len(prv(m)) != 1 or m in seen:
                break
            path += m[-1]
            seen.add(m)
            seen.add(rc(m))
            cur = m
        contigs.append(path)
    # circular genomes leave a cycle with no start: pick up anything unvisited
    for km in S:
        if km in seen:
            continue
        path, cur = km, km
        seen.add(km)
        seen.add(rc(km))
        while True:
            n = nxt(cur)
            if len(n) != 1 or n[0] in seen:
                break
            cur = n[0]
            path += cur[-1]
            seen.add(cur)
            seen.add(rc(cur))
        contigs.append(path)
    # drop contigs that are the reverse complement of a longer one, and tiny tips
    contigs.sort(key=len, reverse=True)
    out, have = [], set()
    for c in contigs:
        key = min(c, rc(c))
        if key in have or len(c) < k + 20:
            continue
        have.add(key)
        out.append(c)
    return out, {"reads": len(reads), "distinct_kmers": total_kmers, "kept_kmers": kept, "min_coverage": min_cov, "k": k}


def stats(contigs):
    lens = sorted((len(c) for c in contigs), reverse=True)
    total = sum(lens)
    acc, n50 = 0, 0
    for l in lens:
        acc += l
        if acc >= total / 2:
            n50 = l
            break
    return {"contigs": len(lens), "total": total, "longest": lens[0] if lens else 0, "n50": n50}


def check_against(contigs, genome, k=31):
    """share of the reference's k-mers present in the contigs, and whether the longest contigs occur in the reference (either strand)"""
    ref = set()
    for s in (genome, rc(genome)):
        for i in range(len(s) - k + 1):
            ref.add(s[i:i + k])
    have = set()
    for c in contigs:
        for s in (c, rc(c)):
            for i in range(len(s) - k + 1):
                have.add(s[i:i + k])
    cov = len(ref & have) / max(1, len(ref))
    g2 = genome + genome[:200]                      # circular genomes: a contig can span the origin
    exact = [(len(c), (c in g2) or (rc(c) in g2)) for c in contigs[:5]]
    return {"reference_kmers_covered": round(cov, 4), "longest_contigs_exact": exact}

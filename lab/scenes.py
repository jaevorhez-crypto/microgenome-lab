"""Scene preparation for the GPU renderers: reads real data (PDB files, FASTA files), works out geometry and the camera, and writes the small
binary input files the renderer programs read (params.bin plus data files). Pure Python, no GPU needed.
Every renderer takes the same camera block: 24 float32 values, 0-2 camera position, 3-5 forward, 6-8 right, 9-11 up, 12 element count, 15 rise."""
import colorsys
import math
import random
import struct
from pathlib import Path

# ---------------------------------------------------------------------------------------------- camera helpers


def _norm(a):
    l = math.sqrt(sum(x * x for x in a))
    return tuple(x / l for x in a)


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def camera_block(cam, fw, rt, cu, count=0, rise=0.0, width=0, height=0):
    p = list(cam) + list(fw) + list(rt) + list(cu) + [float(count), float(width), float(height), float(rise)] + [0.0] * 8
    return struct.pack("<24f", *p)


# ---------------------------------------------------------------------------------------------- molecules
RAD = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "FE": 1.50, "H": 1.20, "P": 1.80, "ZN": 1.40, "MG": 1.45, "CA": 1.80, "NA": 1.70}
CPK = {"N": (0.25, 0.35, 0.95), "O": (0.95, 0.20, 0.15), "S": (0.96, 0.85, 0.25), "FE": (0.92, 0.45, 0.12), "P": (0.95, 0.60, 0.15), "H": (0.9, 0.9, 0.9)}
CHAINS = [(0.45, 0.75, 0.95), (0.96, 0.66, 0.36), (0.62, 0.86, 0.55), (0.86, 0.56, 0.86), (0.95, 0.85, 0.45), (0.55, 0.80, 0.80)]


def _ribbon_spheres(text, cx, cy, cz):
    """Catmull-Rom spline through the C-alpha trace of every chain; radius by secondary structure from the HELIX / SHEET records."""
    lines = text.splitlines()
    sec = {}
    for ln in lines:
        if ln.startswith("HELIX "):
            ch, a, b = ln[19], int(ln[21:25]), int(ln[33:37])
            for r in range(a, b + 1):
                sec[(ch, r)] = "H"
        elif ln.startswith("SHEET "):
            ch, a, b = ln[21], int(ln[22:26]), int(ln[33:37])
            for r in range(a, b + 1):
                sec.setdefault((ch, r), "E")
    chs = {}
    for ln in lines:
        if ln.startswith("ATOM") and ln[12:16] == " CA " and ln[16] in " A":
            chs.setdefault(ln[21], []).append((int(ln[22:26]), float(ln[30:38]) - cx, float(ln[38:46]) - cy, float(ln[46:54]) - cz))
    out_pos, out_col = [], []
    names = sorted(chs)
    for ci, ch in enumerate(names):
        pts = chs[ch]
        runs, cur = [], [pts[0]]
        for a_, b_ in zip(pts, pts[1:]):
            d = math.sqrt((a_[1] - b_[1]) ** 2 + (a_[2] - b_[2]) ** 2 + (a_[3] - b_[3]) ** 2)
            if d > 4.6:
                runs.append(cur)
                cur = []
            cur.append(b_)
        runs.append(cur)
        lo, hi = pts[0][0], pts[-1][0]
        for run in runs:
            if len(run) < 2:
                continue
            P = [run[0]] + run + [run[-1]]
            for i in range(1, len(P) - 2):
                p0, p1, p2, p3 = P[i - 1], P[i], P[i + 1], P[i + 2]
                for k in range(5):
                    t = k / 5.0
                    t2, t3 = t * t, t * t * t
                    xyz = [0.5 * ((2 * p1[j]) + (-p0[j] + p2[j]) * t + (2 * p0[j] - 5 * p1[j] + 4 * p2[j] - p3[j]) * t2 + (-p0[j] + 3 * p1[j] - 3 * p2[j] + p3[j]) * t3) for j in (1, 2, 3)]
                    kind = sec.get((ch, p1[0]), "C")
                    rad = {"H": 1.75, "E": 1.45, "C": 0.85}[kind]
                    if len(names) == 1:
                        hue = 0.72 * (1.0 - (p1[0] - lo) / max(1, hi - lo))
                        c = colorsys.hsv_to_rgb(hue, 0.78 if kind != "C" else 0.45, 0.95)
                    else:
                        base = CHAINS[ci % len(CHAINS)]
                        f = 1.0 if kind != "C" else 0.72
                        c = (base[0] * f, base[1] * f, base[2] * f)
                    out_pos += [xyz[0], xyz[1], xyz[2], rad]
                    out_col += list(c)
    return out_pos, out_col


def molecule_data(pdb_path, style="atoms"):
    """Atom spheres + colours of a PDB file. Returns (positions [x y z r ...], colours [r g b ...], number of spheres, bounding radius, number of real atoms)."""
    text = Path(pdb_path).read_text()
    atoms = []
    for line in text.splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        if line[17:20].strip() == "HOH" or line[16] not in " A":
            continue
        el_ = (line[76:78].strip() or line[12:14].strip()).upper()
        try:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            continue
        atoms.append((x, y, z, el_, line[21], int(line[22:26]), line.startswith("HETATM")))
    n = len(atoms)
    cx = sum(a[0] for a in atoms) / n
    cy = sum(a[1] for a in atoms) / n
    cz = sum(a[2] for a in atoms) / n
    chains = sorted({a[4] for a in atoms})
    ranges = {c: (min(a[5] for a in atoms if a[4] == c), max(a[5] for a in atoms if a[4] == c)) for c in chains}
    pos, col = [], []
    rmax = 0.0
    for (x, y, z, e, ch, resno, het) in atoms:
        x, y, z = x - cx, y - cy, z - cz
        r = RAD.get(e, 1.7)
        rmax = max(rmax, math.sqrt(x * x + y * y + z * z) + r)
        pos += [x, y, z, r]
        if e == "C":
            if len(chains) == 1:
                lo, hi = ranges[ch]
                c = colorsys.hsv_to_rgb(0.72 * (1.0 - (resno - lo) / max(1, hi - lo)), 0.55, 0.95)
            else:
                c = CHAINS[chains.index(ch) % len(CHAINS)]
                if het:
                    c = (0.55, 0.55, 0.58)
        else:
            c = CPK.get(e, (0.7, 0.7, 0.7))
        col += list(c)
    atom_count = n
    if style == "ribbon":
        pos, col = _ribbon_spheres(text, cx, cy, cz)
        n = len(pos) // 4
        rmax = max(math.sqrt(pos[4 * i] ** 2 + pos[4 * i + 1] ** 2 + pos[4 * i + 2] ** 2) + pos[4 * i + 3] for i in range(n))
    return pos, col, n, rmax, atom_count


def molecule_camera(rmax, az, el, count, width, height):
    dist = (rmax * 1.12) / 0.3226
    cam = (dist * math.cos(el) * math.cos(az), dist * math.cos(el) * math.sin(az), dist * math.sin(el))
    fw = _norm((-cam[0], -cam[1], -cam[2]))
    rt = _norm(_cross(fw, (0, 0, 1)))
    cu = _cross(rt, fw)
    return camera_block(cam, fw, rt, cu, count, 0.0, width, height)


def write_molecule(folder, pdb_path, az_deg, el_deg, width=1280, height=720, style="atoms"):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    pos, col, n, rmax, atoms = molecule_data(pdb_path, style)
    (folder / "atoms.bin").write_bytes(struct.pack("<%df" % len(pos), *pos))
    (folder / "colors.bin").write_bytes(struct.pack("<%df" % len(col), *col))
    (folder / "params.bin").write_bytes(molecule_camera(rmax, math.radians(az_deg), math.radians(el_deg), n, width, height))
    return {"spheres": n, "atoms": atoms, "rmax": rmax}


def write_molecule_camera(folder, rmax, spheres, az_deg, el_deg, width, height):
    """Only the camera changes (used for orbit films): the atom files stay."""
    (Path(folder) / "params.bin").write_bytes(molecule_camera(rmax, math.radians(az_deg), math.radians(el_deg), spheres, width, height))


# ---------------------------------------------------------------------------------------------- DNA double helix
TWIST, RISE, RADIUS, OFFSET = 0.5986479, 3.38, 10.0, 2.199115


def write_helix(folder, fasta_path, start, layers, az_deg, el_deg):
    """Idealised B-DNA (10.5 base pairs per turn, rise 3.38 angstrom, backbone radius 10 A) built from the real sequence."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    fa = Path(fasta_path).read_text().splitlines()
    seq = "".join(l.strip() for l in fa[1:]).upper()[start:start + layers]
    codes = []
    for ch in seq:
        c = {"A": 0, "C": 1, "G": 2, "T": 3, "U": 3}.get(ch)
        if c is not None:
            codes.append(c)
    nl = min(len(codes), 256)
    ax, ay, az_, cx, cy, bs = [[0.0] * 256 for _ in range(5)] + [[0] * 256]
    for i in range(nl):
        t = i * TWIST
        u = t + OFFSET
        ax[i] = int(RADIUS * math.cos(t) * 100.0) / 100.0
        ay[i] = int(RADIUS * math.sin(t) * 100.0) / 100.0
        az_[i] = int(i * RISE * 100.0) / 100.0
        cx[i] = int(RADIUS * math.cos(u) * 100.0) / 100.0
        cy[i] = int(RADIUS * math.sin(u) * 100.0) / 100.0
        bs[i] = codes[i]
    for name, arr in (("ax", ax), ("ay", ay), ("az", az_), ("cx", cx), ("cy", cy)):
        (folder / (name + ".bin")).write_bytes(struct.pack("<256f", *arr))
    (folder / "bs.bin").write_bytes(struct.pack("<256I", *bs))
    az, el = math.radians(az_deg), math.radians(el_deg)
    dist = 12.0 + nl * 3.38 * 1.5
    zc = (nl / 2.0) * 3.38
    tgt = (0.0, 0.0, zc)
    cam = (dist * math.cos(el) * math.cos(az), dist * math.cos(el) * math.sin(az), zc + dist * math.sin(el))
    fw = _norm((tgt[0] - cam[0], tgt[1] - cam[1], tgt[2] - cam[2]))
    rt = _norm(_cross(fw, (0, 0, 1)))
    cu = _cross(rt, fw)
    (folder / "params.bin").write_bytes(camera_block(cam, fw, rt, cu, nl, RISE))
    return {"base_pairs": nl, "sequence": seq[:nl]}


# ---------------------------------------------------------------------------------------------- schematic bacterium
def write_cell(folder, az_deg, el_deg):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    rng = random.Random(7)
    caps = []                                     # x1 y1 z1 x2 y2 z2 radius (micrometres)
    for k in range(3):                            # three helical flagella from the left pole
        phase = k * 2.0 * math.pi / 3.0
        prev = None
        for i in range(0, 70):
            s = i * 0.045
            amp = 0.05 + 0.04 * min(1.0, s / 0.5)
            pnt = (-0.52 - s, amp * math.cos(2 * math.pi * s / 0.9 + phase) - 0.04, amp * math.sin(2 * math.pi * s / 0.9 + phase) + 0.02 * k)
            if prev is not None and pnt[1] < 0.20:
                caps.append(prev + pnt + (0.011,))
            prev = pnt
    for k in range(16):                           # pili (fimbriae) on the kept half
        a = rng.uniform(math.pi, 2 * math.pi)
        x = rng.uniform(-0.45, 0.45)
        y, z = 0.5 * math.sin(a), 0.5 * math.cos(a)
        L = rng.uniform(0.22, 0.38)
        dx, dy, dz = rng.uniform(-0.25, 0.25), math.sin(a), math.cos(a)
        n = math.sqrt(dx * dx + dy * dy + dz * dz)
        caps.append((x, y, z, x + L * dx / n, y + L * dy / n, z + L * dz / n, 0.0065))
    flat = [v for c in caps for v in c]
    (folder / "cell_caps.bin").write_bytes(struct.pack("<%df" % len(flat), *flat))
    az, el = math.radians(az_deg), math.radians(el_deg)
    tgt, dist = (-0.45, 0.0, 0.0), 3.4
    cam = (tgt[0] + dist * math.cos(el) * math.cos(az), tgt[1] + dist * math.sin(el), tgt[2] + dist * math.cos(el) * math.sin(az))
    fw = _norm((tgt[0] - cam[0], tgt[1] - cam[1], tgt[2] - cam[2]))
    rt = _norm(_cross(fw, (0, 1, 0)))
    cu = _cross(rt, fw)
    (folder / "params.bin").write_bytes(camera_block(cam, fw, rt, cu, len(caps)))
    return {"capsules": len(caps)}


# ---------------------------------------------------------------------------------------------- circular genome map
MIN_AA = 150
STARTS = {"ATG", "GTG", "TTG"}
STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGT", "TGCA")
NB = 3600


def _scan(s):
    out = []
    for frame in range(3):
        start = None
        for i in range(frame, len(s) - 2, 3):
            c = s[i:i + 3]
            if c in STOPS:
                if start is not None and (i + 3 - start) // 3 - 1 >= MIN_AA:
                    out.append((start, i + 3))
                start = None
            elif start is None and c in STARTS:
                start = i
    return out


def write_genome(folder, fasta_path):
    """Open reading frames (six frames, >= 150 aa), GC content and smoothed GC skew in 3600 bins around the circle."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    lines = Path(fasta_path).read_text().splitlines()
    name = lines[0][1:].strip()
    seq = "".join(l.strip() for l in lines[1:]).upper()
    L = len(seq)
    fwd = [(a, b, 1) for a, b in _scan(seq)]
    rc = seq.translate(COMP)[::-1]
    rev = [(L - b, L - a, -1) for a, b in _scan(rc)]
    orfs = fwd + rev
    binw = L / NB
    cov = [[0.0] * NB, [0.0] * NB]
    for a, b, st in orfs:
        k = 0 if st == 1 else 1
        for i in range(int(a // binw), min(int((b - 1) // binw), NB - 1) + 1):
            lo, hi = max(a, i * binw), min(b, (i + 1) * binw)
            if hi > lo:
                cov[k][i] += (hi - lo) / binw
    gc, skew, tick = [], [], []
    for i in range(NB):
        w = seq[int(i * binw):int((i + 1) * binw)]
        g, c = w.count("G"), w.count("C")
        gc.append((g + c) / max(1, len(w)))
        skew.append((g - c) / max(1, g + c))
        mb = int(i * binw) // 500000 != int((i + 1) * binw) // 500000
        mn = int(i * binw) // 100000 != int((i + 1) * binw) // 100000
        tick.append(1.0 if mb else (0.5 if mn else 0.0))
    sm = []
    for i in range(NB):
        vals = [skew[j % NB] for j in range(i - 20, i + 21)]
        sm.append(sum(vals) / len(vals))
    data = []
    for i in range(NB):
        data += [min(1.0, cov[0][i]), min(1.0, cov[1][i]), gc[i], sm[i], tick[i]]
    (folder / "bins.bin").write_bytes(struct.pack("<%df" % len(data), *data))
    return {"name": name, "length": L, "orfs": len(orfs), "forward": len(fwd), "reverse": len(rev), "mean_gc": sum(gc) / NB, "min_aa": MIN_AA}

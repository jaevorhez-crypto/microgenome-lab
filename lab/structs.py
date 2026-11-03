"""Structure data for the interactive viewer and the medicine search. Reads real PDB files (public Protein Data Bank) and returns spheres: x y z radius r g b.
styles: atoms (space-filling), ribbon (Catmull-Rom backbone tube, helices and sheets thicker), pocket (ribbon of the protein + atoms of a bound ligand + the residues within a distance of it)."""
import colorsys
import math
from pathlib import Path

from gpu import APP
STRUCT = APP / "data" / "structures"
RAD = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "FE": 1.50, "H": 1.20, "P": 1.80, "ZN": 1.40, "MG": 1.45, "CA": 1.80, "NA": 1.70, "F": 1.47, "CL": 1.75, "BR": 1.85}
CPK = {"N": (0.25, 0.35, 0.95), "O": (0.95, 0.20, 0.15), "S": (0.96, 0.85, 0.25), "FE": (0.92, 0.45, 0.12), "P": (0.95, 0.60, 0.15), "H": (0.9, 0.9, 0.9), "F": (0.5, 0.9, 0.5), "CL": (0.2, 0.9, 0.3), "BR": (0.7, 0.2, 0.1)}
CHAINS = [(0.45, 0.75, 0.95), (0.96, 0.66, 0.36), (0.62, 0.86, 0.55), (0.86, 0.56, 0.86), (0.95, 0.85, 0.45), (0.55, 0.80, 0.80)]
SOLVENT = {"HOH", "WAT", "DOD"}


def parse(pid):
    atoms = []
    for line in (STRUCT / (pid + ".pdb")).read_text().splitlines():
        rec = line[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue
        res = line[17:20].strip()
        if res in SOLVENT or line[16] not in " A":
            continue
        try:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            resno = int(line[22:26])
        except ValueError:
            continue
        el = (line[76:78].strip() or line[12:14].strip()).upper()
        atoms.append({"x": x, "y": y, "z": z, "el": el, "name": line[12:16].strip(), "res": res, "chain": line[21], "resno": resno, "het": rec == "HETATM"})
    return atoms


def secondary(pid):
    sec = {}
    for ln in (STRUCT / (pid + ".pdb")).read_text().splitlines():
        try:
            if ln.startswith("HELIX "):
                ch, a, b = ln[19], int(ln[21:25]), int(ln[33:37])
                for r in range(a, b + 1):
                    sec[(ch, r)] = "H"
            elif ln.startswith("SHEET "):
                ch, a, b = ln[21], int(ln[22:26]), int(ln[33:37])
                for r in range(a, b + 1):
                    sec.setdefault((ch, r), "E")
        except ValueError:
            pass
    return sec


def ligands(pid):
    """bound small molecules (HETATM groups that are not water): {(res, chain, resno): atom count}"""
    out = {}
    for a in parse(pid):
        if a["het"]:
            out[(a["res"], a["chain"], a["resno"])] = out.get((a["res"], a["chain"], a["resno"]), 0) + 1
    return out


def center(atoms):
    n = len(atoms)
    return (sum(a["x"] for a in atoms) / n, sum(a["y"] for a in atoms) / n, sum(a["z"] for a in atoms) / n)


def ribbon(pid, atoms, c, faded=False):
    sec = secondary(pid)
    chs = {}
    for a in atoms:
        if a["name"] == "CA" and not a["het"]:
            chs.setdefault(a["chain"], []).append((a["resno"], a["x"] - c[0], a["y"] - c[1], a["z"] - c[2]))
    out = []
    names = sorted(chs)
    for ci, ch in enumerate(names):
        pts = chs[ch]
        if len(pts) < 2:
            continue
        runs, cur = [], [pts[0]]
        for p_, q_ in zip(pts, pts[1:]):
            if math.dist(p_[1:], q_[1:]) > 4.6:
                runs.append(cur)
                cur = []
            cur.append(q_)
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
                        col = colorsys.hsv_to_rgb(hue, 0.78 if kind != "C" else 0.45, 0.95)
                    else:
                        base = CHAINS[ci % len(CHAINS)]
                        f = 1.0 if kind != "C" else 0.72
                        col = (base[0] * f, base[1] * f, base[2] * f)
                    if faded:
                        col = tuple(0.35 + 0.35 * v for v in col)
                    out += [xyz[0], xyz[1], xyz[2], rad, col[0], col[1], col[2]]
    return out


def atom_spheres(atoms, c, mode="chain", chains=None, lig_color=None):
    out = []
    chains = chains or sorted({a["chain"] for a in atoms})
    rng = {ch: (min(a["resno"] for a in atoms if a["chain"] == ch), max(a["resno"] for a in atoms if a["chain"] == ch)) for ch in chains}
    for a in atoms:
        e = a["el"]
        if e == "C":
            if lig_color is not None and a["het"]:
                col = lig_color
            elif len(chains) == 1:
                lo, hi = rng[a["chain"]]
                col = colorsys.hsv_to_rgb(0.72 * (1.0 - (a["resno"] - lo) / max(1, hi - lo)), 0.55, 0.95)
            else:
                col = (0.55, 0.55, 0.58) if a["het"] else CHAINS[chains.index(a["chain"]) % len(CHAINS)]
        else:
            col = CPK.get(e, (0.7, 0.7, 0.7))
        out += [a["x"] - c[0], a["y"] - c[1], a["z"] - c[2], RAD.get(e, 1.7), col[0], col[1], col[2]]
    return out


def structure(pid, style="ribbon"):
    atoms = parse(pid)
    c = center(atoms)
    if style == "atoms":
        data = atom_spheres(atoms, c)
    else:
        data = ribbon(pid, atoms, c)
    n = len(data) // 7
    rmax = max(math.sqrt(data[7 * i] ** 2 + data[7 * i + 1] ** 2 + data[7 * i + 2] ** 2) + data[7 * i + 3] for i in range(n))
    return {"n": n, "rmax": rmax, "data": data, "title": ""}


def split_ligand(atoms, key):
    """(ligand atoms, protein atoms). key = (residue name, chain, number); residue name '*' means the whole chain is the ligand (a peptide inhibitor)."""
    if key[0] == "*":
        return [a for a in atoms if a["chain"] == key[1]], [a for a in atoms if a["chain"] != key[1] and not a["het"]]
    return [a for a in atoms if a["het"] and (a["res"], a["chain"], a["resno"]) == key], [a for a in atoms if not a["het"]]


def contacts(atoms, lig_key, cutoff=4.0):
    """atom pairs between a bound ligand and the protein within `cutoff` angstrom; N/O pairs under 3.5 A are listed as possible hydrogen bonds"""
    lig, prot = split_ligand(atoms, lig_key)
    rows = []
    for l in lig:
        for p in prot:
            if abs(l["x"] - p["x"]) > cutoff or abs(l["y"] - p["y"]) > cutoff or abs(l["z"] - p["z"]) > cutoff:
                continue
            d = math.dist((l["x"], l["y"], l["z"]), (p["x"], p["y"], p["z"]))
            if d <= cutoff:
                hb = (l["el"] in ("N", "O") and p["el"] in ("N", "O") and d <= 3.5)
                rows.append({"lig_atom": l["name"], "lig_res": "%s%d" % (l["res"], l["resno"]), "res": "%s%d" % (p["res"], p["resno"]), "chain": p["chain"], "prot_atom": p["name"], "dist": round(d, 2), "hbond": hb})
    rows.sort(key=lambda r: r["dist"])
    return rows


def pocket(pid, lig_key, radius=5.0):
    """ribbon of the nearby protein chains (faded) + atoms of the ligand (bright yellow carbons) + all atoms of residues within `radius` of it"""
    atoms = parse(pid)
    lig, prot = split_ligand(atoms, lig_key)
    c = center(lig)
    near = set()
    near_chains = set()
    lpts = [(l["x"], l["y"], l["z"]) for l in lig]
    for p in prot:
        dmin = min(math.dist(q, (p["x"], p["y"], p["z"])) for q in lpts)
        if dmin <= radius:
            near.add((p["chain"], p["resno"]))
        if dmin <= 14.0:
            near_chains.add(p["chain"])
    res_atoms = [a for a in prot if (a["chain"], a["resno"]) in near]
    rib_atoms = [a for a in atoms if a["chain"] in near_chains and not a["het"] and a["chain"] != (lig_key[1] if lig_key[0] == "*" else "")]
    data = ribbon(pid, rib_atoms, c, faded=True)
    data += atom_spheres(res_atoms, c, chains=sorted({a["chain"] for a in res_atoms}))
    data += atom_spheres(lig, c, lig_color=(1.0, 0.85, 0.2))
    if lig_key[0] == "*":                                        # a peptide inhibitor: its backbone and side chains are all 'het' for colouring
        pass
    n = len(data) // 7
    rmax = max(math.sqrt(data[7 * i] ** 2 + data[7 * i + 1] ** 2 + data[7 * i + 2] ** 2) + data[7 * i + 3] for i in range(n))
    return {"n": n, "rmax": rmax, "data": data, "pocket_residues": len(near), "ligand_atoms": len(lig), "residue_list": sorted("%s%d" % (next(a["res"] for a in res_atoms if (a["chain"], a["resno"]) == k), k[1]) for k in near)}

"""MicroGenome Lab: one window for the whole toolkit. A small local web server (127.0.0.1 only) and one page.
The pictures, the sequence search and the aligner run as small GPU/CPU programs from the bin folder (see gpu.py).
Safety: listens on localhost only, a fixed whitelist of actions, inputs are cleaned (sequences to ACGT, PDB ids from a list, numbers clamped), subprocesses use argument lists
(no shell strings from the browser), GPU jobs run one at a time. usage: MicroGenomeLab [port]"""
import hashlib
import json
import os
import re
import subprocess
import sys
import threading
import time
import traceback
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, unquote, urlparse

import numpy as np

import gpu
import scenes

APP = gpu.APP
HERE = Path(__file__).resolve().parent
DB = APP / "data"
OUT = APP / "out"
STATIC = (Path(getattr(sys, "_MEIPASS", HERE)) / "static") if getattr(sys, "frozen", False) else HERE / "static"
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(HERE))
import sketchlib

PORT = next((int(a) for a in sys.argv[1:] if a.isdigit()), 8770)
GPU_LOCK = threading.Lock()
NAMES = json.load((DB / "index" / "names.json").open())
LINEAGE = json.load((DB / "index" / "lineage.json").open())
META = json.load((DB / "index" / "meta.json").open())
GENOME_META = json.loads((OUT / "genome_meta.json").read_text()) if (OUT / "genome_meta.json").exists() else {}
PDBS = {"1LYZ": "Lysozyme (hen egg white): an enzyme that cuts bacterial cell walls, 129 residues",
        "4HHB": "Hemoglobin (human): oxygen transport, four chains with heme groups",
        "1EMA": "Green fluorescent protein (jellyfish): a barrel that glows, 238 residues",
        "6VXX": "SARS-CoV-2 spike glycoprotein (closed state): 3 chains with sugars"}
SEARCH_DIR = DB / "index"
NQMAX = 4096


def clean_seq(text, maxlen=200000):
    lines = [l for l in (text or "").splitlines() if not l.startswith(">")]
    s = re.sub(r"[^A-Za-z]", "", "".join(lines)).upper().replace("U", "T")
    return re.sub(r"[^ACGTN]", "N", s)[:maxlen]


def clampf(v, lo, hi, default):
    try:
        return max(lo, min(hi, float(v)))
    except (TypeError, ValueError):
        return default


def clampi(v, lo, hi, default):
    return int(clampf(v, lo, hi, default))


def cache_key(*parts):
    return hashlib.sha1(json.dumps(parts, sort_keys=True).encode()).hexdigest()[:16]


def film_to_video(frames, out_base, fps=20):
    """Frames (list of RGB arrays) -> mp4 when ffmpeg is installed, otherwise an animated GIF made with Pillow."""
    import shutil
    from PIL import Image
    ff = shutil.which("ffmpeg")
    if ff:
        tmp = OUT / ("frames_" + out_base)
        tmp.mkdir(exist_ok=True)
        for i, a in enumerate(frames):
            Image.fromarray(a, "RGB").save(tmp / ("f%03d.png" % i))
        subprocess.run([ff, "-y", "-loglevel", "error", "-framerate", str(fps), "-i", str(tmp / "f%03d.png"), "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "19", str(OUT / (out_base + ".mp4"))],
                       check=True, timeout=300, creationflags=gpu.NO_WINDOW)
        for f in tmp.glob("*.png"):
            f.unlink()
        tmp.rmdir()
        return out_base + ".mp4"
    ims = [Image.fromarray(a, "RGB").convert("P", palette=Image.ADAPTIVE, colors=128) for a in frames]
    ims[0].save(OUT / (out_base + ".gif"), save_all=True, append_images=ims[1:], duration=int(1000 / fps), loop=0)
    return out_base + ".gif"


def do_identify(body):
    seq = clean_seq(body.get("sequence", ""))
    if len(seq) < 60:
        return {"error": "Give at least 60 bases (a FASTA header line is fine)."}
    sk = sketchlib.sketch(seq.encode())
    if len(sk) < 3:
        return {"error": "That sequence is too short or too ambiguous to make a sketch (it needs at least a few hundred clean bases)."}
    q = np.full(NQMAX, 0xFFFFFFFF, dtype="<u4")
    n = min(len(sk), NQMAX)
    q[:n] = sk[:n]
    t0 = time.time()
    with GPU_LOCK:
        q.tofile(SEARCH_DIR / "query.bin")
        gpu.run("search", SEARCH_DIR, timeout=180, expect="hits.bin")
        hits = np.fromfile(SEARCH_DIR / "hits.bin", dtype="<u4", count=len(NAMES))
    top = np.argsort(-hits.astype(np.int64), kind="stable")[:10]
    res = []
    for i in top:
        h = int(hits[i])
        if h == 0:
            continue
        nm = NAMES[int(i)]
        lin = LINEAGE.get(nm["id"], {})
        res.append({"id": nm["id"], "title": nm["title"], "kind": nm["kind"], "length": nm["len"], "hits": h, "containment": round(h / n, 3),
                    "lineage": lin.get("lineage", []), "taxon": lin.get("name", "")})
    best = res[0]["containment"] if res else 0.0
    verdict = ("No confident match: the best containment is %.3f; an unrelated sequence scores close to 0." % best) if best < 0.10 else ("Best match: %s (containment %.2f)" % (res[0]["title"], best))
    return {"results": res, "sketch_size": int(n), "bases": len(seq), "seconds": round(time.time() - t0, 2), "verdict": verdict,
            "note": "Scaled-MinHash containment against %d RefSeq sequences (%d viral genomes + %d 16S rRNA). Close relatives and strains of one species can swap places; this is a screen, not a diagnosis." % (META["sequences"], META["viruses"], META["sixteen_s"])}


def do_seqtools(body):
    seq = clean_seq(body.get("sequence", ""), 500000)
    if not seq:
        return {"error": "No sequence."}
    L = len(seq)
    cnt = {c: seq.count(c) for c in "ACGTN"}
    gc = (cnt["G"] + cnt["C"]) / max(1, cnt["A"] + cnt["C"] + cnt["G"] + cnt["T"])
    comp = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    table = None
    codon = {}
    bases = "TCAG"
    aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    k = 0
    for a in bases:
        for b in bases:
            for c in bases:
                codon[a + b + c] = aas[k]
                k += 1
    orfs = []
    for strand, s in (("+", seq), ("-", comp)):
        for frame in range(3):
            start = None
            for i in range(frame, len(s) - 2, 3):
                c = s[i:i + 3]
                if c in ("TAA", "TAG", "TGA"):
                    if start is not None and (i + 3 - start) // 3 - 1 >= 50:
                        prot = "".join(codon.get(s[j:j + 3], "X") for j in range(start, i, 3))
                        a, b = (start + 1, i + 3) if strand == "+" else (L - (i + 3) + 1, L - start)
                        orfs.append({"strand": strand, "start": a, "end": b, "aa": (i - start) // 3, "protein": prot})
                    start = None
                elif start is None and c in ("ATG", "GTG", "TTG"):
                    start = i
    orfs.sort(key=lambda o: -o["aa"])
    return {"length": L, "counts": cnt, "gc_percent": round(gc * 100, 2), "reverse_complement": comp[:600], "orfs": orfs[:8], "orf_total": len(orfs),
            "note": "ORFs are candidates (start codon ATG/GTG/TTG to stop, at least 50 amino acids, six frames), not curated gene annotation."}


def do_align(body):
    ref, qry = clean_seq(body.get("ref", ""), 3500), clean_seq(body.get("qry", ""), 3500)
    if len(ref) < 5 or len(qry) < 5:
        return {"error": "Give two sequences of at least 5 bases (each up to 3,500)."}
    t0 = time.time()
    p = subprocess.run([str(gpu.program("align"))], input=(">ref\n%s\n>qry\n%s\n" % (ref, qry)).encode(), capture_output=True, timeout=170, creationflags=gpu.NO_WINDOW)
    out = p.stdout.decode("utf-8", "replace")
    if not out.startswith("SCORE"):
        return {"error": "The aligner did not run: " + out[-300:]}
    return {"output": out, "seconds": round(time.time() - t0, 2), "note": "Global alignment with affine gaps (match +2, mismatch -3, gap 5 then 2 per extra base). The score is the true optimum."}


def result(kind, png=None, mp4=None, note="", seconds=0.0):
    d = {"kind": kind, "note": note, "seconds": round(seconds, 2)}
    if png:
        d["image"] = "/out/" + png
    if mp4:
        d["video"] = "/out/" + mp4
    return d


def render_job(program, folder, width, height, png):
    gpu.run(program, folder, timeout=300, expect="frame.raw")
    gpu.frame_to_png(folder / "frame.raw", width, height, OUT / png)


def do_render(body):
    what = body.get("what", "")
    t0 = time.time()
    if what == "helix":
        start = clampi(body.get("start"), 0, 5300, 0)
        layers = clampi(body.get("layers"), 8, 36, 16)
        az = clampf(body.get("az"), 0, 360, 50)
        el = clampf(body.get("el"), -30, 60, 10)
        custom = clean_seq(body.get("sequence", ""), 400)
        key = cache_key(what, start, layers, az, el, custom)
        png = key + ".png"
        if not (OUT / png).exists():
            fasta = DB / "genomes" / "phix174.txt"
            if len(custom) >= 8:
                fasta = OUT / (key + ".fa")
                fasta.write_text(">your sequence\n" + custom + "\n")
                start = 0
            job = OUT / ("job_" + key)
            scenes.write_helix(job, fasta, start, layers, az, el)
            render_job("helix", job, 1280, 720, png)
        return result("helix", png, note="Idealised B-DNA (3.38 A rise, 34.3 degrees twist) built from the sequence. %s" % ("Your sequence." if len(custom) >= 8 else "Real phiX174 sequence."), seconds=time.time() - t0)
    if what in ("protein", "protein_orbit"):
        pid = body.get("pdb", "4HHB")
        if pid not in PDBS:
            return {"error": "Unknown structure."}
        style = "ribbon" if body.get("style") == "ribbon" else "atoms"
        az = clampf(body.get("az"), 0, 360, 40)
        el = clampf(body.get("el"), -40, 70, 20)
        pdb = DB / "structures" / (pid + ".pdb")
        if what == "protein":
            key = cache_key(what, pid, style, az, el)
            png = key + ".png"
            if not (OUT / png).exists():
                job = OUT / ("job_" + key)
                scenes.write_molecule(job, pdb, az, el, 1280, 720, style)
                render_job("mol", job, 1280, 720, png)
            return result("protein", png, note="%s. Real coordinates from the public PDB file %s." % (PDBS[pid], pid), seconds=time.time() - t0)
        key = cache_key(what, pid, style, el)
        done = [f for f in (key + ".mp4", key + ".gif") if (OUT / f).exists()]
        if not done:
            job = OUT / ("job_" + key)
            info = scenes.write_molecule(job, pdb, 0, el, 640, 360, style)
            frames = []
            for i in range(60):
                scenes.write_molecule_camera(job, info["rmax"], info["spheres"], i * 6.0, el, 640, 360)
                gpu.run("mol", job, timeout=120, expect="frame.raw")
                frames.append(gpu.frame_to_png(job / "frame.raw", 640, 360, job / "f.png"))
            done = [film_to_video(frames, key)]
        return result("protein_orbit", mp4=done[0], note="%s. 60 frames, one full turn." % PDBS[pid], seconds=time.time() - t0)
    if what == "cell":
        az = clampf(body.get("az"), 0, 360, 60)
        el = clampf(body.get("el"), 5, 70, 28)
        key = cache_key(what, az, el)
        png = key + ".png"
        if not (OUT / png).exists():
            job = OUT / ("job_" + key)
            scenes.write_cell(job, az, el)
            render_job("cell", job, 1280, 720, png)
        return result("cell", png, note="SCHEMATIC rod-shaped bacterium (about 2 x 0.5 micrometre, E. coli-like): the organisation is illustrative, the envelope layers are drawn about 2x thicker than life. Not a measured structure.", seconds=time.time() - t0)
    if what == "genome":
        png = "genome_ecoli.png"
        global GENOME_META
        if not (OUT / png).exists() or not GENOME_META:
            job = OUT / "job_genome"
            GENOME_META = scenes.write_genome(job, DB / "structures" / "NC_000913.3.fasta")
            (OUT / "genome_meta.json").write_text(json.dumps(GENOME_META))
            render_job("genome", job, 1200, 1200, png)
        gm = GENOME_META
        return result("genome", png, note="E. coli K-12 MG1655 (NC_000913.3), %s bases: %s candidate ORFs of at least %s amino acids (E. coli has about 4,300 annotated genes, so some are false). Rings from outside: ruler (ticks every 100 kb, long ticks every 500 kb), forward ORFs, reverse ORFs, GC content, GC skew." % (
            format(gm.get("length", 0), ","), format(gm.get("orfs", 0), ","), gm.get("min_aa", "")), seconds=time.time() - t0)
    return {"error": "Unknown action."}


def do_tree(body):
    import evo
    term = str(body.get("term", ""))[:60].strip()
    kind = body.get("kind", "")
    kind = kind if kind in ("virus", "16S") else ""
    limit = clampi(body.get("limit"), 4, 40, 18)
    extra = None
    seq = clean_seq(body.get("sequence", ""))
    if len(seq) >= 200:
        sk = sketchlib.sketch(seq.encode())
        if len(sk) >= 3:
            extra = ("YOUR SEQUENCE", sk)
    idxs = evo.members(term, kind, limit) if term else []
    if len(idxs) + (1 if extra else 0) < 3:
        return {"error": "Fewer than 3 matching sequences. Try a broader name such as Betacoronavirus, Herpesviridae, Escherichia, Enterobacteriaceae or Bacillus."}
    tree, D, names = evo.tree_for(idxs, extra)
    far = sum(1 for i in range(len(D)) for j in range(i + 1, len(D)) if D[i][j] >= 0.999)
    note = "Neighbour-joining tree from Mash-style distances between scaled sketches (k=21). Branch length is estimated substitutions per site. Reliable for relatives within a family or genus; pairs with no shared sketch saturate at 1.0 (%d such pairs here)." % far
    if kind == "16S" or (idxs and evo._names[idxs[0]]["kind"] == "16S"):
        note += " 16S sketches hold only about 30 hashes, so this tree is coarse."
    return {"tree": tree, "leaves": len(names), "note": note}


def do_assemble(body):
    import asm
    t0 = time.time()
    mode = body.get("mode", "simulate")
    reference = None
    if mode == "simulate":
        ref = body.get("reference", "phix174")
        files = {"phix174": "phix174.txt", "lambda": "lambda_phage.txt", "sars": "sars_cov_2.txt", "hbv": "hbv_genotype_d.txt"}
        if ref not in files:
            return {"error": "Unknown reference."}
        data = DB / "genomes" / files[ref]
        reference = "".join(l.strip() for l in data.read_text().splitlines() if not l.startswith(">")).upper()
        cov = clampi(body.get("coverage"), 10, 100, 30)
        err = clampf(body.get("error"), 0.0, 3.0, 0.3) / 100.0
        reads = asm.simulate_reads(reference, 150, cov, err, seed=clampi(body.get("seed"), 1, 99999, 1))
        src = "%d simulated reads of 150 bases from the real %s genome (%d bases) at %dx coverage with %.1f%% sequencing errors" % (len(reads), ref, len(reference), cov, err * 100)
    else:
        reads = asm.parse_reads(str(body.get("reads", ""))[:6_000_000])
        if len(reads) < 20:
            return {"error": "Give at least 20 reads (FASTA or FASTQ text, up to about 6 MB)."}
        src = "%d reads you provided" % len(reads)
    contigs, info = asm.assemble(reads, k=clampi(body.get("k"), 21, 41, 31))
    st = asm.stats(contigs)
    out = {"source": src, "info": info, "stats": st, "contigs": [{"length": len(c), "start": c[:90], "gc": round(100 * (c.count("G") + c.count("C")) / len(c), 1)} for c in contigs[:12]],
           "seconds": round(time.time() - t0, 2)}
    if reference:
        out["check"] = asm.check_against(contigs, reference)
        out["verdict"] = "%d contig(s); %.2f%% of the real genome's k-mers are recovered; the longest contig %s the real genome (either strand)." % (st["contigs"], 100 * out["check"]["reference_kmers_covered"], "matches exactly" if out["check"]["longest_contigs_exact"] and out["check"]["longest_contigs_exact"][0][1] else "does NOT exactly match")
    else:
        out["verdict"] = "%d contig(s), total %d bases, longest %d, N50 %d." % (st["contigs"], st["total"], st["longest"], st["n50"])
    out["note"] = "De Bruijn graph of canonical k-mers (k=%d), k-mers seen fewer than %d time(s) dropped as likely errors, then unbranched paths. This build runs the assembler in Python." % (info["k"], info["min_coverage"])
    return out


def do_structure(query):
    import base64
    import struct as _st
    import structs
    pid = query.get("pdb", ["4HHB"])[0]
    style = query.get("style", ["ribbon"])[0]
    if pid not in PDBS and pid not in MEDS_PDB:
        return {"error": "Unknown structure."}
    if style == "pocket":
        key = query.get("lig", [""])[0].split("|")
        if len(key) != 3:
            return {"error": "Give lig=RES|CHAIN|NUMBER."}
        d = structs.pocket(pid, (key[0], key[1], int(key[2])))
    else:
        d = structs.structure(pid, "atoms" if style == "atoms" else "ribbon")
    raw = _st.pack("<%df" % len(d["data"]), *d["data"])
    d["data"] = base64.b64encode(raw).decode()
    return d


MEDS = {
    "6LU7": {"name": "SARS-CoV-2 main protease (Mpro, 3CL protease)", "disease": "COVID-19", "lig": ["*", "C", 0], "lig_name": "peptide inhibitor N3 (a research compound, bound covalently to the catalytic cysteine)",
             "what": "The virus makes its proteins as one long chain; this enzyme cuts the chain into working pieces. Block it and the virus cannot copy itself.",
             "medicines": [("nirmatrelvir", "the protease inhibitor in Paxlovid: a peptide-like molecule that binds this same pocket and blocks the catalytic cysteine")],
             "chembl": "CHEMBL3927"},
    "1HVR": {"name": "HIV-1 protease", "disease": "HIV/AIDS", "lig": ["XK2", "A", 263], "lig_name": "DMP 323, a research inhibitor (cyclic urea); not one of the marketed drugs",
             "what": "HIV builds its core from a chain that this dimeric enzyme cuts apart. Two aspartates, one from each half, do the cutting.",
             "medicines": [("ritonavir", "approved HIV protease inhibitor"), ("lopinavir", "approved HIV protease inhibitor (given with ritonavir)"), ("darunavir", "approved HIV protease inhibitor"), ("atazanavir", "approved HIV protease inhibitor")],
             "chembl": "CHEMBL243"},
    "2HU4": {"name": "Influenza N1 neuraminidase (H5N1 strain)", "disease": "influenza", "lig": ["G39", "A", 800], "lig_name": "oseltamivir (its active form), the medicine sold as Tamiflu",
             "what": "New flu particles are stuck to the infected cell by sugars; neuraminidase cuts them free so they can spread. Block it and the virus stays trapped.",
             "medicines": [("oseltamivir", "approved neuraminidase inhibitor (Tamiflu)"), ("zanamivir", "approved neuraminidase inhibitor (inhaled)"), ("peramivir", "approved neuraminidase inhibitor (injected)")],
             "chembl": ""},
    "1M17": {"name": "EGFR tyrosine kinase domain", "disease": "some lung and other cancers", "lig": ["AQ4", "A", 999], "lig_name": "erlotinib, an approved cancer medicine",
             "what": "EGFR is a switch on the cell surface that tells cells to grow. In some cancers it is stuck on. Erlotinib sits in the kinase's ATP pocket and keeps the switch off.",
             "medicines": [("erlotinib", "approved EGFR inhibitor (the one in this structure)"), ("gefitinib", "approved EGFR inhibitor"), ("afatinib", "approved EGFR inhibitor"), ("osimertinib", "approved EGFR inhibitor for resistant mutations")],
             "chembl": "CHEMBL203"},
}
MEDS_PDB = set(MEDS)


def do_med(query):
    import structs
    pid = query.get("pdb", [""])[0]
    m = MEDS.get(pid)
    if not m:
        return {"error": "Unknown target."}
    key = (m["lig"][0], m["lig"][1], int(m["lig"][2]))
    atoms = structs.parse(pid)
    rows = structs.contacts(atoms, key)
    pk = structs.pocket(pid, key)
    hb = [r for r in rows if r["hbond"]]
    return {"pdb": pid, "name": m["name"], "disease": m["disease"], "what": m["what"], "lig_name": m["lig_name"], "lig_key": "|".join(str(x) for x in m["lig"]),
            "medicines": [{"name": a, "note": b} for a, b in m["medicines"]], "pocket_residues": pk["residue_list"], "ligand_atoms": pk["ligand_atoms"], "contacts_total": len(rows),
            "hbonds": hb[:14], "closest": rows[:10], "chembl": m["chembl"],
            "note": "Every number here is measured from the real structure file (distances in angstrom; a nitrogen or oxygen pair within 3.5 A is listed as a possible hydrogen bond). This shows what is known. It does not predict new medicines, doses or safety: finding a drug needs laboratory and clinical testing."}


def do_potency(body):
    """live ChEMBL lookup of the most potent measured compounds for a target; honest about outages"""
    import urllib.request
    pid = body.get("pdb", "")
    m = MEDS.get(pid)
    if not m:
        return {"error": "Unknown target."}
    pack_file = DB / "medpack" / "medicine_pack.json"
    if pack_file.exists() and not body.get("refresh"):
        pack = json.load(pack_file.open())
        t = pack.get("targets", {}).get(pid)
        if t and t.get("potent"):
            src = pack["source"]
            return {"rows": t["potent"], "local": True, "note": "Measured IC50 values (lower is more potent) for ChEMBL target %s, from the lab's local snapshot taken %s (no internet needed). Source: %s, licence %s, attribution required. Different assays are not directly comparable; a few very low values may be reporting errors in the original papers." % (t["chembl_target"], src["fetched"], src["name"], src["license"])}
    tid = m["chembl"]
    if not tid:
        return {"error": "No ChEMBL target is linked for this structure yet."}
    url = "https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id=%s&standard_type=IC50&standard_units=nM&pchembl_value__gte=8&order_by=-pchembl_value&limit=12" % tid
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "MicroGenomeLab"}), timeout=20) as r:
            acts = json.load(r).get("activities", [])
    except Exception as e:
        return {"error": "ChEMBL (the public database of measured compound activity) is not answering right now (%s). Nothing is shown rather than guessing numbers; try again later." % type(e).__name__}
    rows = [{"molecule": a.get("molecule_chembl_id"), "name": a.get("molecule_pref_name") or "", "ic50_nM": a.get("standard_value"), "pchembl": a.get("pchembl_value"), "assay": a.get("assay_description", "")[:110]} for a in acts]
    return {"rows": rows, "note": "Measured IC50 values (lower is more potent) from ChEMBL target %s, as published in the literature. Different assays are not directly comparable." % tid}


def example(name):
    import random
    rng = random.Random(20261102)
    data = DB / "genomes"
    def fa(f):
        return "".join(l.strip() for l in (data / f).read_text().splitlines() if not l.startswith(">")).upper()
    def mut(sq, rate):
        return "".join(rng.choice([x for x in "ACGT" if x != c]) if (rng.random() < rate and c in "ACGT") else c for c in sq)
    nl = chr(10)
    if name == "ecoli16s":
        return {"sequence": ">E. coli 16S rRNA with 3% random substitutions" + nl + mut(fa("16s_escherichia_coli.txt"), 0.03)}
    if name == "sars":
        sq = fa("sars_cov_2.txt")
        return {"sequence": ">piece of SARS-CoV-2 spike region with 2% substitutions" + nl + mut(sq[22400:24400], 0.02)}
    if name == "hbv":
        return {"sequence": ">piece of a hepatitis B genome with 2% substitutions" + nl + mut(fa("hbv_genotype_d.txt")[500:2500], 0.02)}
    return {"sequence": ">random DNA (should match nothing)" + nl + "".join(rng.choice("ACGT") for _ in range(1500))}

ROUTES = {"/api/identify": do_identify, "/api/seqtools": do_seqtools, "/api/align": do_align, "/api/render": do_render, "/api/tree": do_tree, "/api/assemble": do_assemble, "/api/potency": do_potency}
MIME = {".html": "text/html; charset=utf-8", ".js": "text/javascript", ".css": "text/css", ".png": "image/png", ".mp4": "video/mp4", ".json": "application/json", ".bmp": "image/bmp", ".gif": "image/gif"}


class H(BaseHTTPRequestHandler):
    def log_message(self, *a):
        pass

    def _send(self, code, data, ctype="application/json"):
        if isinstance(data, (dict, list)):
            data = json.dumps(data).encode()
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        path = unquote(urlparse(self.path).path)
        if path == "/api/info":
            return self._send(200, {"index": META, "pdbs": PDBS, "genome": GENOME_META, "edition": gpu.edition(), "gpu": gpu.gpu_name})
        if path.endswith(".js") and "/" not in path[1:] and (STATIC / path[1:]).is_file():
            return self._send(200, (STATIC / path[1:]).read_bytes(), MIME[".js"])
        if path == "/api/meds":
            return self._send(200, [{"pdb": k, "name": v["name"], "disease": v["disease"]} for k, v in MEDS.items()])
        if path == "/api/med":
            return self._send(200, do_med(parse_qs(urlparse(self.path).query)))
        if path == "/api/structure":
            return self._send(200, do_structure(parse_qs(urlparse(self.path).query)))
        if path == "/api/example":
            name = parse_qs(urlparse(self.path).query).get("name", [""])[0]
            return self._send(200, example(name))
        if path in ("/", "/index.html"):
            return self._send(200, (STATIC / "index.html").read_bytes(), MIME[".html"])
        for base, root in (("/out/", OUT),):
            if path.startswith(base):
                f = (root / path[len(base):]).resolve()
                if root.resolve() in f.parents and f.is_file() and f.suffix in MIME:
                    return self._send(200, f.read_bytes(), MIME[f.suffix])
        self._send(404, {"error": "not found"})

    def do_POST(self):
        fn = ROUTES.get(urlparse(self.path).path)
        if not fn:
            return self._send(404, {"error": "not found"})
        try:
            n = int(self.headers.get("Content-Length", "0"))
            body = json.loads(self.rfile.read(min(n, 2_000_000)) or b"{}")
            if fn is do_render:
                with GPU_LOCK:
                    res = fn(body)
            else:
                res = fn(body)
            self._send(200, res)
        except subprocess.CalledProcessError as e:
            self._send(200, {"error": "A build step failed: " + ((e.stderr or b"")[-300:].decode("latin-1") if isinstance(e.stderr, bytes) else str(e.stderr)[-300:])})
        except Exception as e:
            traceback.print_exc()
            self._send(200, {"error": "%s: %s" % (type(e).__name__, e)})


if __name__ == "__main__":
    if not (DB / "index" / "hashes.bin").exists():
        print("The data pack is missing: unpack MicroGenomeLab-data.zip next to this program so that the folder 'data' exists.")
        sys.exit(1)
    srv = ThreadingHTTPServer(("127.0.0.1", PORT), H)
    url = "http://127.0.0.1:%d" % PORT
    print("MicroGenome Lab (%s edition) is running at %s  (this PC only). Close this window to stop it." % (gpu.edition(), url), flush=True)
    if "--no-browser" not in sys.argv:
        import webbrowser
        threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    srv.serve_forever()

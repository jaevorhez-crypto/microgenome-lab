"""GPU programs of MicroGenome Lab. The picture makers, the sketch search and the aligner are small stand-alone programs (bin folder) that read their input
from files in a working folder and write their result there. This module starts them and turns their raw pictures into PNG files."""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

APP = Path(sys.executable).parent if getattr(sys, "frozen", False) else Path(__file__).resolve().parent.parent
BIN = Path(os.environ.get("MGL_BIN") or (APP / "bin"))
NO_WINDOW = 0x08000000 if os.name == "nt" else 0
gpu_name = ""


def edition():
    f = BIN / "edition.txt"
    return f.read_text().strip() if f.exists() else "unknown"


def program(name):
    p = BIN / (name + (".exe" if os.name == "nt" else ""))
    if not p.exists():
        raise RuntimeError("The program '%s' is missing from %s." % (name, BIN))
    return p


def run(name, cwd, timeout=300, expect=None, stdin=None):
    """Run one program in `cwd`. `expect` is the output file that must exist afterwards."""
    global gpu_name
    cwd = Path(cwd)
    if expect and (cwd / expect).exists():
        (cwd / expect).unlink()
    p = subprocess.run([str(program(name))], cwd=str(cwd), input=stdin, capture_output=True, timeout=timeout, creationflags=NO_WINDOW)
    out = p.stdout.decode("utf-8", "replace") if p.stdout else ""
    for line in out.splitlines():
        if line.startswith("GPU: "):
            gpu_name = line[5:].strip()
    if expect and not (cwd / expect).exists():
        raise RuntimeError("The GPU program '%s' produced no result (exit code %s). Check that a supported graphics card and its current driver are installed. %s" % (name, p.returncode, (p.stderr or b"")[-200:].decode("utf-8", "replace")))
    return out


def frame_to_png(raw, width, height, png):
    from PIL import Image
    a = np.fromfile(raw, dtype="<u4", count=width * height)
    rgb = np.stack([(a >> 16) & 255, (a >> 8) & 255, a & 255], axis=-1).astype(np.uint8).reshape(height, width, 3)
    Image.fromarray(rgb, "RGB").save(png)
    return rgb

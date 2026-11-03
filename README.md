# MicroGenome Lab

**A free, offline, one-window workbench for looking at DNA, proteins, microbes and the medicines that act on them, accelerated by your own graphics card.**

MicroGenome Lab runs entirely on your own Windows PC. Nothing you paste in is uploaded anywhere, it needs no account and no internet connection once downloaded, and it is free.

It is an **analysis and learning tool**. It reads and draws real, public data. It does not design, build or modify organisms, and it is not a medical device (see [Limits and responsible use](#limits-and-responsible-use)).

---

## What is in the window

Open the program and your browser shows one page with ten tabs:

| Tab | What it does |
|---|---|
| **Identify** | Paste a DNA or RNA sequence (FASTA is fine). Your GPU compares it with 47,271 reference sequences (19,625 viral genomes and 27,646 bacterial 16S rRNA genes from NCBI RefSeq) and lists the closest matches with their full taxonomic lineage. Answers in about one second. |
| **DNA** | Draws an idealised B-DNA double helix from a real sequence (the phiX174 genome, or a sequence you type). |
| **Proteins** | Real structures from the Protein Data Bank (hemoglobin, lysozyme, green fluorescent protein, the SARS-CoV-2 spike). Space-filling or ribbon style, a still picture, a 60-frame orbit movie, or a fully interactive 3D view you can rotate with the mouse. |
| **Bacterium** | A schematic cutaway of a rod-shaped bacterium (envelope layers, DNA nucleoid, ribosomes, flagella and pili). It is an illustration, not a measured structure. |
| **Genome map** | The whole 4.6-million-base *E. coli* K-12 genome as a circle: candidate genes on both strands, GC content, GC skew and a ruler. |
| **Evolution** | Type a group name (for example *Betacoronavirus* or *Enterobacteriaceae*) and get a family tree of the matching sequences in the database, built from sequence distances. |
| **Assembly** | Shows how short sequencing reads are stitched back into a genome: simulate reads from a real small genome with a chosen coverage and error rate, or paste your own reads. |
| **Medicine** | A viewer for four well-studied drug targets (HIV-1 protease, SARS-CoV-2 main protease, influenza neuraminidase, the EGFR cancer switch). It shows the real 3D structure with the real bound molecule, lists the contacts and the binding pocket measured from the structure file, and shows published potency measurements. It shows what is *known*. It does not predict or suggest new drugs. |
| **Align** | Exact pairwise alignment of two sequences (score, CIGAR string and every substitution, insertion and deletion). |
| **Sequence tools** | Base counts, GC percent, reverse complement and a six-frame open-reading-frame finder with protein translation. |

A longer walk-through of every tab is in [docs/USER_GUIDE.md](docs/USER_GUIDE.md), and an explanation of how each part works is in [docs/HOW_IT_WORKS.md](docs/HOW_IT_WORKS.md).

## Two editions

Both editions give identical results. They differ only in which kind of graphics card they support.

| | **Vulkan edition** | **CUDA edition** |
|---|---|---|
| Graphics cards | NVIDIA, AMD and Intel cards with a Vulkan 1.1 driver | NVIDIA only, GeForce RTX 20-series / GTX 16-series or newer |
| Needs | Your normal, current graphics driver | A current NVIDIA driver (version 580 or newer) |
| Tested on | NVIDIA GeForce RTX 3060 | NVIDIA GeForce RTX 3060 |

> **Honest note on testing.** Both editions were developed and tested on one machine (Windows 11, NVIDIA RTX 3060). The Vulkan edition is designed to run on AMD and Intel graphics, but that has **not** been tested yet. If you try it, please open an issue and tell us what happened. If in doubt, pick the edition that matches your card, and use Vulkan if you have a non-NVIDIA card.

## Install (Windows 10/11, 64-bit)

You download **two** zip files from the [Releases](../../releases) page:

1. **One edition**: `MicroGenomeLab-Vulkan-win64.zip` *or* `MicroGenomeLab-CUDA-win64.zip`
2. **The data pack**: `MicroGenomeLab-data.zip` (about 46 MB; the reference database, structures and medicine facts. The same pack works with both editions.)

Then:

1. Unzip the edition anywhere you like (a normal folder in your user area is best, for example `Documents`).
2. Unzip the data pack **into the folder that contains `MicroGenomeLab.exe`**, so that a folder called `data` sits next to `bin` and `MicroGenomeLab.exe`.
3. Double-click `MicroGenomeLab.exe`. A console window opens and your browser opens the lab at `http://127.0.0.1:8770`. Close the console window to stop it.

If Windows SmartScreen warns about an unrecognised app, that is because the programs are not code-signed. Choose *More info*, then *Run anyway*. You can compare the file checksums listed on the Releases page with the files you downloaded.

The lab only listens on your own PC (127.0.0.1). It cannot be reached from your network.

## Folder layout

```
MicroGenomeLab/
  MicroGenomeLab.exe      the lab (starts the page in your browser)
  bin/                    the GPU and CPU programs for your edition
  data/                   the data pack (database, structures, medicine facts)
  out/                    pictures and movies you create (made on first use)
```

## Running from the Python source (optional)

The lab's window, web server and all of its Python analysis code are in the [`lab/`](lab/) folder under the MIT licence. To run it from source you still need the `bin/` and `data/` folders from a release:

```
pip install -r requirements.txt
set MGL_BIN=C:\path\to\MicroGenomeLab\bin
python lab\lab.py
```

The data folder is expected at `<parent of lab>\data`.

## Limits and responsible use

- **It is a screening and teaching tool, not a diagnosis.** *Identify* ranks the closest reference sequences by a fast sketch comparison. Close relatives and strains of one species can swap places. Never use it to make medical, legal or biosafety decisions.
- **The database is a fixed snapshot** (NCBI RefSeq viral genomes and the NCBI 16S rRNA set, downloaded in November 2026). It contains no sequence from after that date and not every organism on Earth.
- **The bacterium picture is schematic.** Envelope layers are drawn about twice as thick as real ones so that they can be seen.
- **Assembly is for small genomes** (viruses and phages) and shows the idea; it is not a replacement for a research assembler.
- **The Medicine tab shows known data only.** Potency numbers are published laboratory measurements from different experiments and are not directly comparable. It cannot tell you what to take, and it does not suggest new compounds.
- **Analysis only.** The tools are built to look at and understand sequences and structures. Do not use them to design, enhance or modify pathogens.

## Licence

- The Python source in [`lab/`](lab/) and the documents are under the [MIT licence](LICENSE).
- The ready-built programs in the release downloads (`bin/` and `MicroGenomeLab.exe`) are free to use and to share unchanged, under [LICENSE-BINARIES.txt](LICENSE-BINARIES.txt).
- The data keeps the licences of its sources. See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). In particular the medicine facts come from ChEMBL and are shared under CC BY-SA 3.0.

## Credits

Reference sequences: NCBI RefSeq and the NCBI 16S rRNA database. Structures: RCSB Protein Data Bank. Medicine facts: ChEMBL (EMBL-EBI). Full details in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

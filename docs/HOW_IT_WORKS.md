# How MicroGenome Lab works

This document explains what each part of the lab does, which methods it uses, and what the numbers mean. It describes the science and the architecture. The Python source of the lab is in [`lab/`](../lab/).

## The big picture

```
 your browser  <──HTTP, this PC only──>  MicroGenomeLab.exe  (small local web server, Python)
                                              │
              ┌───────────────────────────────┼─────────────────────────────┐
              │                               │                             │
        analysis in Python           GPU programs in bin/             data/ folder
   (sketches, trees, assembly,     (pictures, sequence search)     (reference database,
    structures, sequence tools)     + a CPU program (aligner)        structures, medicine facts)
```

- **The page** (`lab/static/index.html`, `viewer3d.js`) is plain HTML and JavaScript. The interactive 3D viewer draws with WebGL2 inside your browser.
- **The server** (`lab/lab.py`) listens on `127.0.0.1` only, so no other computer can reach it. It accepts a fixed list of actions, cleans every input (sequences are reduced to A, C, G, T, N; structure names must come from a list; numbers are clamped), starts helper programs with argument lists rather than shell strings, and runs one GPU job at a time.
- **The helper programs** in `bin/` are small stand-alone programs. Each reads its input from a few files in a working folder, does its job and writes a result file. The Python side prepares the inputs (for example the camera and the atoms) and turns the raw result into a PNG. The programs are:

  | Program | Job |
  |---|---|
  | `search` | compare a query fingerprint with every reference sequence on the GPU |
  | `mol` | draw a protein structure |
  | `helix` | draw the DNA double helix |
  | `cell` | draw the schematic bacterium |
  | `genome` | draw the circular genome map |
  | `align` | exact pairwise alignment (CPU) |

- **Two editions.** The Vulkan edition and the CUDA edition contain the same set of programs, built for the two GPU interfaces. They produce the same results (the pictures agree to within rounding of the last colour step; the search results are bit-identical).

## Identify: scaled MinHash containment on the GPU

**Fingerprints.** A sequence is cut into every overlapping 21-base word (a *k-mer*, k = 21). To make a reverse-complement of a sequence look the same as the sequence, each word is replaced by the smaller of itself and its reverse complement (*canonical* k-mers). Each word is turned into a 32-bit number by a hash function, and only the words whose hash falls in the lowest 1/50 of the range are kept. This is a *scaled MinHash* (also called FracMinHash): a random 2% sample of the words that is the same wherever the sequence is seen. A 30,000-base virus keeps roughly 600 numbers; a 1,500-base 16S gene keeps about 30.

**The database.** The same fingerprint was computed once for 47,271 reference sequences: 19,625 viral genomes (NCBI RefSeq viral) and 27,646 bacterial and archaeal 16S rRNA sequences (NCBI 16S RefSeq set). Together they hold 623,757,420 bases but only 12,345,067 fingerprint numbers (about 49 MB), stored sorted per sequence.

**The search.** For your query the lab computes its fingerprint (up to 4,096 numbers) and sends it to the GPU. The GPU assigns one thread to each reference sequence; each thread looks up every query number in its sequence's sorted list by binary search and counts the hits. All 47,271 comparisons run in parallel in well under a second.

**The score.** *Containment* = hits ÷ (query fingerprint size). If your sequence is a piece or a close relative of a reference, most of its numbers are found there, so containment is close to 1. Mutations break words (one substitution destroys 21 overlapping words) so containment falls smoothly with divergence, which is what makes it useful for finding relatives. An unrelated sequence scores near 0, and the lab calls a result "no confident match" below 0.10.

**The lineage.** Each reference's NCBI taxonomy ID is looked up in the NCBI taxonomy to produce its full lineage (for example Viruses → Riboviria → … → Betacoronavirus → SARS-CoV-2). About 79% of viral and 99% of 16S entries resolve to a complete lineage.

**Limits.** This is a screen. Sketches see shared words, not alignments: strains of one species can swap places, short or very ambiguous queries give weak answers, and anything not in the snapshot cannot be matched.

## DNA double helix

The coordinates are textbook B-DNA: 10.5 base pairs per turn (34.3° twist), 3.38 Å rise, backbone radius 10 Å, and the second strand's backbone offset by 126° (this gives the narrow minor groove on one side and the wide major groove on the other). For every base pair the lab places two backbone beads, a rod for each base (coloured by base) and hydrogen-bond dots in the middle.

The picture is made by **ray marching**: for each pixel the GPU follows a ray into the scene and steps along it by the distance to the nearest surface (a *signed distance field*), merging nearby shapes smoothly. The shading adds a key light with soft shadows, ambient occlusion (how enclosed a point is), a rim highlight, a procedural surface texture and a filmic tone-mapping curve. Every pixel is independent, which is why a GPU draws a 1280 × 720 picture in a fraction of a second.

## Proteins

The lab reads the atom lines (`ATOM` and `HETATM`) of a PDB file and drops water.

- **Space-filling** uses each atom's van der Waals radius (carbon 1.70 Å, nitrogen 1.55 Å, oxygen 1.52 Å, sulfur 1.80 Å and so on). The GPU intersects each ray exactly with every sphere (a quadratic equation per atom), finds the nearest, and shades it with a light, a hard shadow ray, six ambient-occlusion rays and a procedural texture bump.
- **Ribbon** passes a smooth Catmull-Rom spline through the C-alpha trace of every chain and draws it as a chain of spheres whose radius depends on the secondary structure recorded in the file (helix 1.75 Å, sheet 1.45 Å, loop 0.85 Å).
- **Colours.** A single-chain protein is coloured as a rainbow from its N-terminus to its C-terminus; multi-chain structures get one colour per chain.
- **Movie.** The orbit runs the same program 60 times with the camera turned by 6° each time.
- **Interactive 3D** does not use the GPU programs. The Python side sends sphere or ribbon data to the page, and `viewer3d.js` draws it with WebGL2 so that you can rotate in real time.

## Bacterium (schematic)

An illustrative cutaway built from signed distance fields: a capsule-shaped body with the top half removed, banded envelope layers on the cut face, cytoplasm, a folded DNA nucleoid, ribosome grains, three helical flagella and sixteen pili (drawn as thin capsules). All proportions are indicative, and the envelope is exaggerated. It is meant to teach how a cell is organised.

## Genome map

The lab reads the *E. coli* K-12 MG1655 genome (NCBI NC_000913.3, 4,641,652 bases) and:

1. finds **open reading frames** in all six frames: a start codon (ATG, GTG or TTG) followed in the same frame by no stop (TAA, TAG, TGA) for at least 150 codons, keeping the longest per stop;
2. computes **GC content** (fraction of G and C) in 3,600 equal windows;
3. computes **GC skew** (G − C)/(G + C) in the same windows and smooths it over ±20 windows (about 50 kb). The skew changes sign at the replication origin and terminus because the two strands are replicated differently.

The GPU then draws the rings in polar coordinates, one thread per pixel.

## Evolution

For a chosen taxonomic group the lab collects its sequences, computes the pairwise **Mash distance** from the fingerprints: with Jaccard index *j* of two fingerprints, D = −ln(2j / (1 + j)) / k, an estimate of substitutions per site. A **neighbour-joining** tree (Saitou and Nei, 1987) is built from the distance matrix. Branch lengths are estimated substitutions per site. Pairs that share no fingerprint number saturate at the maximum distance.

## Assembly

1. Simulation (optional): reads of 150 bases are drawn at random positions, on random strands, from a real genome, with substitution errors at the rate you choose.
2. Every read is cut into canonical 31-mers and the lab counts how often each occurs.
3. Words seen fewer times than a threshold derived from the coverage are discarded as probable errors.
4. The remaining words form a **de Bruijn graph** (two words are linked if they overlap by 30 bases). Paths that do not branch are merged into **unitigs** (contigs).
5. The report shows the contigs, their lengths and, for simulations, how much of the true genome's 31-mers the contigs recover. On phiX174, lambda phage and SARS-CoV-2 at 30× coverage and 0.3% error the genome comes back as one contig covering 99.3 to 99.9% of its words.

## Medicine tab

The tab reads a real structure of a drug target with a bound molecule from the Protein Data Bank (HIV-1 protease with a cyclic urea inhibitor, SARS-CoV-2 main protease with a peptide inhibitor, influenza neuraminidase, and EGFR kinase with erlotinib). It then:

- finds the **pocket**: residues with an atom within a set distance of any atom of the bound molecule;
- lists **contacts** and **possible hydrogen bonds** (N or O atom pairs closer than 3.5 Å between the molecule and the protein), with the measured distances;
- shows **potency measurements** (IC50 in nanomolar) from a snapshot of ChEMBL stored in `data/medpack/medicine_pack.json`, for compounds tested against the matching target, plus the recorded mechanism of action.

Everything shown is measured from the structure file or published. Nothing is predicted.

## Align

A **Needleman-Wunsch** global alignment with **affine gap penalties** solved with Gotoh's algorithm, so the score is exactly optimal. Scores: match +2, mismatch −3, a gap of length *L* costs 5 + 2·(L − 1). The program keeps three score tables and a traceback, then prints the score, a CIGAR string and each variant. Memory grows with the product of the two lengths, which is why the limit is 3,500 bases each.

## Sequence tools

Counts, GC percent, reverse complement and a six-frame open-reading-frame finder (start ATG, GTG or TTG; stop TAA, TAG or TGA; at least 50 amino acids) with translation through the standard genetic code. Candidates are not curated gene annotation.

## What is verified

- The GPU search returns identical hit counts for all 47,271 sequences in both editions.
- The pictures from the two editions agree (average difference below 0.03 out of 255 per colour channel).
- The aligner was checked against hand-worked cases and reference outputs.
- The assembler was run on three real genomes and rebuilt each as one contig.
- Both editions pass the same end-to-end test of every action of the lab on the development machine (Windows 11, NVIDIA RTX 3060).

## What is not yet verified

- Vulkan edition on AMD and Intel graphics cards.
- CUDA edition on cards other than the RTX 3060.
- Windows 10. The packaged lab was tested on the development machine with a minimal system path (no Python, compilers or SDKs on the path), but not yet on a second, separate PC.

# MicroGenome Lab: user guide

This guide walks through every tab. Start the program (`MicroGenomeLab.exe`), and your browser opens the lab. The footer at the bottom left shows how many sequences are indexed and which edition (Vulkan or CUDA) and graphics card are in use.

Everything you type stays on your PC.

## Starting, stopping and settings

- **Start:** double-click `MicroGenomeLab.exe`. A console window opens (it prints the address) and the browser opens `http://127.0.0.1:8770`.
- **Stop:** close the console window.
- **Different port:** run `MicroGenomeLab.exe 9000` from a command prompt to use port 9000 instead. Add `--no-browser` to stop it opening the browser by itself.
- **Pictures:** the pictures and movies you make are kept in the `out` folder, and the same request is answered instantly the second time. Delete the folder's contents any time to free space.
- **One GPU job at a time:** requests are queued, so clicking several buttons quickly is safe.

## Identify a sequence

1. Paste DNA or RNA into the box. A FASTA header line (`>name`) is allowed. Letters other than A, C, G, T/U are treated as unknown bases. You need at least 60 bases, and a few hundred or more give much better answers.
2. Press **Identify**. In about a second you get up to ten ranked matches.

For each match you see the reference name, whether it is a virus or a 16S rRNA gene, its length, how many of your sequence's fingerprint pieces it contains (**containment**, 0 to 1), and its taxonomic lineage from the NCBI taxonomy.

How to read it:

- **Containment near 1** means almost everything in your sequence is found in that reference. A 3% mutated copy of an *E. coli* 16S gene still scores high; an unrelated random sequence scores close to 0.
- If the best containment is below 0.10 the lab says there is **no confident match**.
- Close relatives and strains of one species can swap places. Treat the answer as a screen that tells you where to look, not as a diagnosis.

Four example buttons fill the box for you: a mutated *E. coli* 16S gene, a piece of the SARS-CoV-2 spike region, a piece of a hepatitis B genome, and random DNA that should match nothing.

## DNA

Draws an idealised B-DNA double helix (3.38 Å rise per base pair, 10.5 base pairs per turn, backbone radius 10 Å) from a sequence.

- Leave the sequence box empty to use the real phiX174 genome, starting at the base number you set in **Start base**; or paste your own sequence (at least 8 bases).
- **Base pairs:** how many base pairs to draw (8 to 36).
- **Angle** and **Elevation:** the camera position.
- The four bases are coloured differently, so you can read your sequence along the helix. A link opens the interactive 3D viewer.

The model is idealised textbook geometry. A real DNA molecule bends and wiggles.

## Proteins

Pick a structure and a view:

- **Structures:** lysozyme (1LYZ), hemoglobin (4HHB), green fluorescent protein (1EMA) and the SARS-CoV-2 spike (6VXX). The coordinates are the real, public Protein Data Bank files.
- **Ribbon** draws the backbone as a tube, thicker for helices and sheets. **Space-filling atoms** draws every atom as a sphere (carbon coloured by chain, nitrogen blue, oxygen red, sulfur yellow).
- **Render picture** makes a still image with shadows and ambient occlusion on your GPU. **Render orbit movie** makes a 60-frame full turn (about 15 to 25 seconds). If `ffmpeg` is installed on your PC you get an MP4, otherwise an animated GIF.
- **Interactive 3D** opens a live viewer. Drag to rotate, wheel to zoom, shift + drag to move, double-click to start or stop the spin.

## Bacterium

A schematic cutaway of a rod-shaped bacterium about 2 × 0.5 micrometres: layered envelope, folded DNA nucleoid, ribosome grains, flagella and pili. It is an illustration of how a cell is organised. It is not a measured structure, and the envelope layers are drawn about twice as thick as life so that you can see them.

## Genome map

Press **Draw the map** to see the whole 4.6-million-base *E. coli* K-12 genome (NCBI NC_000913.3) as a circle. From the outside in:

1. a ruler (a tick every 100,000 bases, a long tick every 500,000);
2. candidate genes on the forward strand;
3. candidate genes on the reverse strand;
4. GC content;
5. GC skew, smoothed over about 50,000 bases. Its sign change marks the replication origin and terminus.

Candidate genes are open reading frames of at least 150 amino acids found in all six frames. *E. coli* has about 4,300 annotated genes, so some candidates are false: it is a map of what the sequence alone suggests.

## Evolution (trees)

Type a group name, choose viral genomes or 16S rRNA, choose how many leaves (4 to 40) and press **Build tree**. Any name from the NCBI taxonomy works (for example *Betacoronavirus*, *Herpesviridae*, *Escherichia*, *Enterobacteriaceae*, *Bacillus*). You can also paste your own sequence to add it as one more leaf and see where it falls.

The tree is a neighbour-joining tree from distances between sequence fingerprints. It is reliable for relatives within a family or genus; very distant pairs saturate. 16S fingerprints hold only about 30 pieces each, so 16S trees are coarse.

## Assembly

Sequencing machines read short pieces (reads). Assembly puts the pieces back together.

- **Simulate and assemble:** choose a real genome (phiX174, hepatitis B, SARS-CoV-2 or lambda phage), a coverage (10 to 100×) and an error rate (0 to 3%). The lab makes random 150-base reads with errors, assembles them and reports how much of the true genome was rebuilt.
- **Assemble my reads:** paste FASTA or FASTQ reads (up to about 6 MB; at least 20 reads).

The assembler counts overlapping 31-base words, discards rare ones (probably errors) and joins the unambiguous paths into contigs. It is meant for small genomes.

## Medicine search

Four well-studied drug targets: HIV-1 protease, the SARS-CoV-2 main protease, influenza neuraminidase and the EGFR cancer switch. For each one the lab shows:

- the real 3D structure with the real bound molecule from the Protein Data Bank file;
- the residues that touch the bound molecule and the possible hydrogen bonds, measured from the structure file (a nitrogen or oxygen pair closer than 3.5 Å is listed);
- published potency measurements (IC50 values, lower is more potent) for compounds tested against the target, from a snapshot of ChEMBL included in the data pack, so it works offline.

This tab **shows what is known**. It does not predict or suggest new compounds, doses or safety. Finding a medicine needs laboratory and clinical testing. Different assays are not directly comparable.

## Align two sequences

Paste a reference and a query (up to 3,500 bases each) and press **Align**. The result is an exact global alignment with affine gaps: match +2, mismatch −3, a gap costs 5 for the first base and 2 for each further base. You get:

```
SCORE 20
CIGAR 11M1D4M          M = aligned base pair, I = extra in the query, D = missing from the query
SNP 5 A T              position 5 on the reference: A in the reference, T in the query
DEL 12 1               one base missing starting at reference position 12
```

Positions are 1-based on the reference. Many alignments can share the best score; the score is always the true optimum.

## Sequence tools

Paste a sequence and press **Analyse** to get its length, A/C/G/T/N counts, GC percent, the reverse complement (first 600 bases) and the longest candidate open reading frames (at least 50 amino acids, six frames) with their protein translation.

## Troubleshooting

| Message or symptom | What to do |
|---|---|
| "The data pack is missing" | Unzip `MicroGenomeLab-data.zip` into the folder that contains `MicroGenomeLab.exe`, so that `data` is next to `bin`. |
| "The GPU program ... produced no result" | Update your graphics driver. For the CUDA edition you need an NVIDIA card with a current driver; try the Vulkan edition if you are unsure. |
| Windows says the app is unrecognised | The programs are not code-signed. Choose *More info*, then *Run anyway*. |
| Port already in use | Another copy is running. Close it, or start with another port number. |
| Movies are GIFs, not MP4 | Install `ffmpeg` and put it on your PATH if you want MP4. |

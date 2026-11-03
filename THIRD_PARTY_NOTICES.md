# Third-party notices and data sources

MicroGenome Lab uses public data and open-source components. This file lists them and their licences. Please keep it when you redistribute the lab.

## Data in the data pack (`MicroGenomeLab-data.zip`)

| Data | Source | Terms |
|---|---|---|
| 19,625 viral genome sequences and their names | NCBI RefSeq viral genomes (`viral.1.1.genomic.fna`), https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/ | NCBI data are in the public domain in the United States. NCBI asks that you cite the source and notes that some submitters may claim rights; see https://www.ncbi.nlm.nih.gov/home/about/policies/ |
| 27,646 16S rRNA sequences | NCBI 16S ribosomal RNA RefSeq set (BLAST database `16S_ribosomal_RNA`), https://ftp.ncbi.nlm.nih.gov/blast/db/ | same as above |
| Taxonomic lineages | NCBI Taxonomy (as shipped with the BLAST taxonomy files) | same as above |
| Escherichia coli K-12 MG1655 genome (NC_000913.3) and the five small reference genomes in `data/genomes` (phiX174, lambda phage, SARS-CoV-2, hepatitis B genotype D, E. coli 16S) | NCBI GenBank/RefSeq | same as above |
| Protein structures 1LYZ, 4HHB, 1EMA, 6VXX, 6LU7, 1HVR, 2HU4, 1M17 | RCSB Protein Data Bank, https://www.rcsb.org | PDB archive data are available without restriction (CC0). Please cite the original structure papers; the citation of every entry is on its RCSB page. |
| Medicine facts (`data/medpack/medicine_pack.json`): measured potency values and mechanism rows for four targets | ChEMBL, EMBL-EBI, https://www.ebi.ac.uk/chembl/ (snapshot fetched November 2026) | **CC BY-SA 3.0.** If you share this file or a work derived from it you must give attribution and share it under the same licence. Cite: Zdrazil et al., *Nucleic Acids Research* 2024, "The ChEMBL Database in 2023". |

The reference fingerprint index (`data/index`) is derived from the sequence data above. The lab does not redistribute the original FASTA files except for the six genomes listed.

## Components inside the programs

| Component | Used for | Licence |
|---|---|---|
| Python 3 | the lab's server and analysis code | PSF licence |
| NumPy | array computation | BSD 3-clause |
| Pillow | writing PNG/GIF pictures | HPND (Pillow) licence |
| PyInstaller | packaging `MicroGenomeLab.exe` | GPL-2.0 with the bootloader exception, which permits distributing the packaged program under any licence |
| Khronos Vulkan headers and loader interface | the Vulkan edition | Apache-2.0. The programs call the Vulkan loader (`vulkan-1.dll`) that is installed with your graphics driver; it is not shipped. |
| Microsoft Visual C++ runtime (statically linked) | all programs | Microsoft redistributable terms |
| NVIDIA CUDA runtime (statically linked) | the CUDA edition | NVIDIA CUDA Toolkit EULA (redistributable runtime). NVIDIA, CUDA and GeForce are trademarks of NVIDIA Corporation. |
| MinGW-w64 C runtime (statically linked) | the aligner | permissive licences of the MinGW-w64 project, https://www.mingw-w64.org |

In the interactive 3D view the page uses your browser's built-in WebGL2. No third-party JavaScript library is loaded.

## Methods

- Scaled MinHash / containment: Broder 1997; Ondov et al., "Mash", *Genome Biology* 2016; Irber et al., "Lightweight compositional analysis of metagenomes with FracMinHash and minimizer sketches", 2022.
- Neighbour-joining: Saitou and Nei, *Molecular Biology and Evolution* 1987.
- Affine-gap alignment: Needleman and Wunsch 1970; Gotoh 1982.
- Assembly: de Bruijn graphs for sequence assembly (Pevzner, Tang and Waterman 2001).

All trademarks belong to their owners. Use of a name does not imply endorsement.

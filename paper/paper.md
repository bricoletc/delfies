---
title: 'delfies: a Python package for the detection of DNA elimination breakpoints with neo-telomere addition'
tags:
  - Python
  - Bioinformatics
  - Genomics
  - Programmed DNA Elimination
  - Soma/germline differentiation
authors:
  - name: Brice Letcher
    orcid: 0000-0002-8921-6005
    affiliation: 1 # (Multiple affiliations must be quoted)
affiliations:
  - name: Laboratory of Biology and Modeling of the Cell, Ecole Normale Supérieure de Lyon, CNRS UMR 5239, Inserm U1293, University Claude Bernard Lyon 1, Lyon, France
    index: 1
bibliography: paper.bib
---

# Summary

In multicellular organisms, cells generally carry an identical genome,
transmitted through successive cell divisions from the founding zygote. This is not the
case in species that undergo Programmed DNA Elimination (PDE), the systematic
destruction of portions of the genome in somatic cells during early development.
In these species, germline cells maintain an intact genome, while somatic cells 
maintain a reduced genome.

PDE was first documented in 1887 in a parasitic nematode [@Boveri1887], and
since then various forms of PDE have been found in a wide variety of organisms,
including birds, fish, insects, mammals, crustaceans, other nematodes
and ciliates [@Dedukh2021; @Drotos2022]. Some species eliminate entire 
chromosomes (birds, fish, insects, mammals) while others eliminate portions of chromosomes, 
with or without changes in chromosome number (copepod crustaceans, nematodes, ciliates).

While of considerable interest, the function of PDE remains largely unknown.
One of the best-studied clade so far is ciliates, unicellular eukaryotes that
maintain distinct 'germline' and 'somatic' nuclei within the same cellular
membrane. Ciliates eliminate both small sequences that are spliced out of the
genome, called internal eliminated sequences (IESs), and large fragments of
chromosomes, including at the subtelomeres. In ciliates, the latter has
received less attention, both genomically and functionally.

Copepods and nematodes also eliminate large fragments of DNA, including at the
subtelomeres. The main commonality so far between these three clades is the
elimination of all telomeric sequences during PDE, followed by the addition of
new telomeres at the extermities of newly-maintained mini-chromosomes in
somatic cells.

Here, we present a tool called `delfies` to systematically detect these sites
of chromosome breakage and neo-telomere addition. `delfies` enables rapidly and
comprehensively mapping the locations of elimination breakpoints in all
species in which this form of PDE occurs.


# Statement of need

Several other tools for the detection of DNA elimination breakpoints have been
developed and validated in the context of ciliates: `parTIES` [@parties:2015],
`SIGAR` [@sigar:2020], `ADFinder` [@adfinder:2020] and `bleTIES`
[@bleties:2021]. Of these, `parTIES`, `ADFinder` and `SIGAR` allow the
detection of IESs only, not sites of chromosome breakage with de-novo telomere
addition, and were primarily designed for short-read sequencing data. `bleTIES`
was designed to detect and reconstruct IESs in ciliates in the context of
long-read sequencing data, and also includes a module for detecting chomosome
breakage sites with telomere addition, called MILTEL [@bleties:2021].

`delfies` presents several new features compared to MILTEL. Both tools output
the locations of breakpoints, in standard bioinformatics formats: MILTEL in a
GFF3-formatted file, `delfies` in a BED-formatted file. While MILTEL expresses
breakpoints in isolation, `delfies` can merge multiple breakpoints occurring in
close proximity, in a user-configurable way. This allows for simplified output,
and for more or less sharply-defined breakpoints. `delfies` also outputs the
'strand' of breakpoints in the appropriate BED column, enabling easily
classifying the genome in 'retained' and 'eliminated' compartments (see
repository docs for how). 

`delfies` additionally extracts and outputs the sequences around the
breakpoints in a Fasta-formatted file, which MILTEL does not. This enables
searching for motifs specifying breakpoints, e.g. using MEME [@Bailey2015].

MILTEL was designed to detect telomere addition in the reads only. `delfies`
additionally detects breakpoints in which telomeres have been assembled in the
genome and reads containing non-telomeric sequence from the 'complete' genome.
These breakpoints can be used to re-assemble the genome past the 'reduced
genome' (documented in the software repository).

In practical terms, `delfies` has a highly configurable command-line interface,
enabling specifying how much to filter read alignments, which regions of the
genome to analyse and the types of breakpoints to look for. On a nematode
genome of size 240Mbp sequenced at 85X average coverage with PacBio HiFi data,
`delfies` finds all breakpoints in less than 2 minutes, using a single thread.
For further speed, `delfies` also supports multi-threading.

`delfies` has already been used to successfully characterise the breakpoints,
motifs, and retained/eliminated genomes of nematode genera in the family
*Rhabditidae*, supporting two upcoming publications (Letcher *et al.* and
Stevens *et al.*, both in preparation). As more and more genomic data from
species undergoing Programmed DNA Elimination become available, we anticipate
this tool can be of broad use to the research community as a whole.

# Acknowledgements

We acknowledge the many interactions with Lewis Stevens and Pablo Manuel Gonzalez de la Rosa 
at the Wellcome Sanger Institute, which helped foster the development of `delfies`.

# References

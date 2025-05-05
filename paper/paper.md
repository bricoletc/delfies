---
title: 'delfies: a Python package for the detection of DNA breakpoints with neo-telomere addition'
tags:
  - Python
  - Bioinformatics
  - Genomics
  - Programmed DNA Elimination
  - Soma/germline differentiation
authors:
  - name: Brice Letcher
    corresponding: true
    orcid: 0000-0002-8921-6005
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Marie Delattre
    corresponding: true
    orcid: 0000-0003-1640-0300
    affiliation: 1
affiliations:
  - name: Laboratory of Biology and Modelling of the Cell, Ecole Normale Supérieure de Lyon, CNRS UMR 5239, Inserm U1293, University Claude Bernard Lyon 1, Lyon, France
    index: 1
bibliography: paper.bib
---

# Summary

In multicellular organisms, all cells generally carry an identical genome,
faithfully transmitted through cell divisions from the founding zygote. This is not the
case in species that undergo Programmed DNA Elimination (PDE), the systematic
destruction of portions of the genome in somatic cells during early development.
In these species, somatic cells carry a reduced genome, while germline cells maintain an 
intact genome.

PDE was first documented in 1887 in a parasitic nematode [@Boveri1887], and
since then various forms of PDE have been found in a wide variety of organisms
including birds, fish, insects, mammals, crustaceans, other nematodes and
ciliates [@Dedukh2021; @Drotos2022]. Some species eliminate entire chromosomes
(birds, fish, insects, mammals) while others eliminate portions of chromosomes,
with or without changes in chromosome number (copepod crustaceans, nematodes,
ciliates).

In species that eliminate portions of chromosomes, two main types of
elimination have been documented. The first is the elimination of small
sequences (~100s of bp) called 'IESs', by a splicing process: a double-strand
break is produced at each IES extremity, the IES is excised, and the two
extremities are rejoined. This form has so far been documented in ciliates
only. The second type is the elimination of large fragments of chromosomes (up
to >1Mbp): a single double-strand break is produced, one side is eliminated,
and telomeres on the retained side allow the new 'mini-chromosome' to be
maintained in the soma. This form occurs in ciliates [@Yu1991], nematodes
[@GonzalezdelaRosa2020; @Rey2023], and probably also in copepods
[@Beermann1977]. While IES elimination in ciliates has been well-characterised
genomically and functionally, chromosome fragmentation with neo-telomere
addition has not.

Here, we present a tool called `delfies` to systematically detect sites of
chromosome breakage and neo-telomere addition. `delfies` enables rapidly and
comprehensively mapping the locations of elimination breakpoints in all species
in which this form of DNA elimination occurs.

# Statement of need

Several other tools for the detection of DNA elimination breakpoints have been
developed and tested, all in the context of ciliates: `parTIES` [@parties:2015],
`SIGAR` [@sigar:2020], `ADFinder` [@adfinder:2020] and `bleTIES`
[@bleties:2021]. Of these, `parTIES`, `ADFinder` and `SIGAR` allow the
detection of IESs only, not sites of chromosome breakage with neo-telomere
addition, and were primarily designed for short-read sequencing data. `bleTIES`
was designed to detect and reconstruct IESs in the context of
long-read sequencing data, and also includes a module for detecting chromosome
breakage sites with telomere addition called MILTEL [@bleties:2021].

`delfies` was developed when studying PDE in nematodes, and presents several
new features compared to MILTEL. Both tools output the locations of breakpoints
in standard bioinformatics formats: MILTEL in a GFF3-formatted file, `delfies`
in a BED-formatted file. While MILTEL expresses each putative breakpoint in
isolation, `delfies` can merge multiple breakpoints occurring in close
proximity in a user-configurable way. This allows for directly detecting more
or less sharply-defined breakpoints, a feature that is known to vary in both
ciliates and nematodes [@Betermier2023; @GonzalezdelaRosa2020;
@Dockendorff2022; @Estrem2023]. `delfies` also outputs the strand of
breakpoints in the appropriate BED column, enabling subsequently classifying
the genome into 'retained' and 'eliminated' compartments (details in the
software repository). 

`delfies` also explicitly models and outputs two types of breakpoints: in the
first, the assembled genome is 'complete' and reads from the 'reduced' genome
contain telomeres after the breakpoint. In the second, the assembled genome is
'reduced' and contains telomeres, and reads from the 'complete' genome contain
unique non-telomeric sequence. These two types can be treated separately. For
example, in the case of a reduced assembled genome, reads coming from the
reduced genome can be specifically depleted at breakpoints and the new read-set
used to assemble the complete genome (instructions provided in the software
repository).

In addition to breakpoint locations, `delfies` extracts and outputs the
sequences around the breakpoints in a Fasta-formatted file. This enables
searching for motifs specifying breakpoints, e.g. using MEME [@Bailey2015].

In practical terms, `delfies` has a highly configurable command-line interface,
enabling specifying how much to filter read alignments, which regions of the
genome to analyse and the types of breakpoints to look for. On a nematode
genome of size 240Mbp sequenced at 85X average coverage with PacBio HiFi data,
`delfies` finds all breakpoints in less than 2 minutes, using a single thread.
For further speed, `delfies` also supports multi-threading.

`delfies` has already been used to successfully characterise the breakpoints,
motifs, and retained/eliminated genomes of several nematode genera in the
family *Rhabditidae*, supporting two upcoming publications (Letcher *et al.*
and Stevens *et al.*, in preparation). For testing purposes, the author has
prepared a subset of publicly-available data from the nematode *Oscheius
onirici*, whose elimination breakpoint motif has been previously described
[@Estrem2023]. The data are available on Zenodo [@ZenodoData2024], and consist
of a small genome region containing a single elimination breakpoint and alignment
files for reads sequenced using three distinct technologies: Illumina NovaSeq,
Oxford Nanopore Technologies PromethION and Pacific Biosciences Sequel II. The
reads span a range of average lengths (151bp to 11.9kbp) and per-base qualities
(Q11 to Q28). `delfies` recovered a single, identical breakpoint across all
three datasets.

We anticipate this tool can be of broad use to researchers studying Programmed
DNA Elimination, to characterise species known to eliminate but also to
discover or screen for elimination in new species. This is especially relevant
as new long-read and high-coverage sequencing data (of both germline and
somatic cells) of eukaryotic species become increasingly available [@DTOL2022;
@EBGP2022]. `delfies` may also be useful in other fields of research in which
modified chromosomes with neo-telomeres are formed and maintained, such as
cancer biology.

# Acknowledgements

We acknowledge the many interactions with Lewis Stevens and Pablo Manuel
Gonzalez de la Rosa at the Wellcome Sanger Institute and Marie Delattre at the
Laboratory of Biology and Modelling of the Cell, which helped foster the
development of `delfies`.

This work was supported by a grant from the Agence Nationale de la Rercherche: ANR-22-
CE12-0027.

# References

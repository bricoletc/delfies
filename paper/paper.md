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

In multicellular organisms, all cells generally carry an identical genome,
transmitted through successive cell divisions from the founding zygote. This is not the
case in species that undergo Programmed DNA Elimination (PDE), the systematic
destruction of portions of the genome in somatic cells during early development.
In these species, germline cells maintain an intact genome. 

While PDE was first documented in 1887 in a parasitic nematode [@Boveri1887],
since then various forms of PDE have been found in a wide variety of organisms,
including birds, fish, insects, mammals, crustaceans, non-parasitic nematodes
and ciliates [@Dedukh2021; @Droto2022]. In particular, some species eliminate entire 
chromosomes (birds, fish, insects, mammals) while others eliminate portions of chromosomes, 
with or without changes in chromosome number (copepod crustaceans, nematodes, ciliates).

While of considerable interest, the function of PDE remains largely unknown. 
One of the best-studied clade so far is ciliates, unicellular eukaryotes which maintain distinct 
'germline' and 'somatic' nuclei within the same cellular membrane. Ciliates eliminate both 
small sequences that are spliced out of the genome, called internal eliminated sequences (IESs), 
and large fragments of chromosomes, including at the subtelomeres. The latter has 
been less well-studied functionally. 

Copepods and nematodes also eliminate large fragments of DNA including at the subtelomeres. 
A commonality between the three clades is the elimination of all telomeric sequences during PDE, 
followed by the addition of new telomeres at the extermities of the new maintained mini-chromosomes 
in somatic cells.

Here, we present a tool to systematically detect these sites of chromosome breakage and neo-telomere 
addition, called `delfies`. We anticipate `delfies` will help researchers easily map the locations 
of elimination breakpoints, enabling them to detect sequence motifs that occur at breakpoints, delineate 
and characterise the eliminated and retained genomes, and possibly produce more contiguous assemblies.

# Statement of need

Several other tools for the detection of DNA elimination breakpoints have been developed, 
but mainly for the characterisation of IESs in ciliates. These are `parTIES` [@parties:2015], 
`SIGAR` [@sigar:2020], `ADFinder` [@adfinder:2020] and `bleTIES` [@bleties:2021].


<!-- looks for germline-specific reads mapping to somatically-assembled genomic sequences, in ciliates. -->

# Acknowledgements

We acknowledge the many interactions with Lewis Stevens and Pablo Manuel Gonzalez de la Rosa 
at the Wellcome Sanger Institute, which helped foster the development of `delfies`.

# References

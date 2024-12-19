<p align="center">
<img height="175" src="img/delfies_logo.png">
<h1 align="center"></h1>
</p>


[![PyPI](https://img.shields.io/pypi/v/delfies)](https://pypi.org/project/delfies/)
[![codecov](https://codecov.io/github/bricoletc/delfies/graph/badge.svg?token=7GP56CS6NU)](https://codecov.io/github/bricoletc/delfies)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![JOSS paper status](https://joss.theoj.org/papers/e9a54b5b34327be54050d6796cf9a31b/status.svg)](https://joss.theoj.org/papers/e9a54b5b34327be54050d6796cf9a31b)

`delfies` is a tool for the detection of DNA breakpoints with de-novo telomere addition.

It identifies genomic locations where double-strand breaks have occurred followed by telomere addition.
It was initially designed and validated for studying the process of Programmed DNA Elimination
in [nematodes][Meso_paper], but should work for other clades and applications too.

# <a name="started"></a> Getting started

`delfies` takes as input a genome fasta (gzipped supported) and an indexed SAM/BAM of 
sequencing reads aligned to the genome.

```sh
delfies --help
samtools index <aligned_reads>.bam
delfies <genome>.fa.gz <aligned_reads>.bam <output_dir>
cat <output_dir>/breakpoint_locations.bed
```

For how to obtain a suitable SAM/BAM, see [input data](#input_data), and for 
downloading a real genome and BAMs for a test run of `delfies`, see [test run](#test_run).

# Table of Contents

- [Installation](#installation)
- [Input data](#input_data)
- [Test run with real data](#test_run)
- [User Manual](#manual)
    - [CLI options](#CLI)
    - [Outputs](#outputs)
    - [Validation](#validation)
    - [Applications](#applications)
    - [Detailed documentation](#detailed_docs)
- [Contributing](#contributing)

# Installation
Using `pip` (or equivalent - poetry, etc.): 
```sh
# Install latest release from PyPI
pip install delfies

# Or install a specific release from PyPI:
pip install delfies==0.8.0

# Or clone and install tip of main
git clone https://github.com/bricoletc/delfies/
pip install ./delfies
```

# <a name="input_data"></a> Input data

## <a name="seq_tech"></a> Sequencing technologies

`delfies` is designed to work with both Illumina short reads and ONT or PacBio
long reads. Long reads are better for finding breakpoints in more repetitive
regions of the genome. A high fraction of sequenced bases with a quality \>Q20
is desirable (e.g. \>70%). I found `delfies` worked on recent data from all three
sequencing technologies: see [test run below](#test_run).

## Aligners

To produce a SAM/BAM with which you can find breakpoints, you need to use a read
aligner that reports soft clips (parts of a reads that are not aligned to the
reference). Both `bowtie2` (in `--local` mode) and `minimap2` (by default) do this. 
Use `minimap2` for long reads (>300bp), with the appropriate preset (e.g. `-x map-ont` 
for Nanopore data).

# <a name="test_run"></a> Test run with real data

I provide a processed subset of publicly-available data here: 
https://doi.org/10.5281/zenodo.14101797. 

The data consist of a 2kbp region of the assembled genome of *Oscheius onirici*
and three alignment BAMs from sequencing data produced using Illumina, ONT and
PacBio. The data were aligned to the 2kbp region using `minimap2`. See the
Zenodo link for details on the sequencing data (read lengths, error rates) and
public links to the raw data.

You can run `delfies` on the inputs in this archive to make sure it is properly 
installed and produces the expected outputs:

```sh
wget https://zenodo.org/records/14282333/files/delfies_zenodo_test_data.tar.gz
tar xf delfies_zenodo_test_data.tar.gz
# Run delfies; for example, having defined genome, bam and odirname variables:
delfies --threads 16 \
    --telo_forward_seq TTAGGC \
    --breakpoint_type all \
    --min_mapq 20 \
    --min_supporting_reads 6 \
    ${genome} ${bam} ${odirname}
# Compare with the expected outputs:
find delfies_zenodo_test_data -name "*breakpoint_locations.bed" | xargs cat
```

# <a name="manual"></a> User Manual

## <a name="CLI"></a> CLI options

```sh
delfies --help
```

* Do use the `--threads` option if you have multiple cores/CPUs available.
* [Breakpoints]
   * There are two types of breakpoints: see [detailed docs][detailed_docs].
   * Nearby breakpoints can be clustered together to account for variability in breakpoint location (`--clustering_threshold`).
* [Region selection]: You can select a specific region to focus on, specified as a string or as a BED file.
* [Telomeres] 
    * Specify the telomere sequence for your organism using `--telo_forward_seq`. 
      If you're unsure, I recommend the tool [telomeric-identifier](https://github.com/tolkit/telomeric-identifier) for finding out.
* [Aligned reads]
    * To analyse confidently-aligned reads only, you can filter reads by MAPQ (`--min_mapq`) and by bitwise flag (`--read_filter_flag`).
    * You can tolerate more or less mutations in the assembly telomeres (and in the sequencing reads) using `--telo_max_edit_distance` and `--telo_array_size`.

## Outputs

The two main outputs of `delfies` are:

- `breakpoint_locations.bed`: a BED-formatted file containing the location of identified 
   elimination breakpoints. 
- `breakpoint_sequences.fasta`: a FASTA-formatted file containing the sequences 
   of identified elimination breakpoints

## <a name="validation"></a> Validating breakpoints

**I highly recommend visualising your results**! E.g., by loading your input
fasta and BAM and output `delfies`' output `breakpoint_locations.bed` in
[IGV](https://github.com/igvteam/igv).

Confident/true breakpoints will typically have:

- Good read support. Note that breakpoints are ordered by read support in the `delfies` 
  output file `breakpoint_locations.bed`, and you can require a minimum number of 
  supporting reads using the CLI option `--min_supporting_reads`.
- A difference in read coverage before and after the breakpoint. The nature of this 
  difference depends on the ratio between cells with and without the breakpoint. As
  an example, in organisms that eliminate parts of their genome in the soma, if
  most sequenced cells are from the soma, expect more reads before the breakpoint 
  than after it ('before' and 'after' defined relative to the reported breakpoint strand).

Ultimately though, only biological experiments can truly validate identified breakpoints.

## Applications

* The fasta output enables looking for sequence motifs that occur at breakpoints, e.g. using [MEME](https://meme-suite.org/meme/).
* The BED output enables classifying a genome into retained and eliminated regions. 
  The 'strand' of breakpoints is especially useful for this: see [detailed docs][detailed_docs].
* The BED output also enables assembling past somatic telomeres: for how to do this, see [detailed docs][detailed_docs].

## <a name="detailed_docs"></a> Detailed documentation

For more details on `delfies`, including outputs and applications, see [detailed_docs][detailed_docs].

# Contributing

Contributions always welcome! 

Please see [CONTRIBUTING.md](CONTRIBUTING.md) for how (reporting issues, requesting
features, contributing code). This document includes instructions on how to run 
`delfies`' unit and functional tests.

[detailed_docs]: docs/detailed_manual.md
[Meso_paper]: https://doi.org/10.1016/j.cub.2023.07.058

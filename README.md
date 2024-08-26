# Delfies

`delfies` is a tool for the detection of DNA breakpoints with de-novo telomere addition.

It identifies genomic locations where double-strand breaks have occurred followed by telomere addition.
It was initially designed and validated for studying the process of Programmed DNA Elimination
in [nematodes](https://doi.org/10.1016/j.cub.2023.07.058), but should work for other clades and applications too.

# <a name="started"></a> Getting started

`delfies` takes as input a genome fasta (gzipped supported) and an indexed SAM/BAM of 
sequencing reads aligned to the genome.

```sh
delfies --help
samtools index <aligned_reads>.bam
delfies <genome>.fa.gz <aligned_reads>.bam <output_dir>
cat <output_dir>/breakpoint_locations.bed
```

# <a name="manual"></a> User Manual

## Installation
Using `pip` (or equivalent - poetry, etc.): 
```sh
# Download and install a specific release
DELFIES_VERSION=0.4.0
wget "https://github.com/bricoletc/delfies/archive/refs/tags/${DELFIES_VERSION}.tar.gz"
tar -xf "delfies-${DELFIES_VERSION}.tar.gz
pip install ./delfies-"${DELFIES_VERSION}"/

# OR clone and install tip of main
git clone https://github.com/bricoletc/delfies/
pip install ./delfies
```

## CLI options

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
    * You can tolerate more or less mutations in the telomere sequences (and in the reads) using `--telo_max_edit_distance` and `--telo_array_size`.

## Outputs

The two main outputs of `delfies` are:

- `breakpoint_locations.bed`: a BED-formatted file containing the location of identified 
   elimination breakpoints.
- `breakpoint_sequences.fasta`: a FASTA-formatted file containing the sequences 
   of identified elimination breakpoints

For more details on outputs, see [detailed docs][detailed_docs].

## Applications

* The fasta output enables looking for sequence motifs that occur at breakpoints, e.g. using [MEME](https://meme-suite.org/meme/).
* The BED output enables classifying a genome into retained and eliminated regions. 
  The 'strand' of breakpoints is especially useful for this: see [detailed docs][detailed_docs].
* The BED output also enables assembling past somatic telomeres: for how to do this, see [detailed docs][detailed_docs].

## Visualising your results

**I highly recommend visualising your results**!
E.g., by loading your input fasta and BAM and output `delfies`' output `breakpoint_locations.bed` in [IGV](https://github.com/igvteam/igv).

[detailed_docs]: docs/detailed_manual.md

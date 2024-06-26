# Delfies

Delfies is a tool for the detection of DNA Elimination breakpoints 

# Getting started

## Installation
Using `pip` (or equivalent - poetry, etc.): 
```sh
git clone https://github.com/bricoletc/delfies/
pip install ./delfies
```

## Usage

`delfies` takes as input a genome fasta (gzipped supported) and a BAM of sequencing reads 
aligned to the genome. 

For the full list of options, see:

```sh
delfies --help
```

## Outputs

The main outputs of `delfies` are:

- `breakpoint_maxima.bed`: a BED-formatted file containing the location of identified 
   elimination breakpoints. 
   - The location itself is provided as an interval of size one; however, 
     all putative breakpoint positions are merged into a larger interval 
     that is provided in the name column (column 4), as `maximal_focus_window: <start>-<end>`.
    - The strand column (column 6) is a '+' if the telomere-containing reads (mainly) occur 
      3' of the identified breakpoint, and '-' if they occur 5' of the identified breakpoint. 
      Thus a '+' suggests the eliminated DNA is 3' of the breakpoint, and vice-versa.
- `breakpoint_sequences.fasta`: a FASTA-formatted file containing the sequences 
   of identified elimination breakpoints. The header field contains:
   - A sequence ID as `<detection_mode>_<breakpoint_type>_<chrom>`

     `detection_mode` is currently always S2G, so ignore this (TODO: add G2S)

     `breakpoint_type` is '3prime' for a '+' breakpoint, and '5prime' for a '-' breakpoint (see BED file above)

     `chrom` is the contig/scaffold/chromosome name
    - Some additional information is provided, e.g. the position of the breakpoint 
      and the number of reads supporting the breakpoint ('num_telo_containing_softclips')

**Importantly, visualise the results yourself.** 
E.g., by loading the input fasta and BAM and output `breakpoint_maxima.bed` in [IGV](https://github.com/igvteam/igv).

# What's next

## Tool constraints

* Currently `delfies` looks for telomeric units repeated at least `--telo_array_size` times 
  at the beginning of the soft-clipped region of a read, with no mismatches in the telomeric 
  repeat unit.

## Detection mode 
* Add 'G2S' detection mode

## Benchmark

* Possibly benchmark against one or more tools below.

### Other tools

| Tool name   | Paper | Code |
| ----------- | ----- | ---- |
| ParTIES     | https://doi.org/10.1093/bioinformatics/btv691  | https://github.com/oarnaiz/ParTIES |
| SIGAR       | https://doi.org/10.1093/gbe/evaa147 | https://github.com/yifeng-evo/SIGAR |
| ADFinder    | https://doi.org/10.1093/bioinformatics/btaa226 | https://github.com/weibozheng/ADFinder |
| BleTIES     | https://doi.org/10.1093/bioinformatics/btab613 | https://github.com/Swart-lab/bleties |

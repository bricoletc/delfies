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

Do use the `--threads` option if you have multiple cores/CPUs available.

For the full list of options, see:

```sh
delfies --help
```

## Outputs

### Terminology: breakpoint types and strandedness

First, we provide two basic definitions:

- Breakpoints can be of two `breakpoint_type`s:
  - S2G: telomere-containing softclipped-reads align to a location in the genome. 
  - G2S: non-telomere-containing softclipped-reads align to a location in the genome 
    that contains telomeres
  
  These two types both describe elimination breakpoints at which telomeres have been 
  added to the retained fragments. In the case of `S2G` breakpoints, the assembled 
  genome is the unbroken genome, and breakpoint-supporting reads come from cells 
  with a broken genome. In the case of `G2S` breakpoints, the assembled genome is 
  the reduced genome (with telomeres), and breakpoint-supporting reads come from cells 
  with an unbroken genome

- Breakpoint strand. The strand is defined as '+', or also called '3prime', if 
  the softclips on reads occur 3' of the assembled genome, and '-' or '3prime', if 
  they occur 5' of the assembled genome. For both `S2G` and `G2S` breakpoints, '+' suggests the 
  eliminated genome occurs 3' of the identified breakpoint, and vice-versa.

### Files 
The main outputs of `delfies` are:

- `breakpoint_locations.bed`: a BED-formatted file containing the location of identified 
   elimination breakpoints. 
   - The location itself is provided as an interval of size one, and corresponds to the 
     first base past the putative breakpoint. 
   - The name column (column 4) records the breakpoint type, as defined above. 
     It also stores larger window, in the format `maximal_focus_window: <start>-<end>`,
     containing several adjacent putative breakpoint positions that have been merged together.
   - The score column (column 5) stores the number of sequencing reads that 
     support the sharply-defined identified breakpoint (=one-base interval).
    - The strand column (column 6) specifies the likely direction of eliminated DNA, 
      as defined above.
      
- `breakpoint_sequences.fasta`: a FASTA-formatted file containing the sequences 
   of identified elimination breakpoints. The header field contains:
   - A sequence ID as `<breakpoint_type>_<breakpoint_direction>_<chrom>`

     `breakpoint_type`: 'S2G' or 'G2S', as defined above

     `breakpoint_type`: '5prime' or '3prime', as defined above

     `chrom`: the contig/scaffold/chromosome name

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
* G2S mode: add germline sequence sequence reconstruction at breakpoints
* S2G mode: remove breakpoints that occur at locations of the genome containing telomeres?

## Benchmark

* Possibly benchmark against one or more tools below.

### Other tools

| Tool name   | Paper | Code |
| ----------- | ----- | ---- |
| ParTIES     | https://doi.org/10.1093/bioinformatics/btv691  | https://github.com/oarnaiz/ParTIES |
| SIGAR       | https://doi.org/10.1093/gbe/evaa147 | https://github.com/yifeng-evo/SIGAR |
| ADFinder    | https://doi.org/10.1093/bioinformatics/btaa226 | https://github.com/weibozheng/ADFinder |
| BleTIES     | https://doi.org/10.1093/bioinformatics/btab613 | https://github.com/Swart-lab/bleties |

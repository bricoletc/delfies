# Delfies

Delfies is a tool for the detection of DNA Elimination breakpoints 

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

## Tool constraints

* Currently `delfies` looks for telomeric units repeated at least `--telo_array_size` times 
  at the beginning of the soft-clipped region of a read, with no mismatches in the telomeric 
  repeat unit.

# Other tools

| Tool name   | Paper | Code |
| ----------- | ----- | ---- |
| ParTIES     | https://doi.org/10.1093/bioinformatics/btv691  | https://github.com/oarnaiz/ParTIES |
| SIGAR       | https://doi.org/10.1093/gbe/evaa147 | https://github.com/yifeng-evo/SIGAR |
| ADFinder    | https://doi.org/10.1093/bioinformatics/btaa226 | https://github.com/weibozheng/ADFinder |
| BleTIES     | https://doi.org/10.1093/bioinformatics/btab613 | https://github.com/Swart-lab/bleties |

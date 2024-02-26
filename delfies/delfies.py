import itertools as it
import multiprocessing as mp
import os
from glob import glob
from pathlib import Path

import click
from pybedtools import BedTool
from pysam import AlignmentFile

from delfies import __version__
from delfies.seq_utils import (
    REGION_CLICK_HELP,
    TELOMERE_SEQS,
    Orientation,
    rev_comp
)
from delfies.breakpoint_foci import find_breakpoint_foci_row_based


@click.command()
@click.argument("bam_fname", type=click.Path(exists=True))
@click.argument("ofname")
@click.option("--seq_region", type=str, help=REGION_CLICK_HELP)
@click.option(
    "--bed",
    type=click.Path(exists=True),
    help="Path to bed of regions to analyse. Overrides 'seq_region'",
)
@click.option(
    "--telomere_forward_seq",
    type=str,
    default=TELOMERE_SEQS["Nematoda"][Orientation.forward],
    help="The telomere sequence used by your organism! Please make sure this is provided in 'forward' orientation (i.e. 5'->3')",
    show_default=True,
)
@click.option(
    "--telo_array_size",
    type=int,
    default=10,
    help="Minimum number of telomeric repeats for a read to be recorded",
    show_default=True,
)
@click.option(
    "--cov_window_size",
    type=int,
    default=5,
    help="Number of positions either side of a soft-clipped telomere-containing read to record. MUST be >=1!! [This constraint could be removed later.]",
    show_default=True,
)
@click.option("--threads", type=int, default=1)
@click.help_option("--help", "-h")
@click.version_option(__version__, "--version", "-V")
def main(bam_fname, ofname, seq_region, bed, telomere_forward_seq, telo_array_size, cov_window_size, threads):
    """
    Looks for DNA Elimination breakpoints from a bam of reads aligned to a genome.

    OFNAME is the output tsv file. In multiprocessing mode, its containing directory
    will be used to store per-contig intermediate results.
    """
    ofpath = Path(ofname)
    ofpath.parent.mkdir(parents=True, exist_ok=True)
    ofname_base = ofpath.parent / ofpath.stem
    bam_fstream = AlignmentFile(bam_fname)

    contig_names = bam_fstream.references
    seq_regions = it.repeat(seq_region)
    if seq_region is not None:
        seq_regions = [seq_region]
        threads = 1
    if bed is not None:
        intervals = BedTool(bed)
        seq_regions = list()
        for interval in intervals:
            seq_regions.append(f"{interval.chrom}:{interval.start}-{interval.end}")

    telomere_seqs = {Orientation.forward : telomere_forward_seq, Orientation.reverse: rev_comp(telomere_forward_seq)}

    with mp.Pool(processes=threads) as pool:
        pool.starmap(
            find_breakpoint_foci_row_based,
            zip(
                it.repeat(bam_fname),
                it.repeat(ofname_base),
                contig_names,
                it.repeat(telomere_seqs),
                it.repeat(telo_array_size),
                it.repeat(cov_window_size),
                seq_regions,
            ),
        )
    all_files = glob(f"{ofname_base}_*.tsv")
    with ofpath.open("w") as ofstream:
        for i, fname in enumerate(all_files):
            with open(fname) as infile:
                if i > 0:
                    next(infile)  # Skips header
                for line in infile:
                    if line != "\n":
                        ofstream.write(line)
            os.remove(fname)


if __name__ == "__main__":
    main()

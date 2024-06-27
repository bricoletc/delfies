import itertools as it
import multiprocessing as mp
import os
from glob import glob
from pathlib import Path

import rich_click as click
from datasci import Tents
from pybedtools import BedTool
from pysam import AlignmentFile

from delfies import REGION_CLICK_HELP, __version__, all_breakpoint_types
from delfies.breakpoint_foci import (
    cluster_breakpoint_foci,
    extract_breakpoint_sequences,
    find_breakpoint_foci_row_based,
    BreakpointDetectionParams,
)
from delfies.SAM_utils import (
    DEFAULT_MIN_MAPQ,
    DEFAULT_READ_FILTER_FLAG,
    DEFAULT_READ_FILTER_NAMES,
)
from delfies.seq_utils import TELOMERE_SEQS, Orientation, rev_comp

click.rich_click.OPTION_GROUPS = {
    "delfies": [
        {
            "name": "Generic",
            "options": ["--help", "--version", "--threads"],
        },
        {
            "name": "Region selection",
            "options": ["--seq_region", "--bed"],
        },
        {
            "name": "Breakpoint detection",
            "options": [
                "--breakpoint_type",
                "--telomere_forward_seq",
                "--telo_array_size",
                "--cov_window_size",
                "--min_mapq",
                "--read_filter_flag",
            ],
        },
        {
            "name": "Breakpoint extraction",
            "options": [
                "--seq_window_size",
                "--min_supporting_reads",
            ],
        },
    ]
}


@click.command()
@click.argument("genome_fname", type=click.Path(exists=True))
@click.argument("bam_fname", type=click.Path(exists=True))
@click.argument("odirname")
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
    help="The telomere sequence used by your organism. Please make sure this is provided in 'forward' orientation (i.e. 5'->3')",
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
@click.option(
    "--min_mapq",
    type=int,
    default=DEFAULT_MIN_MAPQ,
    help="Reads below this MAPQ will be filtered out",
    show_default=True,
)
@click.option(
    "--read_filter_flag",
    type=int,
    default=DEFAULT_READ_FILTER_FLAG,
    help=f"Reads with any of the component bitwise flags will be filtered out (see SAM specs for details)."
    f"   [default: {DEFAULT_READ_FILTER_FLAG} (reads with any of {DEFAULT_READ_FILTER_NAMES} are filtered out)]",
)
@click.option(
    "--min_supporting_reads",
    type=int,
    default=10,
    help="Minimum number of reads supporting a breakpoint",
    show_default=True,
)
@click.option(
    "--seq_window_size",
    type=int,
    default=350,
    help="Number of nucleotides to extract either side of each identified breakpoint",
    show_default=True,
)
@click.option(
    "--breakpoint_type",
    "-b",
    type=click.Choice(list(map(str, all_breakpoint_types))),
    help="The type of breakpoint to look for. If not provided, will look for all types of breakpoints",
)
@click.option("--threads", type=int, default=1)
@click.help_option("--help", "-h")
@click.version_option(__version__, "--version", "-V")
def main(
    genome_fname,
    bam_fname,
    odirname,
    seq_region,
    bed,
    telomere_forward_seq,
    telo_array_size,
    cov_window_size,
    min_mapq,
    read_filter_flag,
    min_supporting_reads,
    seq_window_size,
    threads,
):
    """
    Looks for DNA Elimination breakpoints from a bam of reads aligned to a genome.

    odirname is the directory to store outputs in. It is also used to store per-contig intermediate results.
    """
    odirname = Path(odirname)
    odirname.mkdir(parents=True, exist_ok=True)
    ofname_base = odirname / "breakpoint_foci"
    bam_fstream = AlignmentFile(bam_fname)

    seq_regions = bam_fstream.references
    if seq_region is not None:
        seq_regions = [seq_region]
        threads = 1
    if bed is not None:
        intervals = BedTool(bed)
        seq_regions = list()
        for interval in intervals:
            seq_regions.append(f"{interval.chrom}:{interval.start}-{interval.end}")

    telomere_seqs = {
        Orientation.forward: telomere_forward_seq,
        Orientation.reverse: rev_comp(telomere_forward_seq),
    }

    detection_params = BreakpointDetectionParams(
        bam_fname=bam_fname,
        ofname_base=ofname_base,
        telomere_seqs=telomere_seqs,
        telo_array_size=telo_array_size,
        cov_window_size=cov_window_size,
        min_mapq=min_mapq,
        read_filter_flag=read_filter_flag,
    )

    with mp.Pool(processes=threads) as pool:
        pool.starmap(
            find_breakpoint_foci_row_based,
            zip(
                it.repeat(detection_params),
                seq_regions,
            ),
        )
    all_files = glob(f"{ofname_base}_*.tsv")
    foci_tsv = f"{ofname_base}.tsv"
    with open(foci_tsv, "w") as ofstream:
        for i, fname in enumerate(all_files):
            with open(fname) as infile:
                if i > 0:
                    next(infile)  # Skips header
                for line in infile:
                    if line != "\n":
                        ofstream.write(line)
            os.remove(fname)

    all_foci = Tents.from_tsv(foci_tsv)
    clustered_foci = cluster_breakpoint_foci(all_foci)
    maximal_foci = map(
        lambda cluster: cluster.find_peak_softclip_focus(), clustered_foci
    )
    maximal_foci = filter(lambda m: m.max_value >= min_supporting_reads, maximal_foci)
    maximal_foci = sorted(maximal_foci, key=lambda e: e.max_value, reverse=True)
    breakpoint_bed = odirname / "breakpoint_maxima.bed"
    with breakpoint_bed.open("w") as ofstream:
        for maximal_focus in maximal_foci:
            strand = "+" if maximal_focus.orientation is Orientation.forward else "-"
            breakpoint_name = f"maximal_focus_window:{maximal_focus.interval[0]}-{maximal_focus.interval[1]}"
            out_line = [
                maximal_focus.focus.contig,
                maximal_focus.focus.start,
                maximal_focus.focus.end,
                breakpoint_name,
                ".",
                strand,
            ]
            ofstream.write("\t".join(map(str, out_line)) + "\n")

    breakpoint_sequences = extract_breakpoint_sequences(
        maximal_foci, genome_fname, seq_window_size
    )
    breakpoint_fasta = odirname / "breakpoint_sequences.fasta"
    with breakpoint_fasta.open("w") as ofstream:
        for breakpoint_sequence in breakpoint_sequences:
            ofstream.write(str(breakpoint_sequence))


if __name__ == "__main__":
    main()

import itertools as it
import multiprocessing as mp
import os
from glob import glob
from pathlib import Path
from typing import List

import rich_click as click
from datasci import Tents
from pybedtools import BedTool
from pyfastx import Fasta
from pysam import AlignmentFile

from delfies import (
    ID_DELIM,
    REGION_CLICK_HELP,
    BreakpointType,
    __version__,
    all_breakpoint_types,
)
from delfies.breakpoint_foci import (
    READ_SUPPORT_PREFIX,
    BreakpointDetectionParams,
    MaximalFocus,
    cluster_breakpoint_foci,
    find_breakpoint_foci_row_based,
)
from delfies.interval_utils import Interval, Intervals
from delfies.SAM_utils import (
    DEFAULT_MIN_MAPQ,
    DEFAULT_READ_FILTER_FLAG,
    DEFAULT_READ_FILTER_NAMES,
)
from delfies.seq_utils import (
    TELOMERE_SEQS,
    FastaRecord,
    Orientation,
    find_all_occurrences_in_genome,
    rev_comp,
)

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
                "--clustering_threshold",
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


MaximalFoci = List[MaximalFocus]


def extract_breakpoint_sequences(
    maximal_foci: MaximalFoci, genome_fasta: str, seq_window_size: int
) -> List[FastaRecord]:
    genome = Fasta(genome_fasta)
    result = list()
    for max_focus in maximal_foci:
        focus = max_focus.focus
        start = int(focus.start)
        if start < 0:
            continue
        breakpoint_sequence = (
            genome.fetch(focus.contig, (start - seq_window_size, start))
            + "N"
            + genome.fetch(focus.contig, (start, start + seq_window_size))
        )
        strand_name = "3prime"
        if max_focus.orientation is Orientation.reverse:
            breakpoint_sequence = rev_comp(breakpoint_sequence)
            strand_name = "5prime"
        breakpoint_name = f"{max_focus.breakpoint_type}_{strand_name}_{focus.contig} breakpoint_pos:{start} {READ_SUPPORT_PREFIX}:{max_focus.max_value} next_best_value_on_same_strand:{max_focus.next_max_value} best_value_on_other_strand:{max_focus.max_value_other_orientation}"
        result.append(FastaRecord(breakpoint_name, breakpoint_sequence))
    return result


def run_breakpoint_detection(
    detection_params: BreakpointDetectionParams, seq_regions: Intervals, threads
) -> MaximalFoci:
    with mp.Pool(processes=threads) as pool:
        pool.starmap(
            find_breakpoint_foci_row_based,
            zip(
                it.repeat(detection_params),
                seq_regions,
            ),
        )
    all_files = glob(f"{detection_params.ofname_base}_*.tsv")
    foci_tsv = f"{detection_params.ofname_base}.tsv"
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
    clustered_foci = cluster_breakpoint_foci(
        all_foci, tolerance=detection_params.clustering_threshold
    )
    maximal_foci = map(
        lambda cluster: cluster.find_peak_softclip_focus(), clustered_foci
    )
    maximal_foci = sorted(maximal_foci, key=lambda e: e.max_value, reverse=True)
    for m_f in maximal_foci:
        m_f.breakpoint_type = detection_params.breakpoint_type
    return maximal_foci


def write_breakpoint_sequences(
    genome_fname: str, maximal_foci: MaximalFoci, odirname: str, seq_window_size: int
) -> None:
    breakpoint_sequences = extract_breakpoint_sequences(
        maximal_foci, genome_fname, seq_window_size
    )
    breakpoint_fasta = odirname / "breakpoint_sequences.fasta"
    with breakpoint_fasta.open("w") as ofstream:
        for breakpoint_sequence in breakpoint_sequences:
            ofstream.write(str(breakpoint_sequence))


def write_breakpoint_bed(maximal_foci: MaximalFoci, odirname: str) -> None:
    breakpoint_bed = odirname / "breakpoint_locations.bed"
    with breakpoint_bed.open("w") as ofstream:
        for maximal_focus in maximal_foci:
            strand = "+" if maximal_focus.orientation is Orientation.forward else "-"
            breakpoint_name = f"Type:{maximal_focus.breakpoint_type};breakpoint_window:{maximal_focus.interval[0]}-{maximal_focus.interval[1]}"
            out_line = [
                maximal_focus.focus.contig,
                maximal_focus.focus.start,
                maximal_focus.focus.end,
                breakpoint_name,
                maximal_focus.max_value,
                strand,
            ]
            ofstream.write("\t".join(map(str, out_line)) + "\n")


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
    "--clustering_threshold",
    type=int,
    default=5,
    help=f"Any breakpoints within this value (in bp) of each other will be merged. "
    f"A larger threshold allows for more imprecise breakpoint locations",
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
    type=click.Choice(list(map(str, all_breakpoint_types)) + ["all"]),
    help="The type of breakpoint to look for. By default, looks for all types of breakpoints",
    default="all",
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
    clustering_threshold,
    min_mapq,
    read_filter_flag,
    min_supporting_reads,
    seq_window_size,
    breakpoint_type,
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

    seq_regions: Intervals = list()
    if bed is not None:
        intervals = BedTool(bed)
        for interval in intervals:
            seq_regions.append(Interval.from_pybedtools_interval(interval))
    elif seq_region is not None:
        seq_regions.append(Interval.from_region_string(seq_region))
        threads = 1
    else:
        # Analyse the entire genome
        for contig in bam_fstream.references:
            seq_regions.append(Interval(contig))

    telomere_seqs = {
        Orientation.forward: telomere_forward_seq,
        Orientation.reverse: rev_comp(telomere_forward_seq),
    }

    clustering_threshold = max(clustering_threshold, 0)
    detection_params = BreakpointDetectionParams(
        bam_fname=bam_fname,
        telomere_seqs=telomere_seqs,
        telo_array_size=telo_array_size,
        clustering_threshold=clustering_threshold,
        min_mapq=min_mapq,
        read_filter_flag=read_filter_flag,
        min_supporting_reads=min_supporting_reads,
    )

    try:
        breakpoint_types_to_analyse = [BreakpointType(breakpoint_type)]
    except ValueError:
        breakpoint_types_to_analyse = all_breakpoint_types

    maximal_foci = []
    for breakpoint_type_to_analyse in breakpoint_types_to_analyse:
        detection_params.breakpoint_type = breakpoint_type_to_analyse
        detection_params.ofname_base = (
            f"{ofname_base}{ID_DELIM}{breakpoint_type_to_analyse}"
        )
        if breakpoint_type_to_analyse is BreakpointType.G2S:
            telomere_query = (
                detection_params.telomere_seqs[Orientation.forward]
                * detection_params.telo_array_size
            )
            genome_fasta = Fasta(genome_fname, build_index=True)
            seq_regions = find_all_occurrences_in_genome(
                telomere_query,
                genome_fasta,
                seq_regions,
                interval_window_size=20,
            )
        maximal_foci += run_breakpoint_detection(detection_params, seq_regions, threads)

    write_breakpoint_bed(maximal_foci, odirname)
    write_breakpoint_sequences(genome_fname, maximal_foci, odirname, seq_window_size)


if __name__ == "__main__":
    main()

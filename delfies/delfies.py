import itertools as it
import multiprocessing as mp
from pathlib import Path
from typing import Tuple

import rich_click as click
from datasci import Tents
from pybedtools import BedTool
from pysam import AlignmentFile

from delfies import (
    ID_DELIM,
    REGION_CLICK_HELP,
    BreakpointType,
    Orientation,
    PutativeBreakpoints,
    __version__,
    all_breakpoint_types,
)
from delfies.breakpoint_foci import (
    BreakpointDetectionParams,
    cluster_breakpoint_foci,
    find_breakpoint_foci,
    setup_breakpoint_tents,
)
from delfies.breakpoint_sequences import write_breakpoint_sequences
from delfies.interval_utils import Interval, Intervals
from delfies.SAM_utils import (
    DEFAULT_MIN_MAPQ,
    DEFAULT_READ_FILTER_FLAG,
    DEFAULT_READ_FILTER_NAMES,
    SoftclippedReads,
    filter_out_read_intervals,
)
from delfies.seq_utils import find_all_occurrences_in_genome, rev_comp
from delfies.telomere_utils import TELOMERE_SEQS, record_telomere_additions

click.rich_click.OPTION_GROUPS = {
    "delfies": [
        {
            "name": "Generic",
            "options": ["--help", "--version", "--threads", "--sample_name"],
        },
        {
            "name": "Region selection",
            "options": ["--seq_region", "--bed"],
        },
        {
            "name": "Breakpoint detection",
            "options": [
                "--breakpoint_type",
                "--telo_forward_seq",
                "--telo_array_size",
                "--telo_max_edit_distance",
                "--min_mapq",
                "--read_filter_flag",
                "--keep_telomeric_breakpoints",
            ],
        },
        {
            "name": "Output breakpoints",
            "options": [
                "--seq_window_size",
                "--min_supporting_reads",
                "--clustering_threshold",
            ],
        },
    ]
}


def filter_out_tents_intervals(input_tents: Tents, target_intervals: Intervals) -> None:
    to_remove = set()
    for i, input_tent in enumerate(input_tents):
        input_interval = Interval(input_tent.contig, input_tent.start, input_tent.end)
        for target_interval in target_intervals:
            if target_interval.overlaps_or_touches(input_interval):
                to_remove.add(i)
                break
    input_tents.drop_elements(to_remove)


def run_breakpoint_detection(
    detection_params: BreakpointDetectionParams,
    seq_regions: Intervals,
    threads: int,
    telomere_array_regions: Intervals,
) -> Tuple[PutativeBreakpoints, SoftclippedReads]:
    with mp.Pool(processes=threads) as pool:
        pooled_results = pool.starmap(
            find_breakpoint_foci,
            zip(
                it.repeat(detection_params),
                seq_regions,
            ),
        )
    all_reads_for_telomere_addition = list()
    all_foci = setup_breakpoint_tents()
    for foci_pool, read_pool in pooled_results:
        all_foci.extend(foci_pool)
        all_reads_for_telomere_addition.extend(read_pool)
    if (
        detection_params.breakpoint_type is BreakpointType.S2G
        and not detection_params.keep_telomeric_breakpoints
    ):
        filter_out_tents_intervals(all_foci, telomere_array_regions)
        all_reads_for_telomere_addition = filter_out_read_intervals(
            all_reads_for_telomere_addition, telomere_array_regions
        )
    foci_tsv = f"{detection_params.ofname_base}.tsv"
    with open(foci_tsv, "w") as ofstream:
        print(all_foci, file=ofstream)
    clustered_foci = cluster_breakpoint_foci(
        all_foci, tolerance=detection_params.clustering_threshold
    )
    putative_breakpoints = map(
        lambda cluster: cluster.find_peak_softclip_focus(), clustered_foci
    )
    putative_breakpoints = sorted(
        putative_breakpoints, key=lambda e: e.max_value, reverse=True
    )
    for m_f in putative_breakpoints:
        m_f.breakpoint_type = detection_params.breakpoint_type
    return putative_breakpoints, all_reads_for_telomere_addition


def write_breakpoint_bed(
    putative_breakpoints: PutativeBreakpoints, sample_name: str, odirname: Path
) -> None:
    breakpoint_bed = odirname / "breakpoint_locations.bed"
    with breakpoint_bed.open("w") as ofstream:
        for putative_breakpoint in putative_breakpoints:
            strand = putative_breakpoint.orientation.value
            breakpoint_name = f"Type:{putative_breakpoint.breakpoint_type};breakpoint_window:{putative_breakpoint.interval[0]}-{putative_breakpoint.interval[1]}"
            out_line = [
                putative_breakpoint.focus.contig,
                putative_breakpoint.focus.start,
                putative_breakpoint.focus.end,
                breakpoint_name,
                putative_breakpoint.max_value,
                strand,
                sample_name,
            ]
            ofstream.write("\t".join(map(str, out_line)) + "\n")


@click.command()
@click.argument("genome_fname", type=click.Path(exists=True))
@click.argument("bam_fname", type=click.Path(exists=True))
@click.argument("odirname")
@click.option(
    "--sample_name",
    type=str,
    help="The name of your sample; will appear in a dedicated column in all output TSV files",
    default=".",
    show_default=True,
)
@click.option("--seq_region", type=str, help=REGION_CLICK_HELP)
@click.option(
    "--bed",
    type=click.Path(exists=True),
    help="Path to bed of regions to analyse. Overrides 'seq_region'",
)
@click.option(
    "--telo_forward_seq",
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
    "--telo_max_edit_distance",
    type=int,
    default=3,
    help="Maximum number of mutations allowed in the searched telomere array",
    show_default=True,
)
@click.option(
    "--clustering_threshold",
    type=int,
    default=5,
    help="Any identified breakpoints within this value (in bp) of each other will be merged. "
    "A larger threshold allows for more imprecise breakpoint locations",
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
    help="The type of breakpoint to look for. By default, looks for all",
    default="all",
)
@click.option(
    "--keep_telomeric_breakpoints",
    is_flag=True,
    help="Forces delfies to keep breakpoints occurring inside telomeric arrays. As these are often false positives, they are discarded by default.",
    show_default=True,
)
@click.option("--threads", type=int, default=1)
@click.help_option("--help", "-h")
@click.version_option(__version__, "--version", "-V")
def main(
    genome_fname,
    bam_fname,
    odirname,
    sample_name,
    seq_region,
    bed,
    telo_forward_seq,
    telo_array_size,
    telo_max_edit_distance,
    clustering_threshold,
    min_mapq,
    read_filter_flag,
    min_supporting_reads,
    seq_window_size,
    breakpoint_type,
    keep_telomeric_breakpoints,
    threads,
):
    """
    Looks for DNA Elimination breakpoints from a bam of reads aligned to a genome.

    odirname is the directory to store outputs in.
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
        Orientation.forward: telo_forward_seq,
        Orientation.reverse: rev_comp(telo_forward_seq),
    }

    clustering_threshold = max(clustering_threshold, 0)
    detection_params = BreakpointDetectionParams(
        bam_fname=bam_fname,
        sample_name=sample_name,
        telomere_seqs=telomere_seqs,
        telo_array_size=telo_array_size,
        max_edit_distance=telo_max_edit_distance,
        clustering_threshold=clustering_threshold,
        min_mapq=min_mapq,
        read_filter_flag=read_filter_flag,
        min_supporting_reads=min_supporting_reads,
        keep_telomeric_breakpoints=keep_telomeric_breakpoints,
    )

    try:
        breakpoint_types_to_analyse = [BreakpointType(breakpoint_type)]
    except ValueError:
        breakpoint_types_to_analyse = all_breakpoint_types

    identified_breakpoints = []
    searched_telo_unit = detection_params.telomere_seqs[Orientation.forward]
    searched_telo_array = searched_telo_unit * detection_params.telo_array_size
    interval_window_size = len(searched_telo_array)
    all_reads_for_telomere_addition = list()
    telomere_array_regions = find_all_occurrences_in_genome(
        searched_telo_array,
        genome_fname,
        seq_regions,
        interval_window_size,
    )
    for breakpoint_type_to_analyse in breakpoint_types_to_analyse:
        detection_params.breakpoint_type = breakpoint_type_to_analyse
        detection_params.ofname_base = (
            f"{ofname_base}{ID_DELIM}{breakpoint_type_to_analyse}"
        )
        if breakpoint_type_to_analyse is BreakpointType.G2S:
            # Restrict regions to analyse to those containing telomere arrays
            seq_regions = telomere_array_regions
        candidate_breakpoints, reads_for_telomere_addition = run_breakpoint_detection(
            detection_params, seq_regions, threads, telomere_array_regions
        )
        all_reads_for_telomere_addition.extend(reads_for_telomere_addition)
        identified_breakpoints.extend(candidate_breakpoints)

    write_breakpoint_bed(identified_breakpoints, sample_name, odirname)
    seq_window_size = max(seq_window_size, 1)
    write_breakpoint_sequences(
        genome_fname, identified_breakpoints, odirname, seq_window_size
    )

    telomere_additions = record_telomere_additions(
        all_reads_for_telomere_addition,
        len(telo_forward_seq),
        genome_fname,
        sample_name,
    )
    telomere_addition_ofname = odirname / "telomere_additions.tsv"
    with telomere_addition_ofname.open("w") as ofstream:
        print(telomere_additions, file=ofstream)


if __name__ == "__main__":
    main()

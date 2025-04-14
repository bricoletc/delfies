from pyfastx import Fasta

from delfies import BreakpointDetectionParams
from delfies.breakpoint_foci import MaximalFoci
from delfies.interval_utils import Interval, Intervals
from delfies.seq_utils import Orientation, find_all_occurrences_in_genome

TELOMERE_SEQS = {
    "Nematoda": {Orientation.forward: "TTAGGC", Orientation.reverse: "GCCTAA"}
}


def find_telomere_arrays(
    genome_fasta: Fasta,
    detection_params: BreakpointDetectionParams,
    seq_regions: Intervals,
) -> Intervals:
    """
    seq_regions: regions in which to look for the arrays
    """
    telomere_query = (
        detection_params.telomere_seqs[Orientation.forward]
        * detection_params.telo_array_size
    )
    result = find_all_occurrences_in_genome(
        telomere_query,
        genome_fasta,
        seq_regions,
        interval_window_size=10,
    )
    return result


def remove_breakpoints_in_telomere_arrays(
    genome_fasta: Fasta,
    detection_params: BreakpointDetectionParams,
    maximal_foci: MaximalFoci,
) -> MaximalFoci:
    result = list()
    telo_array_size = len(
        detection_params.telomere_seqs[Orientation.forward]
        * detection_params.telo_array_size
    )
    for maximal_focus in maximal_foci:
        region_to_search = Interval(
            maximal_focus.focus.contig,
            max(maximal_focus.interval[0] - telo_array_size, 0),
            maximal_focus.interval[1] + telo_array_size,
        )
        telomere_arrays_overlapping_breakpoint = find_telomere_arrays(
            genome_fasta, detection_params, [region_to_search]
        )
        if len(telomere_arrays_overlapping_breakpoint) == 0:
            result.append(maximal_focus)
    return result

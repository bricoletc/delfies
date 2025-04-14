from edlib import align as edlib_align
from pyfastx import Fasta

from delfies import BreakpointDetectionParams, Orientation, PutativeBreakpoints
from delfies.interval_utils import Interval, Intervals
from delfies.SAM_utils import SoftclippedRead
from delfies.seq_utils import find_all_occurrences_in_genome

TELOMERE_SEQS = {
    "Nematoda": {Orientation.forward: "TTAGGC", Orientation.reverse: "GCCTAA"}
}


def has_softclipped_telo_array(
    read: SoftclippedRead,
    orientation: Orientation,
    telomere_seqs,
    min_telo_array_size: int,
    max_edit_distance: int,
) -> bool:
    """
    Note: we allow for the softclipped telo array to start with any cyclic shift
    of the telomeric repeat unit.
    """
    telo_unit = telomere_seqs[orientation]
    telo_array = telo_unit * min_telo_array_size
    subseq_clip_end = len(telo_array) + len(telo_unit)
    if orientation is Orientation.forward:
        end = read.sc_query + subseq_clip_end
        subseq = read.sequence[read.sc_query : end]
    else:
        start = max(read.sc_query + 1 - subseq_clip_end, 0)
        subseq = read.sequence[start : read.sc_query + 1]
    result = edlib_align(
        telo_array, subseq, mode="HW", task="distance", k=max_edit_distance
    )
    found_telo_array = result["editDistance"] != -1
    return found_telo_array


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
    maximal_foci: PutativeBreakpoints,
) -> PutativeBreakpoints:
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

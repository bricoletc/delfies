from edlib import align as edlib_align

from delfies import Orientation, PutativeBreakpoints
from delfies.interval_utils import Interval
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
    searched_telo_array = telo_unit * min_telo_array_size
    subseq_clip_end = len(searched_telo_array) + len(telo_unit)
    if orientation is Orientation.forward:
        end = read.sc_query + subseq_clip_end
        subseq = read.sequence[read.sc_query : end]
    else:
        start = max(read.sc_query + 1 - subseq_clip_end, 0)
        subseq = read.sequence[start : read.sc_query + 1]
    result = edlib_align(
        searched_telo_array, subseq, mode="HW", task="distance", k=max_edit_distance
    )
    found_telo_array = result["editDistance"] != -1
    return found_telo_array


def remove_breakpoints_in_telomere_arrays(
    genome_fname: str,
    searched_telo_array: str,
    interval_window_size: int,
    putative_breakpoints: PutativeBreakpoints,
) -> PutativeBreakpoints:
    result = list()
    telo_array_size = len(searched_telo_array)
    for putative_breakpoint in putative_breakpoints:
        region_to_search = Interval(
            putative_breakpoint.focus.contig,
            max(
                putative_breakpoint.interval[0]
                - telo_array_size
                - interval_window_size,
                0,
            ),
            putative_breakpoint.interval[1] + telo_array_size + interval_window_size,
        )
        telomere_arrays_overlapping_breakpoint = find_all_occurrences_in_genome(
            searched_telo_array, genome_fname, [region_to_search], interval_window_size
        )
        if len(telomere_arrays_overlapping_breakpoint) == 0:
            result.append(putative_breakpoint)
    return result

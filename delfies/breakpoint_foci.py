from collections import defaultdict
from dataclasses import dataclass
from itertools import chain as it_chain
from typing import Dict, List, Tuple

from datasci import Tent, Tents
from pysam import AlignedSegment, AlignmentFile

from delfies import ID_DELIM, BreakpointType
from delfies.interval_utils import Interval, get_contiguous_ranges
from delfies.SAM_utils import (
    find_softclip_at_extremity,
    has_softclipped_telo_array,
    read_flag_matches,
)
from delfies.seq_utils import ORIENTATIONS, Orientation

READ_SUPPORT_PREFIX = "num_supporting_reads"
READ_SUPPORTS = [f"{READ_SUPPORT_PREFIX}{ID_DELIM}{o}" for o in ORIENTATIONS]


@dataclass
class BreakpointDetectionParams:
    bam_fname: str
    telomere_seqs: dict
    telo_array_size: int
    max_edit_distance: int
    clustering_threshold: int
    min_mapq: int
    read_filter_flag: int
    min_supporting_reads: int
    breakpoint_type: str = ""
    ofname_base: str = None


def setup_breakpoint_tents() -> Tents:
    tents_header = [
        "contig",
        "start",
        "end",
        "read_depth",
        "breakpoint_type",
    ] + READ_SUPPORTS
    tents = Tents(header=tents_header, required_header=tents_header[:5], unset_value=0)
    return tents


####################
## Foci detection ##
####################
def focus_has_enough_support(focus_tent, min_support: int) -> bool:
    result = False
    for read_support in READ_SUPPORTS:
        result |= int(focus_tent[read_support]) >= min_support
    return result


def record_softclips(
    aligned_read: AlignedSegment,
    breakpoint_foci: Tents,
    breakpoint_foci_positions,
    detection_params: BreakpointDetectionParams,
    seq_region: Interval,
) -> None:
    for read_support in READ_SUPPORTS:
        orientation = Orientation[read_support.split(ID_DELIM)[1]]
        softclipped_read = find_softclip_at_extremity(aligned_read, orientation)
        if softclipped_read is None:
            continue
        if detection_params.breakpoint_type is BreakpointType.G2S:
            softclipped_telo_array_found = has_softclipped_telo_array(
                softclipped_read,
                orientation,
                detection_params.telomere_seqs,
                min_telo_array_size=3,
                max_edit_distance=detection_params.max_edit_distance,
            )
            softclips_start_inside_target_region = seq_region.spans(
                softclipped_read.sc_ref
            )
            keep_read = (
                softclips_start_inside_target_region
                and not softclipped_telo_array_found
            )
        else:
            softclipped_telo_array_found = has_softclipped_telo_array(
                softclipped_read,
                orientation,
                detection_params.telomere_seqs,
                detection_params.telo_array_size,
                max_edit_distance=detection_params.max_edit_distance,
            )
            keep_read = softclipped_telo_array_found
        if not keep_read:
            continue
        pos_to_commit = softclipped_read.sc_ref
        ref_name = aligned_read.reference_name
        match_tent_key = f"{ref_name}{ID_DELIM}{pos_to_commit}"
        if match_tent_key in breakpoint_foci_positions:
            breakpoint_foci_positions[match_tent_key][read_support] += 1
        else:
            new_tent = breakpoint_foci.new()
            new_tent.update(
                contig=ref_name,
                start=pos_to_commit,
                end=pos_to_commit + 1,
                breakpoint_type=str(detection_params.breakpoint_type),
            )
            new_tent[read_support] += 1
            breakpoint_foci_positions[match_tent_key] = new_tent


def find_breakpoint_foci(
    detection_params: BreakpointDetectionParams,
    seq_region: Interval,
) -> Tents:
    breakpoint_foci = setup_breakpoint_tents()
    breakpoint_foci_positions = {}
    contig_name = seq_region.name
    if seq_region.has_coordinates():
        fetch_args = dict(
            contig=contig_name, start=seq_region.start, stop=seq_region.end
        )
    else:
        fetch_args = dict(contig=contig_name)
    bam_fstream = AlignmentFile(detection_params.bam_fname)
    for aligned_read in bam_fstream.fetch(**fetch_args):
        if aligned_read.mapping_quality < detection_params.min_mapq:
            continue
        if read_flag_matches(aligned_read, detection_params.read_filter_flag):
            continue
        record_softclips(
            aligned_read,
            breakpoint_foci,
            breakpoint_foci_positions,
            detection_params,
            seq_region,
        )
    # Filter for minimum support
    breakpoint_foci_positions = {
        key: val
        for key, val in breakpoint_foci_positions.items()
        if focus_has_enough_support(val, detection_params.min_supporting_reads)
    }
    # Expand to a few positions before and after putative breakpoints: allows users to
    # assess changes in coverage around breakpoints (using the corresponding output tsv)
    positions_to_commit = set()
    for match_key, focus_tent in breakpoint_foci_positions.items():
        committed_position = focus_tent["start"]
        positions_to_commit.update(
            range(committed_position - 2, committed_position + 3)
        )
    record_read_depth_at_breakpoint_foci(
        positions_to_commit,
        breakpoint_foci_positions,
        contig_name,
        breakpoint_foci,
        bam_fstream,
        detection_params,
    )
    return breakpoint_foci


def record_read_depth_at_breakpoint_foci(
    positions_to_commit,
    breakpoint_foci_positions,
    contig_name,
    breakpoint_foci,
    bam_fstream,
    detection_params,
):
    """
    Adds read depth at each position in :positions_to_commit:, 
    using pysam 'pileup'.
    Special cases: 
        - When a breakpoint occurs at first position of a contig, 'pileup' will not 
          record any read depth as the breakpoint position is set to -1. 
          So we manually commit that position to :breakpoint_foci: to ensure it is output.
        - [TODO] When a breakpoint occurs at last position of a contig
    """
    for start, end in get_contiguous_ranges(positions_to_commit):
        if start < 0:
            negative_tent_key = f"{contig_name}{ID_DELIM}-1"
            if negative_tent_key in breakpoint_foci_positions:
                breakpoint_foci.add(breakpoint_foci_positions[negative_tent_key])
        # +1 for `end` because `end` needs to be exclusive in pysam `pileup`
        pileup_args = dict(
            contig=contig_name,
            start=max(start, 0),
            end=max(end + 1, 0),
            flag_filter=detection_params.read_filter_flag,
            min_mapping_quality=detection_params.min_mapq,
            ignore_orphans=False,
            truncate=True,
        )
        for pileup_column in bam_fstream.pileup(**pileup_args):
            read_depth = pileup_column.nsegments
            ref_pos = pileup_column.reference_pos
            tent_key = f"{contig_name}{ID_DELIM}{ref_pos}"
            if tent_key in breakpoint_foci_positions:
                breakpoint_foci_positions[tent_key]["read_depth"] = read_depth
                breakpoint_foci.add(breakpoint_foci_positions[tent_key])
            else:
                new_tent = breakpoint_foci.new()
                new_tent.update(
                    contig=contig_name,
                    start=ref_pos,
                    end=ref_pos + 1,
                    read_depth=read_depth,
                )
                breakpoint_foci.add(new_tent)


#####################
## Foci clustering ##
#####################
@dataclass
class MaximalFocus:
    orientation: Orientation
    max_value: int
    next_max_value: int
    max_value_other_orientation: int
    interval: Tuple[int, int]
    focus: Tent
    breakpoint_type: str = ""

    def update(self, query_focus: Tent):
        query_focus_value = int(
            query_focus[f"{READ_SUPPORT_PREFIX}{ID_DELIM}{self.orientation.name}"]
        )
        if query_focus_value > self.max_value:
            self.next_max_value = self.max_value
            self.max_value = query_focus_value
            self.focus = query_focus
        elif query_focus_value > self.next_max_value:
            self.next_max_value = query_focus_value


MaximalFoci = List[MaximalFocus]


class FociWindow:
    def __init__(self, focus):
        self.foci = [focus]
        self.Min = int(focus.start)
        self.Max = int(focus.end)

    def includes(self, focus: Tent, tolerance: int):
        focus_start_past_end = int(focus.start) > self.Max + tolerance
        focus_end_before_start = int(focus.end) < self.Min - tolerance
        return not focus_start_past_end and not focus_end_before_start

    def add(self, focus):
        self.foci.append(focus)
        start = int(focus.start)
        end = int(focus.end)
        if end > self.Max:
            self.Max = end
        if start < self.Min:
            self.Min = start

    def find_peak_softclip_focus(self) -> MaximalFocus:
        forward_maximum = MaximalFocus(
            Orientation.forward, 0, 0, 0, (self.Min, self.Max), None
        )
        reverse_maximum = MaximalFocus(
            Orientation.reverse, 0, 0, 0, (self.Min, self.Max), None
        )
        for focus in self.foci:
            forward_maximum.update(focus)
            reverse_maximum.update(focus)
        if forward_maximum.max_value > reverse_maximum.max_value:
            max_maximum = forward_maximum
            max_maximum.max_value_other_orientation = reverse_maximum.max_value
        else:
            max_maximum = reverse_maximum
            max_maximum.max_value_other_orientation = forward_maximum.max_value
        return max_maximum

    def __repr__(self):
        return f"[{self.Min},{self.Max}]"


def cluster_breakpoint_foci(foci: Tents, tolerance: int) -> List[FociWindow]:
    """
    Developer note:
        foci without any softclipped-reads are ignored for the purpose of clustering,
        as they are only present in the output tsv to assess coverage changes near breakpoints.
    """
    result: Dict[str, List[FociWindow]] = defaultdict(list)
    for focus in foci:
        if not focus_has_enough_support(focus, 1):
            continue
        contig_windows = result[focus.contig]
        found_window = False
        for elem in contig_windows:
            if elem.includes(focus, tolerance=tolerance):
                elem.add(focus)
                found_window = True
                break
        if not found_window:
            contig_windows.append(FociWindow(focus))
    return list(it_chain(*result.values()))

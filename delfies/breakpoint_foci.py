from collections import defaultdict
from dataclasses import dataclass
from itertools import chain as it_chain
from tempfile import NamedTemporaryFile
from typing import Dict, List, Set, Tuple

from datasci import Tent, Tents
from pysam import AlignedSegment, AlignmentFile

from delfies import ID_DELIM, BreakpointType
from delfies.interval_utils import get_contiguous_ranges, parse_region_string
from delfies.SAM_utils import (
    find_softclip_at_extremity,
    has_softclipped_telo_array,
    read_flag_matches,
)
from delfies.seq_utils import ORIENTATIONS, Orientation, rev_comp

READ_SUPPORT_PREFIX = "num_supporting_reads"
READ_SUPPORTS = [READ_SUPPORT_PREFIX + ID_DELIM + o for o in ORIENTATIONS]


@dataclass
class BreakpointDetectionParams:
    bam_fname: str
    telomere_seqs: dict
    telo_array_size: int
    cov_window_size: int
    min_mapq: int
    read_filter_flag: int
    min_supporting_reads: int
    parse_seq_region: bool
    breakpoint_type: str = ""
    ofname_base: str = None


def setup_tents() -> Dict:
    tents_headers = [
        "contig",
        "start",
        "end",
        "read_depth",
        "breakpoint_type",
    ] + READ_SUPPORTS
    tents = Tents(
        headers=tents_headers, required_headers=tents_headers[:5], unset_value=0
    )
    return tents


def write_tents(ofname_base: str, tents: Tents) -> None:
    ofpath = NamedTemporaryFile(
        prefix=f"{ofname_base}_", suffix=".tsv", dir=".", delete=False
    )
    with open(ofpath.name, "w") as ofstream:
        print(tents, file=ofstream)


PositionTents = Dict[str, Tent]


####################
## Foci detection ##
####################
def record_softclips(
    aligned_read: AlignedSegment,
    tents: Tents,
    position_tents: PositionTents,
    positions_to_commit: Set[int],
    detection_params: BreakpointDetectionParams,
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
            )
            keep_read = not softclipped_telo_array_found
        else:
            softclipped_telo_array_found = has_softclipped_telo_array(
                softclipped_read,
                orientation,
                detection_params.telomere_seqs,
                detection_params.telo_array_size
            )
            keep_read = softclipped_telo_array_found
        if not keep_read:
            continue
        pos_to_commit = softclipped_read.sc_ref
        ref_name = aligned_read.reference_name
        match_tent_key = f"{ref_name}{ID_DELIM}{pos_to_commit}"
        if match_tent_key in position_tents:
            position_tents[match_tent_key][read_support] += 1
        else:
            new_tent = tents.new()
            new_tent.update(
                contig=ref_name,
                start=pos_to_commit,
                end=pos_to_commit + 1,
                breakpoint_type=str(detection_params.breakpoint_type),
            )
            new_tent[read_support] += 1
            position_tents[match_tent_key] = new_tent
        positions_to_commit.update(
            range(
                pos_to_commit - detection_params.cov_window_size,
                pos_to_commit + detection_params.cov_window_size,
            )
        )


def find_breakpoint_foci_row_based(
    detection_params: BreakpointDetectionParams,
    seq_region: str,
) -> None:
    tents = setup_tents()
    position_tents: PositionTents = {}
    positions_to_commit = set()
    if detection_params.parse_seq_region:
        contig_name, start, stop = parse_region_string(seq_region)
        fetch_args = dict(contig=contig_name, start=start, stop=stop)
    else:
        contig_name = seq_region
        fetch_args = dict(contig=contig_name)
    bam_fstream = AlignmentFile(detection_params.bam_fname)
    for aligned_read in bam_fstream.fetch(**fetch_args):
        if aligned_read.mapping_quality < detection_params.min_mapq:
            continue
        if read_flag_matches(aligned_read, detection_params.read_filter_flag):
            continue
        record_softclips(
            aligned_read, tents, position_tents, positions_to_commit, detection_params
        )
    # Record read depth at breakpoint foci
    for start, stop in get_contiguous_ranges(positions_to_commit):
        # Special case: sometimes sotfclipped telomere arrays extend 5' from the first position of a contig/scaffold.
        if start < 0:
            negative_tent_key = f"{contig_name}{ID_DELIM}-1"
            if negative_tent_key in position_tents:
                tents.add(position_tents[negative_tent_key])
                position_tents.pop(negative_tent_key)
        pileup_args = dict(
            contig=contig_name,
            start=max(start, 0),
            stop=max(stop, 0),
            flag_filter=detection_params.read_filter_flag,
            min_mapping_quality=detection_params.min_mapq,
            ignore_orphans=False,
            truncate=True,
        )
        for pileup_column in bam_fstream.pileup(**pileup_args):
            read_depth = pileup_column.nsegments
            ref_pos = pileup_column.reference_pos
            tent_key = f"{contig_name}{ID_DELIM}{ref_pos}"
            if tent_key in position_tents:
                position_tents[tent_key]["read_depth"] = read_depth
                tents.add(position_tents[tent_key])
            else:
                new_tent = tents.new()
                new_tent.update(
                    contig=contig_name,
                    start=ref_pos,
                    end=ref_pos + 1,
                    read_depth=read_depth,
                )
                tents.add(new_tent)
    write_tents(detection_params.ofname_base, tents)


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


def cluster_breakpoint_foci(foci: Tents, tolerance: int = 10) -> List[FociWindow]:
    result: Dict[str, List[FociWindow]] = defaultdict(list)
    for focus in foci:
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

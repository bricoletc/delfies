from collections import defaultdict
from dataclasses import dataclass
from itertools import chain as it_chain
from tempfile import NamedTemporaryFile
from typing import Dict, List, Set

from datasci import Tent, Tents
from pyfastx import Fasta
from pysam import AlignedSegment, AlignmentFile

from delfies import ID_DELIM
from delfies.interval_utils import get_contiguous_ranges, parse_region_string
from delfies.SAM_utils import (
    find_softclip_at_extremity,
    has_softclipped_telo_array,
    read_flag_matches,
)
from delfies.seq_utils import ORIENTATIONS, FastaRecord, Orientation, rev_comp

TELO_FEATURES = ["telo_containing_softclips" + ID_DELIM + o for o in ORIENTATIONS]


def setup_tents() -> Dict:
    tents_headers = ["contig", "start", "end", "read_depth"] + TELO_FEATURES
    tents = Tents(
        headers=tents_headers, required_headers=tents_headers[:4], unset_value=0
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
    cov_window_size: int,
    telomere_seqs: Dict,
    telo_array_size: int,
) -> None:
    for telo_feature in TELO_FEATURES:
        orientation = Orientation[telo_feature.split(ID_DELIM)[1]]
        softclipped_read = find_softclip_at_extremity(aligned_read, orientation)
        if softclipped_read is None:
            continue
        if not has_softclipped_telo_array(
            softclipped_read, orientation, telomere_seqs, telo_array_size
        ):
            continue
        pos_to_commit = softclipped_read.sc_ref
        ref_name = aligned_read.reference_name
        match_tent_key = f"{ref_name}{ID_DELIM}{pos_to_commit}"
        if match_tent_key in position_tents:
            position_tents[match_tent_key][telo_feature] += 1
        else:
            new_tent = tents.new()
            new_tent.update(contig=ref_name, start=pos_to_commit, end=pos_to_commit + 1)
            new_tent[telo_feature] += 1
            position_tents[match_tent_key] = new_tent
        positions_to_commit.update(
            range(pos_to_commit - cov_window_size, pos_to_commit + cov_window_size)
        )


def find_breakpoint_foci_row_based(
    bam_fname: str,
    ofname_base: str,
    contig_name: str,
    telomere_seqs: Dict,
    telo_array_size: int,
    cov_window_size: int,
    min_mapq: int,
    read_filter_flag: int,
    seq_region: str = None,
) -> None:
    tents = setup_tents()
    position_tents: PositionTents = {}
    positions_to_commit = set()
    fetch_args = dict(contig=contig_name)
    if seq_region is not None:
        contig_name, start, stop = parse_region_string(seq_region)
        fetch_args.update(contig=contig_name, start=start, stop=stop)
    bam_fstream = AlignmentFile(bam_fname)
    for aligned_read in bam_fstream.fetch(**fetch_args):
        if aligned_read.mapping_quality < min_mapq:
            continue
        if read_flag_matches(aligned_read, read_filter_flag):
            continue
        record_softclips(
            aligned_read,
            tents,
            position_tents,
            positions_to_commit,
            cov_window_size,
            telomere_seqs,
            telo_array_size,
        )
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
            flag_filter=read_filter_flag,
            min_mapping_quality=min_mapq,
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
    write_tents(ofname_base, tents)

from typing import Dict
import re
from tempfile import NamedTemporaryFile

from datasci import Tents
from pysam import AlignmentFile

from delfies.seq_utils import (
    ID_DELIM,
    Orientation,
    ORIENTATIONS,
    REGION_DELIM1,
    REGION_DELIM2,
)
from delfies.SAM_utils import FLAGS, SoftclippedRead, find_softclip_at_extremity

FEATURES = ["telo_containing_softclips" + ID_DELIM + o for o in ORIENTATIONS]
READ_FILTER_OUT = (
    FLAGS["UNMAP"] | FLAGS["SECONDARY"] | FLAGS["DUP"] | FLAGS["SUPPLEMENTARY"]
)
MIN_MAPQ = 20


def has_softclipped_telo_array(
    read: SoftclippedRead, orientation: Orientation, telomere_seqs, telo_array_size: int
) -> bool:
    telo_array = telomere_seqs[orientation] * telo_array_size
    if orientation is Orientation.forward:
        subseq = read.sequence[read.sc_query :]
    else:
        subseq = read.sequence[: read.sc_query + 1]
    return re.search(telo_array, subseq) is not None


def find_breakpoint_foci(
    bam_fname: str,
    ofname_base: str,
    contig_name: str,
    telomere_seqs: Dict,
    telo_array_size: int,
    cov_window_size: int,
    seq_region: str = None,
):
    """
    A note on how the search works:
        - We extract positions in the reference at which telomere-containing soft-clipped reads occur, plus
        a window around that, to capture coverage changes
        - We process 'pileups': one set of aligned reads per genome position, moving 5'->3' on the genome.
          Note that pileups don't include soft/hard-clipped bases.
        - Reads are only processed once
        - Statistics are committed as early as possible - saves memory
    An note on **when** we commit statistics:
        - If softclips occur in reverse orientation, they will occur before the position at which reads containing them are processed -
          because we only process reads once, using the 'is_head' property, which excludes soft/hard clips
        - Therefore we must not commit the stats on a position before such reads are processed. Thus we must at least not commit the
          previous position from the one we are currently processing. In other words we can only commit a position after the position
          after it is committed.
    """
    tents_headers = ["contig", "start", "end", "read_depth"] + FEATURES
    tents = Tents(
        headers=tents_headers, required_headers=tents_headers[:4], unset_value=0
    )
    position_tents: Dict[str, Tent] = {}
    pileup_args = dict(
        flag_filter=READ_FILTER_OUT, min_mapping_quality=MIN_MAPQ, ignore_orphans=False
    )
    if seq_region is not None:
        contig, regs = seq_region.split(REGION_DELIM1)
        start, stop = map(lambda e: int(e.replace(",", "")), regs.split(REGION_DELIM2))
        pileup_args.update(contig=contig, start=start, stop=stop)
    else:
        pileup_args.update(contig=contig_name)
    bam_fstream = AlignmentFile(bam_fname)
    all_pileups = bam_fstream.pileup(**pileup_args)

    positions_to_commit = set()
    for pileup_column in all_pileups:
        ref_name = pileup_column.reference_name
        ref_pos = pileup_column.reference_pos
        if ref_pos % 100000 == 0:
            print(f"Processed pos: {ref_name}:{ref_pos}")
        read_depth = pileup_column.nsegments
        cur_tent_key = f"{ref_name}{ID_DELIM}{ref_pos}"
        if cur_tent_key in position_tents:
            position_tents[cur_tent_key]["read_depth"] = read_depth
        else:
            new_tent = tents.new()
            new_tent.update(
                contig=ref_name,
                start=ref_pos,
                end=ref_pos + 1,
                read_depth=read_depth,
            )
            position_tents[cur_tent_key] = new_tent
        for pileup_read in pileup_column.pileups:
            if not pileup_read.is_head:  # only process a read once
                continue
            aligned_read = pileup_read.alignment
            for feature in FEATURES:
                orientation = Orientation[feature.split(ID_DELIM)[1]]
                softclipped_read = find_softclip_at_extremity(aligned_read, orientation)
                if softclipped_read is None:
                    continue
                if not has_softclipped_telo_array(
                    softclipped_read, orientation, telomere_seqs, telo_array_size
                ):
                    continue
                pos_to_commit = softclipped_read.sc_ref
                match_tent_key = f"{ref_name}{ID_DELIM}{pos_to_commit}"
                if match_tent_key in position_tents:
                    position_tents[match_tent_key][feature] += 1
                else:
                    new_tent = tents.new()
                    new_tent.update(
                        contig=ref_name, start=pos_to_commit, end=pos_to_commit + 1
                    )
                    new_tent[feature] += 1
                    position_tents[match_tent_key] = new_tent
                positions_to_commit.update(
                    range(
                        pos_to_commit - cov_window_size, pos_to_commit + cov_window_size
                    )
                )
        rolling_position = ref_pos - cov_window_size
        tent_key_to_remove = f"{ref_name}{ID_DELIM}{rolling_position}"
        if tent_key_to_remove in position_tents:
            putatively_committed_tent = position_tents.pop(tent_key_to_remove)
            if rolling_position in positions_to_commit:
                positions_to_commit.remove(rolling_position)
                tents.add(putatively_committed_tent)
    ofpath = NamedTemporaryFile(
        prefix=f"{ofname_base}_", suffix=".tsv", dir=".", delete=False
    )
    with open(ofpath.name, "w") as ofstream:
        print(tents, file=ofstream)

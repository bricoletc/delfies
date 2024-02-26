from typing import Dict, Set
import re
from tempfile import NamedTemporaryFile

from datasci import Tent, Tents
from pysam import AlignmentFile, AlignedSegment

from delfies.seq_utils import ID_DELIM, Orientation, ORIENTATIONS, parse_region_string
from delfies.SAM_utils import FLAGS, SoftclippedRead, find_softclip_at_extremity
from delfies.num_utils import get_contiguous_ranges

TELO_FEATURES = ["telo_containing_softclips" + ID_DELIM + o for o in ORIENTATIONS]
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
) -> None:
    """
    A note on how the search works:
        - We extract positions in the reference at which telomere-containing soft-clipped reads occur, plus
        a window around that, to capture coverage changes
        - We process 'pileups': one set of aligned reads per genome position, moving 5'->3' on the genome.
          Note that pileups don't include soft/hard-clipped bases.
        - Reads are only processed once
        - Statistics are committed as early as possible - saves memory
    An note on **when** we commit statistics:
        - We could commit a position as soon as it is processed (5'->3'), but this would only work
          for softclipped reads in the forward orientation, as a softclip-containing position
          will have been recorded before we process it.
        - For softclips occurring in reverse orientation, this strategy fails, as they occur
          before the position at which reads containing them are processed.
          To deal with this, we commit a position only after it has been processed. This way,
          when sofclips occur upstream (5') of a processed read, they can only get committed
          after the condition 'is_head' for the read containing them holds true.
          (Note: 'is_head' holds true when we are at the left-most mapped base of a read,
          regardless of its mapped orientation, or the orientation of its softclips).
        - For that strategy to work, 'cov_window_size' must be > 1.
    """
    tents = setup_tents()
    position_tents: PositionTents = {}
    positions_to_commit = set()
    pileup_args = dict(
        flag_filter=READ_FILTER_OUT,
        min_mapping_quality=MIN_MAPQ,
        ignore_orphans=False,
        truncate=True,
    )
    if seq_region is not None:
        contig, start, stop = parse_region_string(seq_region)
        pileup_args.update(contig=contig, start=start, stop=stop)
    else:
        pileup_args.update(contig=contig_name)
    bam_fstream = AlignmentFile(bam_fname)
    all_pileups = bam_fstream.pileup(**pileup_args)

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
            record_softclips(aligned_read, tents, position_tents, positions_to_commit, cov_window_size, telomere_seqs, telo_array_size)
        putatively_committed_position = ref_pos - cov_window_size
        putatively_committed_tent_key = (
            f"{ref_name}{ID_DELIM}{putatively_committed_position}"
        )
        if putatively_committed_tent_key not in position_tents:
            continue
        putatively_committed_tent = position_tents.pop(putatively_committed_tent_key)
        if putatively_committed_position in positions_to_commit:
            positions_to_commit.remove(putatively_committed_position)
            tents.add(putatively_committed_tent)
    write_tents(ofname_base, tents)

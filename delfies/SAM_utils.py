from dataclasses import dataclass
from typing import Optional

from pysam import CSOFT_CLIP, AlignedSegment

from delfies.seq_utils import Orientation


@dataclass
class SoftclippedRead:
    """
    `sc`: position of softclip start, located at read extremity
    """

    sequence: str
    name: str
    sc_ref: int
    sc_query: int


_ordered_flags = [
    "PAIRED",
    "PROPER_PAIR",
    "UNMAP",
    "MUNMAP",
    "REVERSE",
    "MREVERSE",
    "READ1",
    "READ2",
    "SECONDARY",
    "QCFAIL",
    "DUP",
    "SUPPLEMENTARY",
]
FLAGS = {key: 2**x for x, key in enumerate(_ordered_flags)}


def find_softclip_at_extremity(
    read: AlignedSegment, orientation: Orientation
) -> Optional[SoftclippedRead]:
    """
    pysam (used version 0.20.0): attributes `reference_start` and `reference_end` refer to
    left and rightmost position (on reference; i.e. irrespective of 'reverse' flag) of
    consumed alignment positions (i.e. excluding hard and softclips, etc; see SAM spec).

    `reference_start`: first aligned reside (0-based)
    `reference_end`: 1 past the last aligned residue (0-based)
    """
    result = SoftclippedRead(read.query_sequence, read.query_name, None, None)
    if orientation is Orientation.forward:
        if read.cigartuples[-1][0] == CSOFT_CLIP:
            result.sc_ref = read.reference_end
            result.sc_query = read.query_alignment_end
    else:
        if read.cigartuples[0][0] == CSOFT_CLIP:
            result.sc_ref = read.reference_start - 1
            result.sc_query = read.query_alignment_start - 1
    if result.sc_ref is None:
        return None
    else:
        return result
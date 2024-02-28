import pytest
from pysam import AlignedSegment, CSOFT_CLIP, CMATCH

from delfies.SAM_utils import find_softclip_at_extremity, SoftclippedRead
from delfies.seq_utils import TELOMERE_SEQS, Orientation

DEFAULT_ALIGNED_SEQ = "ATGCAAAAAAAAATTTGGA"
# This could be anything, but I'm using telomeric sequence, for sake of realism :-)
DEFAULT_ADDED_SOFCLIPPED_SEQ = TELOMERE_SEQS["Nematoda"][Orientation.forward] * 3


@pytest.fixture
def pysam_basic_read():
    aligned_seq = DEFAULT_ALIGNED_SEQ
    len_aligned_seq = len(aligned_seq)
    read = AlignedSegment()
    read.query_sequence = aligned_seq
    read.query_name = "test_query"
    read.reference_start = 200
    read.cigartuples = [(CMATCH, len_aligned_seq)]
    return read


@pytest.fixture
def pysam_read_with_3prime_softclips(pysam_basic_read):
    telo_seq_added = DEFAULT_ADDED_SOFCLIPPED_SEQ
    read = pysam_basic_read
    read.query_sequence += telo_seq_added
    read.cigartuples += [(CSOFT_CLIP, len(telo_seq_added))]
    return read


@pytest.fixture
def pysam_read_with_5prime_softclips(pysam_basic_read):
    telo_seq_added = DEFAULT_ADDED_SOFCLIPPED_SEQ
    read = pysam_basic_read
    read.query_sequence = telo_seq_added + read.query_sequence
    read.cigartuples = [(CSOFT_CLIP, len(telo_seq_added))] + read.cigartuples
    return read


class TestSoftclipDetection:
    def test_read_no_softclips(self, pysam_basic_read):
        result_forward = find_softclip_at_extremity(
            pysam_basic_read, Orientation.forward
        )
        result_reverse = find_softclip_at_extremity(
            pysam_basic_read, Orientation.reverse
        )
        assert result_forward is None
        assert result_reverse is None

    def test_read_with_3prime_softclips(self, pysam_read_with_3prime_softclips):
        result_forward = find_softclip_at_extremity(
            pysam_read_with_3prime_softclips, Orientation.forward
        )
        assert result_forward == SoftclippedRead(
            sequence=pysam_read_with_3prime_softclips.query_sequence,
            name=pysam_read_with_3prime_softclips.query_name,
            sc_ref=pysam_read_with_3prime_softclips.reference_start
            + len(DEFAULT_ALIGNED_SEQ),
            sc_query=len(DEFAULT_ALIGNED_SEQ),
        )
        assert (
            find_softclip_at_extremity(
                pysam_read_with_3prime_softclips, Orientation.reverse
            )
            is None
        )

    def test_read_with_5prime_softclips(self, pysam_read_with_5prime_softclips):
        result_reverse = find_softclip_at_extremity(
            pysam_read_with_5prime_softclips, Orientation.reverse
        )
        assert result_reverse == SoftclippedRead(
            sequence=pysam_read_with_5prime_softclips.query_sequence,
            name=pysam_read_with_5prime_softclips.query_name,
            sc_ref=pysam_read_with_5prime_softclips.reference_start - 1,
            sc_query=len(DEFAULT_ADDED_SOFCLIPPED_SEQ) - 1,
        )
        assert (
            find_softclip_at_extremity(
                pysam_read_with_5prime_softclips, Orientation.forward
            )
            is None
        )

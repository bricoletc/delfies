import pytest
from pysam import CMATCH, CSOFT_CLIP, AlignedSegment

from delfies import Orientation
from delfies.SAM_utils import (
    FLAGS,
    SoftclippedRead,
    find_softclip_at_extremity,
    read_flag_matches,
)

DEFAULT_ALIGNED_SEQ = "ATGCAAAAAAAAATTTGGA"
DEFAULT_SOFTCLIPPED_SEQ = "TATAGGTAACATCGCGGCATTCTACGG"


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
    seq_added = DEFAULT_SOFTCLIPPED_SEQ
    read = pysam_basic_read
    read.query_sequence += seq_added
    read.cigartuples += [(CSOFT_CLIP, len(seq_added))]
    return read


@pytest.fixture
def pysam_read_with_5prime_softclips(pysam_basic_read):
    seq_added = DEFAULT_SOFTCLIPPED_SEQ
    read = pysam_basic_read
    read.query_sequence = seq_added + read.query_sequence
    read.cigartuples = [(CSOFT_CLIP, len(seq_added))] + read.cigartuples
    return read


class TestReadFiltering:
    read_flags = FLAGS["UNMAP"] | FLAGS["DUP"]

    def test_read_flag_matches(self, pysam_basic_read):
        pysam_basic_read.flag = self.read_flags
        assert read_flag_matches(pysam_basic_read, FLAGS["PAIRED"] | FLAGS["UNMAP"])

    def test_read_flag_no_matches(self, pysam_basic_read):
        pysam_basic_read.flag = self.read_flags
        assert not read_flag_matches(
            pysam_basic_read, FLAGS["PAIRED"] | FLAGS["SUPPLEMENTARY"]
        )

    def test_read_flag_zero_no_matches(self, pysam_basic_read):
        # Flag of zero == read is mapped in forward orientation
        # On a read with such a flag, no filtering is applicable
        pysam_basic_read.flag = 0
        assert not read_flag_matches(pysam_basic_read, FLAGS["PAIRED"] | FLAGS["UNMAP"])


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
            sc_length=len(DEFAULT_SOFTCLIPPED_SEQ),
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
            sc_query=len(DEFAULT_SOFTCLIPPED_SEQ) - 1,
            sc_length=len(DEFAULT_SOFTCLIPPED_SEQ),
        )
        assert (
            find_softclip_at_extremity(
                pysam_read_with_5prime_softclips, Orientation.forward
            )
            is None
        )

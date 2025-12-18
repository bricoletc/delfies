import pytest
from pysam import CMATCH, CSOFT_CLIP, AlignedSegment, AlignmentHeader

from delfies import Orientation
from delfies.interval_utils import Interval
from delfies.SAM_utils import (
    FLAGS,
    SoftclippedRead,
    filter_out_read_intervals,
    find_softclip_at_extremity,
    read_flag_matches,
)

DEFAULT_ALIGNED_SEQ = "ATGCAAAAAAAAATTTGGA"
DEFAULT_SOFTCLIPPED_SEQ = "TATAGGTAACATCGCGGCATTCTACGG"
DEFAULT_REF_NAME = "test_reference"


@pytest.fixture
def pysam_basic_read():
    header = AlignmentHeader.from_references(
        reference_names=[DEFAULT_REF_NAME], reference_lengths=[100000]
    )
    aligned_seq = DEFAULT_ALIGNED_SEQ
    len_aligned_seq = len(aligned_seq)
    read = AlignedSegment(header=header)
    read.query_sequence = aligned_seq
    read.query_name = "test_query"
    read.reference_name = DEFAULT_REF_NAME
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


class TestReadBitWiseFiltering:
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
            reference_name=pysam_read_with_3prime_softclips.reference_name,
            sc_ref=pysam_read_with_3prime_softclips.reference_start
            + len(DEFAULT_ALIGNED_SEQ),
            sc_query=len(DEFAULT_ALIGNED_SEQ),
            sc_length=len(DEFAULT_SOFTCLIPPED_SEQ),
            orientation=Orientation.forward,
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
            reference_name=pysam_read_with_5prime_softclips.reference_name,
            sc_ref=pysam_read_with_5prime_softclips.reference_start - 1,
            sc_query=len(DEFAULT_SOFTCLIPPED_SEQ) - 1,
            sc_length=len(DEFAULT_SOFTCLIPPED_SEQ),
            orientation=Orientation.reverse,
        )
        assert (
            find_softclip_at_extremity(
                pysam_read_with_5prime_softclips, Orientation.forward
            )
            is None
        )


@pytest.fixture
def basic_softclipped_read():
    read = SoftclippedRead(
        DEFAULT_ALIGNED_SEQ,
        DEFAULT_REF_NAME,
        sc_ref=400,
        sc_query=200,
        sc_length=50,
        orientation=Orientation.forward,
    )
    return read


class TestReadIntervalFiltering:
    default_range = [Interval(DEFAULT_REF_NAME, 100, 150)]

    def test_read_outside_range_not_filtered(self, basic_softclipped_read):
        result = filter_out_read_intervals([basic_softclipped_read], self.default_range)
        expected = [basic_softclipped_read]
        assert result == expected

    def test_inside_range_filtered(self, basic_softclipped_read):
        basic_softclipped_read.sc_ref = 101
        result = filter_out_read_intervals([basic_softclipped_read], self.default_range)
        expected = []
        assert result == expected

    def test_at_range_exclusion_point_not_filtered(self, basic_softclipped_read):
        """
        We want the ref position just before the softclipped position to be considered
        """
        basic_softclipped_read.sc_ref = 100
        result = filter_out_read_intervals([basic_softclipped_read], self.default_range)
        expected = [basic_softclipped_read]
        assert result == expected

    def test_at_range_exclusion_point_reverse_not_filtered(
        self, basic_softclipped_read
    ):
        """
        We want the ref position just after the softclipped position to be considered
        """
        basic_softclipped_read.orientation = Orientation.reverse
        basic_softclipped_read.sc_ref = 149
        result = filter_out_read_intervals([basic_softclipped_read], self.default_range)
        expected = [basic_softclipped_read]
        assert result == expected

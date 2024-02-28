import pytest
from pysam import CMATCH, CSOFT_CLIP, AlignedSegment

from delfies.SAM_utils import (
    FLAGS,
    SoftclippedRead,
    read_flag_matches,
    find_softclip_at_extremity,
    has_softclipped_telo_array,
)
from delfies.seq_utils import TELOMERE_SEQS, Orientation

DEFAULT_ALIGNED_SEQ = "ATGCAAAAAAAAATTTGGA"
DEFAULT_TELO_DICT = TELOMERE_SEQS["Nematoda"]
DEFAULT_MIN_TELO_ARRAY_SIZE = 3
DEFAULT_ADDED_SOFCLIPPED_SEQ = (
    DEFAULT_TELO_DICT[Orientation.forward] * DEFAULT_MIN_TELO_ARRAY_SIZE
)


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


@pytest.fixture
def softclipped_read():
    return SoftclippedRead(
        sequence=DEFAULT_ALIGNED_SEQ,
        name="test_sofclipped_read",
        sc_ref=200,
        sc_query=None,
    )


class TestTeloArrayDetection:
    def test_softclipped_read_with_nontelomeric_3prime_softclips(
        self, softclipped_read
    ):
        # 6 is the length of Nematoda telomeric repeat unit
        added_softclips = "A" * 6 * DEFAULT_MIN_TELO_ARRAY_SIZE
        softclipped_read.sequence += added_softclips
        assert not has_softclipped_telo_array(
            softclipped_read,
            Orientation.forward,
            DEFAULT_TELO_DICT,
            DEFAULT_MIN_TELO_ARRAY_SIZE,
        )

    def test_softclipped_read_with_telomeric_3prime_softclips(self, softclipped_read):
        added_softclips = (
            DEFAULT_TELO_DICT[Orientation.forward] * DEFAULT_MIN_TELO_ARRAY_SIZE
        )
        softclipped_read.sequence += added_softclips
        softclipped_read.sc_query = len(DEFAULT_ALIGNED_SEQ)
        assert has_softclipped_telo_array(
            softclipped_read,
            Orientation.forward,
            DEFAULT_TELO_DICT,
            DEFAULT_MIN_TELO_ARRAY_SIZE,
        )

    def test_softclipped_read_with_telomeric_5prime_softclips(self, softclipped_read):
        added_softclips = (
            DEFAULT_TELO_DICT[Orientation.reverse] * DEFAULT_MIN_TELO_ARRAY_SIZE
        )
        softclipped_read.sequence = added_softclips + softclipped_read.sequence
        softclipped_read.sc_query = len(added_softclips) - 1
        assert has_softclipped_telo_array(
            softclipped_read,
            Orientation.reverse,
            DEFAULT_TELO_DICT,
            DEFAULT_MIN_TELO_ARRAY_SIZE,
        )

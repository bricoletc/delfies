import pytest
from pysam import CMATCH, CSOFT_CLIP, AlignedSegment

from delfies import Orientation
from delfies.SAM_utils import (
    FLAGS,
    SoftclippedRead,
    find_softclip_at_extremity,
    read_flag_matches,
)
from delfies.seq_utils import cyclic_shifts, randomly_substitute, rev_comp
from delfies.telomere_utils import TELOMERE_SEQS, has_softclipped_telo_array

DEFAULT_ALIGNED_SEQ = "ATGCAAAAAAAAATTTGGA"
DEFAULT_TELO_DICT = TELOMERE_SEQS["Nematoda"]
DEFAULT_NON_TELO_UNIT_FORWARD = "AAAAAA"
DEFAULT_MIN_TELO_ARRAY_SIZE = 3
DEFAULT_FORWARD_TELO = DEFAULT_TELO_DICT[Orientation.forward]
DEFAULT_FORWARD_TELO_ARRAY = DEFAULT_FORWARD_TELO * DEFAULT_MIN_TELO_ARRAY_SIZE
DEFAULT_REVERSE_TELO_ARRAY = (
    DEFAULT_TELO_DICT[Orientation.reverse] * DEFAULT_MIN_TELO_ARRAY_SIZE
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
    telo_seq_added = DEFAULT_FORWARD_TELO_ARRAY
    read = pysam_basic_read
    read.query_sequence += telo_seq_added
    read.cigartuples += [(CSOFT_CLIP, len(telo_seq_added))]
    return read


@pytest.fixture
def pysam_read_with_5prime_softclips(pysam_basic_read):
    telo_seq_added = DEFAULT_FORWARD_TELO_ARRAY
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
            sc_length=len(DEFAULT_FORWARD_TELO_ARRAY),
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
            sc_query=len(DEFAULT_FORWARD_TELO_ARRAY) - 1,
            sc_length=len(DEFAULT_FORWARD_TELO_ARRAY),
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
        name="test_softclipped_read",
        sc_ref=200,
        sc_query=len(DEFAULT_ALIGNED_SEQ),
        sc_length=0,
    )


class TestTeloArrayDetection:
    def test_softclipped_read_with_nontelomeric_3prime_softclips_is_not_recognised(
        self, softclipped_read
    ):
        added_softclips = DEFAULT_NON_TELO_UNIT_FORWARD * DEFAULT_MIN_TELO_ARRAY_SIZE
        softclipped_read.sequence += added_softclips
        assert not has_softclipped_telo_array(
            softclipped_read,
            Orientation.forward,
            DEFAULT_TELO_DICT,
            DEFAULT_MIN_TELO_ARRAY_SIZE,
            0,
        )

    def test_softclipped_read_with_telomeric_3prime_softclips_is_recognised(
        self, softclipped_read
    ):
        softclipped_read.sequence += DEFAULT_FORWARD_TELO_ARRAY
        assert has_softclipped_telo_array(
            softclipped_read,
            Orientation.forward,
            DEFAULT_TELO_DICT,
            DEFAULT_MIN_TELO_ARRAY_SIZE,
            0,
        )

    def test_softclipped_read_with_telomeric_5prime_softclips_is_recognised(
        self, softclipped_read
    ):
        softclipped_read.sequence = DEFAULT_REVERSE_TELO_ARRAY
        softclipped_read.sc_query = len(DEFAULT_REVERSE_TELO_ARRAY) - 1
        assert has_softclipped_telo_array(
            softclipped_read,
            Orientation.reverse,
            DEFAULT_TELO_DICT,
            DEFAULT_MIN_TELO_ARRAY_SIZE,
            0,
        )

    def test_softclipped_read_with_shifted_telomeric_softclips_is_recognised(
        self, softclipped_read
    ):
        """
        Telomeres could be added starting with any of the bases in the telomeric
        repeat unit. This tests checks these are tolerated
        """
        read_seq = softclipped_read.sequence
        for orientation in Orientation:
            telo_unit = DEFAULT_TELO_DICT[orientation]
            for shifted_telo_unit in cyclic_shifts(telo_unit):
                added_softclips = (
                    shifted_telo_unit + shifted_telo_unit * DEFAULT_MIN_TELO_ARRAY_SIZE
                )
                if orientation is Orientation.forward:
                    softclipped_read.sequence = read_seq + added_softclips
                else:
                    softclipped_read.sequence = added_softclips + read_seq
                    softclipped_read.sc_query = len(added_softclips) - 1
                assert has_softclipped_telo_array(
                    softclipped_read,
                    orientation,
                    DEFAULT_TELO_DICT,
                    DEFAULT_MIN_TELO_ARRAY_SIZE,
                    0,
                )

    def test_softclipped_read_with_downstream_telomeric_softclips_is_not_recognised(
        self, softclipped_read
    ):
        """
        We only want to match telomere repeats found at the beginning of the
        softclipped portion of the `softclipped_read`.
        """
        read_seq = softclipped_read.sequence
        for orientation in Orientation:
            if orientation is Orientation.forward:
                non_telo_unit = DEFAULT_NON_TELO_UNIT_FORWARD
                added_softclips = (
                    non_telo_unit * DEFAULT_MIN_TELO_ARRAY_SIZE
                    + DEFAULT_FORWARD_TELO_ARRAY
                )
                softclipped_read.sequence = read_seq + added_softclips
            else:
                non_telo_unit = rev_comp(DEFAULT_NON_TELO_UNIT_FORWARD)
                added_softclips = (
                    DEFAULT_REVERSE_TELO_ARRAY
                    + non_telo_unit * DEFAULT_MIN_TELO_ARRAY_SIZE
                )
                softclipped_read.sequence = added_softclips + read_seq
                softclipped_read.sc_query = len(added_softclips) - 1
            assert non_telo_unit != DEFAULT_TELO_DICT[orientation]
            assert not has_softclipped_telo_array(
                softclipped_read,
                orientation,
                DEFAULT_TELO_DICT,
                DEFAULT_MIN_TELO_ARRAY_SIZE,
                0,
            )

    def test_softclipped_read_with_mutated_telomeres_is_recognised(
        self, softclipped_read
    ):
        mutated_telo_array = str()
        for i in range(DEFAULT_MIN_TELO_ARRAY_SIZE):
            mutated_telo_array += randomly_substitute(
                DEFAULT_FORWARD_TELO, num_mutations=1
            )
        for orientation in Orientation:
            if orientation is Orientation.forward:
                softclipped_read.sequence += mutated_telo_array
            else:
                softclipped_read.sequence = (
                    rev_comp(mutated_telo_array) + softclipped_read.sequence
                )
                softclipped_read.sc_query = len(mutated_telo_array) - 1
            assert has_softclipped_telo_array(
                softclipped_read,
                orientation,
                DEFAULT_TELO_DICT,
                DEFAULT_MIN_TELO_ARRAY_SIZE,
                max_edit_distance=DEFAULT_MIN_TELO_ARRAY_SIZE,
            )

    def test_max_edit_distance_threshold_is_applied(self, softclipped_read):
        mutated_telo_array = str()
        for i in range(DEFAULT_MIN_TELO_ARRAY_SIZE):
            mutated_telo_array += randomly_substitute(
                DEFAULT_FORWARD_TELO, num_mutations=1
            )
            softclipped_read.sequence += mutated_telo_array
            assert not has_softclipped_telo_array(
                softclipped_read,
                Orientation.forward,
                DEFAULT_TELO_DICT,
                DEFAULT_MIN_TELO_ARRAY_SIZE,
                max_edit_distance=1,
            )

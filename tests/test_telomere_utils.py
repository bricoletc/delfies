import pytest

from delfies import Orientation
from delfies.SAM_utils import SoftclippedRead
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
def softclipped_read():
    return SoftclippedRead(
        sequence=DEFAULT_ALIGNED_SEQ,
        name="test_softclipped_read",
        sc_ref=200,
        sc_query=len(DEFAULT_ALIGNED_SEQ),
        sc_length=0,
    )


class TestFindSoftClippedTeloArraysInReads:
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

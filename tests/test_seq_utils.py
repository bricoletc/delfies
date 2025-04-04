import pytest

from delfies.interval_utils import Interval
from delfies.seq_utils import (
    cyclic_shifts,
    find_all_occurrences_in_genome,
    randomly_substitute,
    rev_comp,
)
from tests import ClassWithTempFasta


def test_rev_comp():
    seq1 = "AATTAACCGG"
    expected_seq1 = "CCGGTTAATT"
    assert rev_comp(seq1) == expected_seq1


class TestSequenceMutations:
    def test_not_a_nucleotide_fails(self):
        seq_to_mutate = "PPDD"
        with pytest.raises(ValueError):
            randomly_substitute(seq_to_mutate)

    def test_too_many_mutations_fails(self):
        seq_to_mutate = "AATT"
        with pytest.raises(ValueError):
            randomly_substitute(seq_to_mutate, num_mutations=len(seq_to_mutate) + 1)

    def test_various_numbers_of_mutations(self):
        seq_to_mutate = "AATTCCGG"
        num_nucleotides = len(seq_to_mutate)
        for num_mutations in range(num_nucleotides + 1):
            mutated_seq = randomly_substitute(seq_to_mutate, num_mutations)
            num_found_mutations = 0
            for i in range(num_nucleotides):
                if seq_to_mutate[i] != mutated_seq[i]:
                    num_found_mutations += 1
            assert num_found_mutations == num_mutations


def test_cyclic_shifts():
    str_to_shift = "TTAGGC"
    expected_shifts = [
        "TTAGGC",
        "TAGGCT",
        "AGGCTT",
        "GGCTTA",
        "GCTTAG",
        "CTTAGG",
    ]
    result = cyclic_shifts(str_to_shift)
    assert result == expected_shifts


class TestFindOccurrencesInGenomePerfectTeloArrays(ClassWithTempFasta):
    chrom_name = "chr1"
    telo_unit = "TTAGGC"
    default_query_array = telo_unit * 10
    telo_array = f"{rev_comp(default_query_array)}TT{default_query_array}"
    telo_array_record = f">{chrom_name}\n{telo_array}\n"
    search_region = [Interval(chrom_name)]
    targeted_region = [Interval(chrom_name, 0, len(default_query_array))]

    def test_find_all_occs_no_hits(self):
        fasta = self.make_fasta(self.telo_array_record)
        result = find_all_occurrences_in_genome(
            "AATTTTTTAAA", fasta, self.search_region, interval_window_size=0
        )
        assert result == []

    def test_find_all_occs_hits_whole_chrom(self):
        fasta = self.make_fasta(self.telo_array_record)

        result = find_all_occurrences_in_genome(
            self.default_query_array,
            fasta,
            self.search_region,
            interval_window_size=0,
        )
        expected = [
            # Forward hit
            Interval(
                self.chrom_name,
                len(self.default_query_array) + 2,
                len(self.telo_array) - 1,
            ),
            # Reverse hit
            Interval(self.chrom_name, 0, len(self.default_query_array) - 1),
        ]
        assert result == expected

    def test_find_all_occs_hits_in_region(self):
        fasta = self.make_fasta(self.telo_array_record)
        result = find_all_occurrences_in_genome(
            self.default_query_array,
            fasta,
            self.targeted_region,
            interval_window_size=0,
        )
        expected = [
            Interval(
                self.chrom_name,
                self.targeted_region[0].start,
                self.targeted_region[0].end - 1,
            )
        ]
        assert result == expected

    def test_find_all_occs_hits_with_window(self):
        fasta = self.make_fasta(self.telo_array_record)
        overflowing_window_size = 400
        result = find_all_occurrences_in_genome(
            self.default_query_array,
            fasta,
            self.targeted_region,
            interval_window_size=overflowing_window_size,
        )
        expected_end = len(self.telo_array) - 1
        expected = [Interval(self.chrom_name, 0, expected_end)]
        assert result == expected


class TestFindOccurrencesInGenomeImperfectTeloArrays(ClassWithTempFasta):
    chrom_name = "chr1"
    telo_unit = "TTAGGC"
    mutated_telo_unit = "TTTGGC"
    telo_array = f"{telo_unit*3}{mutated_telo_unit*2}{telo_unit*3}"
    telo_array_record = f">{chrom_name}\n{telo_array}\n"
    search_region = [Interval(chrom_name)]

    def test_contiguous_units_are_clustered(self):
        fasta = self.make_fasta(self.telo_array_record)

        result = find_all_occurrences_in_genome(
            self.telo_unit,
            fasta,
            self.search_region,
            interval_window_size=0,
        )
        expected = [
            Interval(self.chrom_name, 0, len(self.telo_unit) * 3 - 1),
            Interval(
                self.chrom_name, len(self.telo_unit) * 5, len(self.telo_unit) * 8 - 1
            ),
        ]
        assert result == expected

    def test_find_separate_arrays(self):
        fasta = self.make_fasta(self.telo_array_record)

        result = find_all_occurrences_in_genome(
            self.telo_unit * 3,
            fasta,
            self.search_region,
            interval_window_size=0,
        )

        expected = [
            Interval(self.chrom_name, 0, len(self.telo_unit) * 3 - 1),
            Interval(
                self.chrom_name, len(self.telo_unit) * 5, len(self.telo_unit) * 8 - 1
            ),
        ]
        assert result == expected

    def test_cluster_arrays_within_window_size(self):
        fasta = self.make_fasta(self.telo_array_record)

        result = find_all_occurrences_in_genome(
            self.telo_unit * 3,
            fasta,
            self.search_region,
            interval_window_size=7,
        )

        expected = [
            Interval(self.chrom_name, 0, len(self.telo_unit) * 8 - 1),
        ]
        assert result == expected

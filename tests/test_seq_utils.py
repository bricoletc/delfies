from pathlib import Path
from tempfile import TemporaryDirectory

import pytest
from pyfastx import Fasta

from delfies.interval_utils import Interval
from delfies.seq_utils import (
    cyclic_shifts,
    find_all_occurrences_in_genome,
    randomly_substitute,
    rev_comp,
)


def test_rev_comp():
    seq1 = "AATTAACCGG"
    expected_seq1 = "CCGGTTAATT"
    assert rev_comp(seq1) == expected_seq1


class TestSequenceMutations:
    def test_not_a_nucleotide_fails(self):
        seq_to_mutate = "PPDD"
        with pytest.raises(ValueError):
            mutated_seq = randomly_substitute(seq_to_mutate)

    def test_too_many_mutations_fails(self):
        seq_to_mutate = "AATT"
        with pytest.raises(ValueError):
            mutated_seq = randomly_substitute(
                seq_to_mutate, num_mutations=len(seq_to_mutate) + 1
            )

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


class TestFindOccurrencesInGenome:
    default_chrom_name = "chr1"
    default_query = "TTAGGC"
    default_query_revcomp_array = rev_comp(default_query) * 9
    len_default_query_revcomp_array = len(default_query_revcomp_array)
    default_reference = (
        f">{default_chrom_name}\n{default_query_revcomp_array}{default_query * 10}\n"
    )
    default_search_regions = [Interval(default_chrom_name)]

    @classmethod
    def setup_class(cls):
        cls.temp_dir = TemporaryDirectory()
        cls.temp_fasta = Path(cls.temp_dir.name) / "temp.fasta"

    @classmethod
    def teardown_class(cls):
        cls.temp_dir.cleanup()

    def make_fasta(self, input_string=default_reference):
        with self.temp_fasta.open("w") as ofstream:
            ofstream.write(input_string)
        return Fasta(str(self.temp_fasta), build_index=True)

    def test_find_all_occs_no_hits(self):
        fasta = self.make_fasta()
        result = find_all_occurrences_in_genome(
            "AATTTTTTAAA", fasta, self.default_search_regions, interval_window_size=0
        )
        assert result == []

    def test_find_all_occs_hits_whole_chrom(self):
        fasta = self.make_fasta()

        # Forward only hits
        result = find_all_occurrences_in_genome(
            self.default_query * 10,
            fasta,
            self.default_search_regions,
            interval_window_size=0,
        )
        expected_forward_start = self.len_default_query_revcomp_array
        expected = [
            Interval(
                self.default_chrom_name, expected_forward_start, expected_forward_start
            )
        ]
        assert result == expected

        # Forward and reverse hits
        result = find_all_occurrences_in_genome(
            self.default_query * 9,
            fasta,
            self.default_search_regions,
            interval_window_size=0,
        )
        expected_reverse_start = self.len_default_query_revcomp_array - 1
        expected = [
            Interval(
                self.default_chrom_name, expected_forward_start, expected_forward_start
            ),
            Interval(
                self.default_chrom_name, expected_reverse_start, expected_reverse_start
            ),
        ]
        assert result == expected

    def test_find_all_occs_hits_in_region(self):
        region_start = self.len_default_query_revcomp_array
        search_regions = [
            Interval(
                self.default_chrom_name,
                start=region_start,
                end=region_start + len(self.default_query),
            )
        ]
        fasta = self.make_fasta()
        result = find_all_occurrences_in_genome(
            self.default_query, fasta, search_regions, interval_window_size=0
        )
        expected = [Interval(self.default_chrom_name, region_start, region_start)]
        assert result == expected

    def test_find_all_occs_hits_with_window_inside_genome(self):
        region_start = self.len_default_query_revcomp_array
        search_regions = [
            Interval(
                self.default_chrom_name,
                start=region_start,
                end=region_start + len(self.default_query),
            )
        ]
        fasta = self.make_fasta()
        contained_window_size = 2
        result = find_all_occurrences_in_genome(
            self.default_query,
            fasta,
            search_regions,
            interval_window_size=contained_window_size,
        )
        expected = [
            Interval(
                self.default_chrom_name,
                region_start - contained_window_size,
                region_start + contained_window_size,
            )
        ]
        assert result == expected

    def test_find_all_occs_hits_with_window_outside_genome(self):
        region_start = self.len_default_query_revcomp_array
        search_regions = [
            Interval(
                self.default_chrom_name,
                start=region_start,
                end=region_start + len(self.default_query),
            )
        ]
        fasta = self.make_fasta()
        overflowing_window_size = 400
        result = find_all_occurrences_in_genome(
            self.default_query,
            fasta,
            search_regions,
            interval_window_size=overflowing_window_size,
        )
        expected_end = len(self.default_query) * 19 - 1
        expected = [Interval(self.default_chrom_name, 0, expected_end)]
        assert result == expected

import pytest

from delfies.interval_utils import get_contiguous_ranges, parse_region_string


class TestRegionStrings:
    def test_invalid_region_strings_fail(self):
        with pytest.raises(ValueError):
            parse_region_string("chr1:2--200")
        with pytest.raises(ValueError):
            parse_region_string("chr1::2-200")

    def test_valid_region_string(self):
        contig, start, stop = parse_region_string("chr1:2-200")
        assert contig == "chr1"
        assert start == 2
        assert stop == 200


class TestContiguousRanges:
    def test_get_repeated_range_for_single_input_num(self):
        input_nums = set([1])
        expected = [(1, 1)]
        assert get_contiguous_ranges(input_nums) == expected

    def test_contiguous_ranges_from_contiguous_inputs(self):
        input_nums = set([1, 2, 3, 4, 3, 2, 1])
        expected = [(1, 4)]
        assert get_contiguous_ranges(input_nums) == expected

    def test_contiguous_ranges_from_noncontiguous_inputs(self):
        input_nums = set([1, 2, 8, 7, 6, 5])
        expected = [(1, 2), (5, 8)]
        assert get_contiguous_ranges(input_nums) == expected

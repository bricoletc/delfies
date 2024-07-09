import pytest
from pybedtools import Interval as pybedtools_Interval

from delfies.interval_utils import Interval, get_contiguous_ranges, parse_region_string

EXAMPLE_INVALID_REGION_STRING = "chr1:2--200"
EXAMPLE_VALID_REGION_STRING = "chr1:2-200"


class TestRegionStrings:
    def test_invalid_region_strings_fail(self):
        with pytest.raises(ValueError):
            parse_region_string(EXAMPLE_INVALID_REGION_STRING)
        with pytest.raises(ValueError):
            parse_region_string("chr1::2-200")

    def test_valid_region_string(self):
        contig, start, stop = parse_region_string(EXAMPLE_VALID_REGION_STRING)
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


@pytest.fixture
def test_interval():
    return Interval("chr1", 2, 200)


class TestIntervalConstruction:
    def test_build_interval_from_invalid_region_string_fails(self):
        with pytest.raises(ValueError):
            Interval.from_region_string(EXAMPLE_INVALID_REGION_STRING)

    def test_build_interval_from_valid_region_string_succeeds(self, test_interval):
        result = Interval.from_region_string(EXAMPLE_VALID_REGION_STRING)
        assert result == test_interval

    def test_build_interval_from_pybedtools_interval(self, test_interval):
        pb_interval = pybedtools_Interval("chr1", 2, 200)
        result = Interval.from_pybedtools_interval(pb_interval)
        assert result == test_interval

    def test_interval_to_region_string(self, test_interval):
        assert test_interval.to_region_string() == EXAMPLE_VALID_REGION_STRING
        test_interval.start = None
        assert test_interval.to_region_string() == "chr1"


class TestIntervalCoordinates:
    def test_interval_with_and_without_coordinates(self, test_interval):
        assert test_interval.has_coordinates()
        test_interval.start = None
        assert not test_interval.has_coordinates()
        test_interval.start = 2
        test_interval.end = None
        assert not test_interval.has_coordinates()
        test_interval.start = None
        assert not test_interval.has_coordinates()

    def test_interval_spanning_invalid_coordinates(self, test_interval):
        with pytest.raises(ValueError):
            test_interval.start = None
            test_interval.spans(2)

    def test_interval_spanning_valid_coordinates(self, test_interval):
        assert test_interval.spans(test_interval.start)
        assert test_interval.spans(test_interval.start + 1)
        assert test_interval.spans(test_interval.end)
        assert not test_interval.spans(test_interval.start - 1)
        assert not test_interval.spans(test_interval.end + 1)

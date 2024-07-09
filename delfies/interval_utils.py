from dataclasses import dataclass
from itertools import groupby
from operator import itemgetter
from typing import List, Set, Tuple

from delfies import REGION_DELIM1, REGION_DELIM2


@dataclass
class Interval:
    name: str
    start: int = None
    end: int = None

    def spans(self, query_pos: int) -> bool:
        if not self.has_coordinates():
            raise ValueError("Interval object has not been assigned coordinates")
        return self.start <= query_pos <= self.end

    def to_region_string(self) -> str:
        result = f"{self.name}"
        if self.has_coordinates():
            result += f"{REGION_DELIM1}{self.start}{REGION_DELIM2}{self.end}"
        return result

    @classmethod
    def from_region_string(cls, region_string: str) -> "Interval":
        contig, start, end = parse_region_string(region_string)
        return Interval(contig, start, end)

    @classmethod
    def from_pybedtools_interval(cls, pybedtools_interval) -> "Interval":
        return Interval(
            pybedtools_interval.chrom,
            pybedtools_interval.start,
            pybedtools_interval.end,
        )

    def has_coordinates(self) -> bool:
        return self.start is not None and self.end is not None


Intervals = List[Interval]


def parse_region_string(region_string: str) -> Tuple[str, int, int]:
    contig, regs = region_string.split(REGION_DELIM1)
    start, stop = map(lambda e: int(e.replace(",", "")), regs.split(REGION_DELIM2))
    return contig, start, stop


def get_contiguous_ranges(input_nums: Set[int]) -> List[Tuple[int, int]]:
    """
    Credit: https://stackoverflow.com/a/2154437/12519542
    """
    result = []

    for k, g in groupby(enumerate(sorted(input_nums)), lambda x: x[0] - x[1]):
        group = map(itemgetter(1), g)
        group = list(map(int, group))
        result.append((group[0], group[-1]))
    return result

from itertools import groupby
from operator import itemgetter
from typing import List, Set, Tuple

from delfies import REGION_DELIM1, REGION_DELIM2


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

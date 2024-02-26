from itertools import groupby
from operator import itemgetter

from typing import Set, List, Tuple

def get_contiguous_ranges(input_nums: Set[int]) -> List[Tuple[int, int]]:
    """
    Credit: https://stackoverflow.com/a/2154437/12519542
    """
    result =[]

    for k,g in groupby(enumerate(input_nums),lambda x:x[0]-x[1]):
        group = (map(itemgetter(1),g))
        group = list(map(int,group))
        result.append((group[0],group[-1]))
    return result

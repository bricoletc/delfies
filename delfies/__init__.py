from importlib import metadata
from typing import Tuple

__version__ = metadata.version("delfies")

ID_DELIM = "__"
REGION_DELIM1 = ":"
REGION_DELIM2 = "-"
REGION_CLICK_HELP = (
    f"Region to focus on, format: 'contig{REGION_DELIM1}start{REGION_DELIM2}stop'."
)


def parse_region_string(region_string: str) -> Tuple[str, int, int]:
    contig, regs = region_string.split(REGION_DELIM1)
    start, stop = map(lambda e: int(e.replace(",", "")), regs.split(REGION_DELIM2))
    return contig, start, stop

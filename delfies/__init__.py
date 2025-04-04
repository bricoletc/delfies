from enum import Enum
from importlib import metadata
from dataclasses import dataclass

__version__ = metadata.version("delfies")

ID_DELIM = "__"
REGION_DELIM1 = ":"
REGION_DELIM2 = "-"
REGION_CLICK_HELP = (
    f"Region to focus on (format: 'contig{REGION_DELIM1}start{REGION_DELIM2}stop')"
)


class BreakpointType(Enum):
    """
    S2G: 'soma-to-germline': telomere-containing softclips aligned to non-telomere-containing genome region
    G2S: 'germline-to-soma': non-telomere-containing softclips aligned to telomere-containing genome region
    """

    S2G = "S2G"
    G2S = "G2S"

    def __str__(self):
        return f"{self.name}"


all_breakpoint_types = list(BreakpointType)


@dataclass
class BreakpointDetectionParams:
    bam_fname: str
    telomere_seqs: dict
    telo_array_size: int
    max_edit_distance: int
    clustering_threshold: int
    min_mapq: int
    read_filter_flag: int
    min_supporting_reads: int
    breakpoint_type: str = ""
    ofname_base: str = None

import re
from dataclasses import dataclass
from enum import Enum

from pyfastx import Fasta

from delfies.interval_utils import Interval, Intervals


class Orientation(Enum):
    forward = 1
    reverse = 2


@dataclass
class FastaRecord:
    ID: str
    sequence: str

    def __repr__(self):
        return f">{self.ID}\n{self.sequence}\n"


ORIENTATIONS = list(map(lambda e: e.name, Orientation))

REVCOMP_TABLE_DNA = dict(A="T", C="G", G="C", T="A", N="N")


def rev_comp(seq: str) -> str:
    result = "".join([REVCOMP_TABLE_DNA[elem] for elem in seq[::-1].upper()])
    return result


TELOMERE_SEQS = {
    "Nematoda": {Orientation.forward: "TTAGGC", Orientation.reverse: "GCCTAA"}
}


def cyclic_shifts(input_str: str):
    result = list()
    for i in range(len(input_str)):
        result.append(input_str[i:] + input_str[:i])
    return result


def find_all_occurrences_in_genome(
    query_sequence: str,
    genome_fasta: Fasta,
    seq_regions: Intervals,
    interval_window_size: int,
) -> Intervals:
    result = list()
    patterns = {
        Orientation.forward: re.compile(query_sequence),
        Orientation.reverse: re.compile(rev_comp(query_sequence)),
    }
    for seq_region in seq_regions:
        if seq_region.has_coordinates():
            relative_to_absolute = seq_region.start
            target_seq = genome_fasta[seq_region.name][
                seq_region.start : seq_region.end
            ]
        else:
            relative_to_absolute = 0
            target_seq = genome_fasta[seq_region.name]
        for orientation, pattern in patterns.items():
            for match in pattern.finditer(str(target_seq)):
                if orientation is Orientation.forward:
                    interval_midpoint = match.start()
                else:
                    interval_midpoint = match.end()
                interval_midpoint += relative_to_absolute
                new_interval = Interval(
                    name=seq_region.name,
                    start=interval_midpoint - interval_window_size,
                    end=interval_midpoint + interval_window_size,
                )
                result.append(new_interval)
    return result

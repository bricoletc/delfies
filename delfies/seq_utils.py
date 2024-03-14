from dataclasses import dataclass
from enum import Enum


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
    result = "".join([REVCOMP_TABLE_DNA[elem] for elem in seq[::-1]])
    return result


TELOMERE_SEQS = {
    "Nematoda": {Orientation.forward: "TTAGGC", Orientation.reverse: "GCCTAA"}
}


def cyclic_shifts(input_str: str):
    result = list()
    for i in range(len(input_str)):
        result.append(input_str[i:] + input_str[:i])
    return result

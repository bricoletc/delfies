from enum import Enum


class Orientation(Enum):
    forward = 1
    reverse = 2

    def from_str(string: str):
        if "5prime" in string:
            return Orientation.reverse
        elif "3prime" in string:
            return Orientation.forward
        else:
            raise ValueError(f"{string} has no '[53]prime' indication")


ORIENTATIONS = list(map(lambda e: e.name, Orientation))

REVCOMP_TABLE_DNA = dict(A="T", C="G", G="C", T="A", N="N")


def rev_comp(seq: str) -> str:
    result = "".join([REVCOMP_TABLE_DNA[elem] for elem in seq[::-1]])
    return result


TELOMERE_SEQS = {
    "Nematoda": {Orientation.forward: "TTAGGC", Orientation.reverse: "GCCTAA"}
}

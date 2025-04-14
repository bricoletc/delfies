import re
from dataclasses import dataclass
from enum import Enum
from random import choice as random_choice

from pyfastx import Fasta

from delfies.interval_utils import Interval, Intervals


class Orientation(Enum):
    forward = "+"
    reverse = "-"


@dataclass
class FastaRecord:
    ID: str
    sequence: str

    def __repr__(self):
        return f">{self.ID}\n{self.sequence}\n"


ORIENTATIONS = list(map(lambda e: e.name, Orientation))

REVCOMP_TABLE_DNA = dict(A="T", C="G", G="C", T="A", N="N")
NUCLEOTIDES = set(REVCOMP_TABLE_DNA.keys())


def rev_comp(seq: str) -> str:
    result = "".join([REVCOMP_TABLE_DNA[elem] for elem in seq[::-1].upper()])
    return result


def randomly_substitute(seq: str, num_mutations: int = 1) -> str:
    if num_mutations > len(seq):
        raise ValueError(
            f"{num_mutations} is greater than the length of the input sequence"
        )
    index_choices = list(range(len(seq)))
    mutated_seq = seq
    while num_mutations != 0:
        index_choice = random_choice(index_choices)
        index_choices = [el for el in index_choices if el != index_choice]
        chosen_nucleotide = seq[index_choice]
        if chosen_nucleotide not in NUCLEOTIDES:
            raise ValueError("Not a nucleotide: {chosen_nucleotide}")
        possible_mutations = [el for el in NUCLEOTIDES if el != chosen_nucleotide]
        chosen_mutation = random_choice(possible_mutations)
        mutated_seq = (
            mutated_seq[:index_choice]
            + chosen_mutation
            + mutated_seq[index_choice + 1 :]
        )
        num_mutations -= 1
    return mutated_seq


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
        chrom_length = len(genome_fasta[seq_region.name])
        for orientation, pattern in patterns.items():
            match = None
            for match in pattern.finditer(str(target_seq)):
                new_interval = Interval(
                    name=seq_region.name,
                    start=max(
                        0, match.start() + relative_to_absolute - interval_window_size
                    ),
                    end=min(
                        chrom_length - 1,
                        match.end() + relative_to_absolute - 1 + interval_window_size,
                    ),
                )
                if len(result) == 0:
                    result.append(new_interval)
                else:
                    if new_interval.overlaps_or_touches(result[-1]):
                        result[-1].end = new_interval.end
                    else:
                        result.append(new_interval)
            if match is not None and (
                len(result) == 0 or not new_interval.overlaps_or_touches(result[-1])
            ):
                result.append(new_interval)
    return result

from typing import List

from pyfastx import Fasta

from delfies.breakpoint_foci import READ_SUPPORT_PREFIX, MaximalFoci
from delfies.seq_utils import FastaRecord, Orientation, rev_comp


def extract_breakpoint_sequences(
    maximal_foci: MaximalFoci, genome: Fasta, seq_window_size: int
) -> List[FastaRecord]:
    result = list()
    for max_focus in maximal_foci:
        focus = max_focus.focus
        breakpoint_pos = max(int(focus.start), 0)
        windowed_start = max(breakpoint_pos - seq_window_size + 1, 0)
        windowed_end = breakpoint_pos + seq_window_size
        breakpoint_sequence = (
            genome.fetch(focus.contig, (windowed_start, breakpoint_pos))
            + "N"
            + genome.fetch(focus.contig, (breakpoint_pos + 1, windowed_end))
        )
        strand_name = "3prime"
        if max_focus.orientation is Orientation.reverse:
            breakpoint_sequence = rev_comp(breakpoint_sequence)
            strand_name = "5prime"
        breakpoint_name = f"{max_focus.breakpoint_type}_{strand_name}_{focus.contig} breakpoint_pos:{breakpoint_pos} {READ_SUPPORT_PREFIX}:{max_focus.max_value} next_best_value_on_same_strand:{max_focus.next_max_value} best_value_on_other_strand:{max_focus.max_value_other_orientation}"
        result.append(FastaRecord(breakpoint_name, breakpoint_sequence))
    return result


def write_breakpoint_sequences(
    genome_fname: str, maximal_foci: MaximalFoci, odirname: str, seq_window_size: int
) -> None:
    genome = Fasta(genome_fname)
    breakpoint_sequences = extract_breakpoint_sequences(
        maximal_foci, genome, seq_window_size
    )
    breakpoint_fasta = odirname / "breakpoint_sequences.fasta"
    with breakpoint_fasta.open("w") as ofstream:
        for breakpoint_sequence in breakpoint_sequences:
            ofstream.write(str(breakpoint_sequence))

from delfies.breakpoint_foci import MaximalFocus, setup_tents
from delfies.breakpoint_sequences import extract_breakpoint_sequences
from delfies.seq_utils import Orientation
from tests import ClassWithTempFasta


class TestExtractBreakpointSequences(ClassWithTempFasta):
    default_reference = f">scaffold_1\nACGTGATACA\n"
    tents = setup_tents()
    default_max_focus = MaximalFocus(Orientation.forward, 0, 0, 0, (0, 0), None)
    breakpoint_location = tents.new()
    breakpoint_location.update(contig="scaffold_1", start=3)
    default_max_focus.focus = breakpoint_location

    def test_extract_breakpoint_within_boundaries_single_nucleotide(self):
        fasta = self.make_fasta(self.default_reference)
        self.breakpoint_location.start = 3
        breakpoint_sequences = extract_breakpoint_sequences(
            [self.default_max_focus], fasta, seq_window_size=1
        )
        assert len(breakpoint_sequences) == 1
        assert breakpoint_sequences[0].sequence == "GNT"

    def test_extract_breakpoint_within_boundaries_multiple_nucleotides(self):
        fasta = self.make_fasta(self.default_reference)
        self.breakpoint_location.start = 3
        breakpoint_sequences = extract_breakpoint_sequences(
            [self.default_max_focus], fasta, seq_window_size=2
        )
        assert len(breakpoint_sequences) == 1
        assert breakpoint_sequences[0].sequence == "CGNTG"

    def test_extract_breakpoint_outside_boundaries_negative_start(self):
        fasta = self.make_fasta(self.default_reference)
        for start in [-2, 0]:
            self.breakpoint_location.start = start
            breakpoint_sequences = extract_breakpoint_sequences(
                [self.default_max_focus], fasta, seq_window_size=2
            )
            assert len(breakpoint_sequences) == 1
            assert breakpoint_sequences[0].sequence == "NAC"

    def test_extract_breakpoint_outside_boundaries_end_past_last_nucleotide(self):
        fasta = self.make_fasta(self.default_reference)
        self.breakpoint_location.start = 9
        breakpoint_sequences = extract_breakpoint_sequences(
            [self.default_max_focus], fasta, seq_window_size=3
        )
        assert len(breakpoint_sequences) == 1
        assert breakpoint_sequences[0].sequence == "TACNA"

    def test_extract_breakpoint_within_boundaries_reverse_orientation(self):
        fasta = self.make_fasta(self.default_reference)
        self.breakpoint_location.start = 3
        self.default_max_focus.orientation = Orientation.reverse
        breakpoint_sequences = extract_breakpoint_sequences(
            [self.default_max_focus], fasta, seq_window_size=2
        )
        assert len(breakpoint_sequences) == 1
        assert breakpoint_sequences[0].sequence == "CANCG"
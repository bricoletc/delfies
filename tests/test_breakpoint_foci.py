import pytest
from pysam import AlignedSegment, AlignmentHeader

from delfies import BreakpointDetectionParams, BreakpointType, PutativeBreakpoint
from delfies.breakpoint_foci import (
    READ_SUPPORTS,
    FociWindow,
    Orientation,
    cluster_breakpoint_foci,
    record_softclips,
    setup_breakpoint_tents,
)
from delfies.interval_utils import Interval
from delfies.SAM_utils import DEFAULT_MIN_MAPQ, DEFAULT_READ_FILTER_FLAG
from delfies.seq_utils import rev_comp
from delfies.telomere_utils import TELOMERE_SEQS


@pytest.fixture
def breakpoint_focus():
    tents = setup_breakpoint_tents()
    new_tent = tents.new()
    new_tent.update(contig="test_contig", start=2, end=200)
    new_tent[f"{READ_SUPPORTS[0]}"] = 15
    new_tent[f"{READ_SUPPORTS[1]}"] = 20
    return new_tent


@pytest.fixture
def multiple_breakpoint_foci(breakpoint_focus):
    tents = setup_breakpoint_tents()
    tents.add(breakpoint_focus)
    new_tent = tents.new()
    new_tent.update(contig="test_contig", start=0, end=202)
    new_tent[f"{READ_SUPPORTS[0]}"] = 15
    tents.add(new_tent)
    new_tent = tents.new()
    new_tent.update(contig="test_contig", start=2000, end=2005)
    new_tent[f"{READ_SUPPORTS[0]}"] = 15
    tents.add(new_tent)
    return tents


@pytest.fixture
def focus_window():
    tents = setup_breakpoint_tents()
    new_tent = tents.new()
    new_tent.update(start=205, end=210)
    new_tent[f"{READ_SUPPORTS[0]}"] = 2
    new_tent[f"{READ_SUPPORTS[1]}"] = 200
    return FociWindow(new_tent)


@pytest.fixture
def putative_breakpoint():
    max_focus = PutativeBreakpoint(
        orientation=Orientation.forward,
        max_value=10,
        next_max_value=2,
        max_value_other_orientation=1,
        interval=(205, 210),
        focus=None,
    )
    return max_focus


class TestPutativeBreakpoint:
    def test_putative_breakpoint_update_no_new_max(
        self, breakpoint_focus, putative_breakpoint
    ):
        putative_breakpoint.max_value = 100
        putative_breakpoint.update(breakpoint_focus)
        assert putative_breakpoint.focus is None
        assert putative_breakpoint.max_value == 100
        assert putative_breakpoint.next_max_value == 15

    def test_putative_breakpoint_update_new_max(
        self, breakpoint_focus, putative_breakpoint
    ):
        prev_max_value = putative_breakpoint.max_value
        putative_breakpoint.update(breakpoint_focus)
        assert putative_breakpoint.focus == breakpoint_focus
        assert putative_breakpoint.max_value == breakpoint_focus[f"{READ_SUPPORTS[0]}"]
        assert putative_breakpoint.next_max_value == prev_max_value


class TestFociWindow:
    def test_inclusion_inside_or_overlapping_or_spanning(
        self, breakpoint_focus, focus_window
    ):
        # Overlapping
        fw_start = focus_window.foci[0].start
        fw_end = focus_window.foci[0].end
        breakpoint_focus.start = fw_start + 1
        breakpoint_focus.end = fw_end + 1
        assert focus_window.includes(breakpoint_focus, tolerance=0)
        breakpoint_focus.start = fw_start - 1
        breakpoint_focus.end = fw_end - 1
        # Inside
        assert focus_window.includes(breakpoint_focus, tolerance=0)
        breakpoint_focus.start = fw_start + 1
        breakpoint_focus.end = fw_end - 1
        assert focus_window.includes(breakpoint_focus, tolerance=0)
        # Spanning
        breakpoint_focus.start = fw_start - 1
        breakpoint_focus.end = fw_end + 1
        assert focus_window.includes(breakpoint_focus, tolerance=0)

    def test_no_inclusion_outside(self, breakpoint_focus, focus_window):
        fw_start = focus_window.foci[0].start
        fw_end = focus_window.foci[0].end
        breakpoint_focus.start = fw_start - 10
        breakpoint_focus.end = fw_start - 5
        assert not focus_window.includes(breakpoint_focus, tolerance=0)
        breakpoint_focus.start = fw_end + 10
        breakpoint_focus.end = fw_end + 15
        assert not focus_window.includes(breakpoint_focus, tolerance=0)

    def test_inclusion_with_tolerance(self, breakpoint_focus, focus_window):
        fw_start = focus_window.foci[0].start
        breakpoint_focus.start = fw_start - 2
        breakpoint_focus.end = fw_start - 1
        assert not focus_window.includes(breakpoint_focus, tolerance=0)
        assert focus_window.includes(breakpoint_focus, tolerance=1)

    def test_add(self, breakpoint_focus, focus_window):
        fw_start = focus_window.foci[0].start
        fw_end = focus_window.foci[0].end
        breakpoint_focus.start = fw_start + 1
        breakpoint_focus.end = fw_end + 2
        focus_window.add(breakpoint_focus)
        assert len(focus_window.foci) == 2
        assert focus_window.foci[1] == breakpoint_focus
        assert focus_window.Min == fw_start
        assert focus_window.Max == fw_end + 2

    def test_find_peak_softclip_focus_reverse_max(self, breakpoint_focus, focus_window):
        focus_window.add(breakpoint_focus)
        max_focus = focus_window.find_peak_softclip_focus()
        assert max_focus == PutativeBreakpoint(
            orientation=Orientation.reverse,
            max_value=200,
            next_max_value=20,
            max_value_other_orientation=15,
            interval=(focus_window.Min, focus_window.Max),
            focus=focus_window.foci[0],
        )

    def test_find_peak_softclip_focus_forward_max(self, breakpoint_focus, focus_window):
        breakpoint_focus[f"{READ_SUPPORTS[0]}"] = 400
        focus_window.add(breakpoint_focus)
        max_focus = focus_window.find_peak_softclip_focus()
        assert max_focus == PutativeBreakpoint(
            orientation=Orientation.forward,
            max_value=400,
            next_max_value=2,
            max_value_other_orientation=200,
            interval=(focus_window.Min, focus_window.Max),
            focus=breakpoint_focus,
        )


class TestBreakpointClustering:
    def test_breakpoint_foci_with_no_read_support_are_filtered_out(
        self, multiple_breakpoint_foci
    ):
        for focus in multiple_breakpoint_foci:
            focus[READ_SUPPORTS[0]] = 0
            focus[READ_SUPPORTS[1]] = 0
        result = cluster_breakpoint_foci(multiple_breakpoint_foci, tolerance=10)
        assert len(result) == 0

    def test_cluster_breakpoint_foci(self, multiple_breakpoint_foci):
        result = cluster_breakpoint_foci(multiple_breakpoint_foci, tolerance=10)
        assert len(result) == 2
        assert len(result[0].foci) == 2
        assert result[0].Min == 0
        assert result[0].Max == 202
        assert len(result[1].foci) == 1
        assert result[1].Min == 2000
        assert result[1].Max == 2005


DEFAULT_TELO_SEQ = TELOMERE_SEQS["Nematoda"][Orientation.forward]
DEFAULT_TELO_ARRAY_SIZE = 3
DEFAULT_NON_TELO_SEQ = "TAACCC"
DEFAULT_CHROM = "chr1"


@pytest.fixture
def read_telo_seq_forward_3prime():
    non_telo_seq = DEFAULT_NON_TELO_SEQ * 3
    telo_seq = DEFAULT_TELO_SEQ * DEFAULT_TELO_ARRAY_SIZE
    sequence = non_telo_seq + telo_seq
    cigar = ((0, len(non_telo_seq)), (4, len(telo_seq)))
    return (sequence, cigar)


@pytest.fixture
def read_telo_seq_reverse_3prime():
    non_telo_seq = DEFAULT_NON_TELO_SEQ * 3
    telo_seq = rev_comp(DEFAULT_TELO_SEQ) * DEFAULT_TELO_ARRAY_SIZE
    sequence = non_telo_seq + telo_seq
    cigar = ((0, len(non_telo_seq)), (4, len(telo_seq)))
    return (sequence, cigar)


@pytest.fixture
def softclipped_aligned_read():
    header = AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.0"},
            "SQ": [{"LN": 1e03, "SN": DEFAULT_CHROM}],
        }
    )
    read = AlignedSegment(header=header)
    read.reference_id = 0
    read.query_name = "simulated_read"
    read.flag = 0
    read.mapping_quality = DEFAULT_MIN_MAPQ
    read.reference_start = 3
    return read


@pytest.fixture
def detection_params():
    telomere_seqs = {
        Orientation.forward: DEFAULT_TELO_SEQ,
        Orientation.reverse: rev_comp(DEFAULT_TELO_SEQ),
    }
    return BreakpointDetectionParams(
        bam_fname="NA",
        telomere_seqs=telomere_seqs,
        telo_array_size=DEFAULT_TELO_ARRAY_SIZE,
        max_edit_distance=0,
        clustering_threshold=1,
        min_mapq=DEFAULT_MIN_MAPQ,
        read_filter_flag=DEFAULT_READ_FILTER_FLAG,
        min_supporting_reads=1,
        keep_telomeric_breakpoints=False,
        breakpoint_type=BreakpointType.G2S,
    )


class TestRecordSoftclips:
    def test_G2S_read_with_3prime_forward_telo_softclips_is_rejected(
        self, detection_params, read_telo_seq_forward_3prime, softclipped_aligned_read
    ):
        softclipped_aligned_read.query_sequence = read_telo_seq_forward_3prime[0]
        softclipped_aligned_read.cigar = read_telo_seq_forward_3prime[1]
        breakpoint_tents = setup_breakpoint_tents()
        seq_region = Interval(DEFAULT_CHROM, 1, 100)
        record_softclips(
            softclipped_aligned_read,
            breakpoint_tents,
            dict(),
            detection_params,
            seq_region,
        )
        assert len(breakpoint_tents) == 0

    def test_G2S_read_with_3prime_reverse_telo_softclips_is_rejected(
        self, detection_params, read_telo_seq_reverse_3prime, softclipped_aligned_read
    ):
        softclipped_aligned_read.query_sequence = read_telo_seq_reverse_3prime[0]
        softclipped_aligned_read.cigar = read_telo_seq_reverse_3prime[1]
        breakpoint_tents = setup_breakpoint_tents()
        seq_region = Interval(DEFAULT_CHROM, 1, 100)
        record_softclips(
            softclipped_aligned_read,
            breakpoint_tents,
            dict(),
            detection_params,
            seq_region,
        )
        assert len(breakpoint_tents) == 0

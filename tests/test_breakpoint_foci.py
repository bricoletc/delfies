import pytest

from delfies import ID_DELIM
from delfies.breakpoint_foci import (
    READ_SUPPORTS,
    FociWindow,
    MaximalFocus,
    Orientation,
    cluster_breakpoint_foci,
    setup_tents,
)


@pytest.fixture
def breakpoint_focus():
    tents = setup_tents()
    new_tent = tents.new()
    new_tent.update(contig="test_contig", start=2, end=200)
    new_tent[f"{READ_SUPPORTS[0]}"] = 15
    new_tent[f"{READ_SUPPORTS[1]}"] = 20
    return new_tent


@pytest.fixture
def multiple_breakpoint_foci(breakpoint_focus):
    tents = setup_tents()
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
    tents = setup_tents()
    new_tent = tents.new()
    new_tent.update(start=205, end=210)
    new_tent[f"{READ_SUPPORTS[0]}"] = 2
    new_tent[f"{READ_SUPPORTS[1]}"] = 200
    return FociWindow(new_tent)


@pytest.fixture
def maximal_focus():
    max_focus = MaximalFocus(
        orientation=Orientation.forward,
        max_value=10,
        next_max_value=2,
        max_value_other_orientation=1,
        interval=(205, 210),
        focus=None,
    )
    return max_focus


class TestMaximalFocus:
    def test_maximal_focus_update_no_new_max(self, breakpoint_focus, maximal_focus):
        maximal_focus.max_value = 100
        maximal_focus.update(breakpoint_focus)
        assert maximal_focus.focus is None
        assert maximal_focus.max_value == 100
        assert maximal_focus.next_max_value == 15

    def test_maximal_focus_update_new_max(self, breakpoint_focus, maximal_focus):
        prev_max_value = maximal_focus.max_value
        maximal_focus.update(breakpoint_focus)
        assert maximal_focus.focus == breakpoint_focus
        assert maximal_focus.max_value == breakpoint_focus[f"{READ_SUPPORTS[0]}"]
        assert maximal_focus.next_max_value == prev_max_value


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
        fw_end = focus_window.foci[0].end
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
        assert max_focus == MaximalFocus(
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
        assert max_focus == MaximalFocus(
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

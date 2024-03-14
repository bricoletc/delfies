from delfies.seq_utils import rev_comp, cyclic_shifts


def test_rev_comp():
    seq1 = "AATTAACCGG"
    expected_seq1 = "CCGGTTAATT"
    assert rev_comp(seq1) == expected_seq1


def test_cyclic_shifts():
    str_to_shift = "TTAGGC"
    expected_shifts = [
            "TTAGGC",
            "TAGGCT",
            "AGGCTT",
            "GGCTTA",
            "GCTTAG",
            "CTTAGG",
    ]
    result = cyclic_shifts(str_to_shift)
    assert result == expected_shifts

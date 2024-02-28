from delfies.seq_utils import rev_comp

def test_rev_comp():
    seq1 = "AATTAACCGG"
    expected_seq1 = "CCGGTTAATT"
    assert rev_comp(seq1) == expected_seq1

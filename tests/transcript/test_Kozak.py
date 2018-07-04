"""
Test get_Kozak_seq, get_consensus_Kozak_seq and get_line_Kozak_matrix methods
"""

from tests.transcript.config import genome


def test_get_consensus_Kozak_seq_numbers(genome):
    res = genome.get_consensus_Kozak_seq()
    assert res == "222211032212212"


def test_get_consensus_Kozak_seq_sequence(genome):
    res = genome.get_consensus_Kozak_seq(True)
    assert res == "GGGGCCATGGCGGCG"

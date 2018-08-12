"""
Test get_stop_codon_context, get_consensus_stop_codon_context
"""

from tests.transcript.config import genome, transcript2

import pandas as pd
import pandas.util.testing as pdt


def test_get_stop_codon_context(transcript2):
    res = transcript2.get_stop_codon_context()
    assert res == "TTTATTTCCACAAAT"


def test_get_consensus_stop_codon_context_numbers(genome):
    res = genome.get_consensus_stop_codon_context()
    assert res == "122222222222112"


def test_get_consensus_stop_codon_context(genome):
    res = genome.get_consensus_stop_codon_context(True)
    assert res == "CGGGGGGGGGGGCCG"

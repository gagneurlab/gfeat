"""
Test get_stop_codon_context, get_consensus_stop_codon_context and get_line_stop_codon_context_matrix
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


def test_get_line_stop_codon_context_matrix(transcript2):
    res = transcript2.get_stop_codon_context_as_df()
    dict_stop_codon_context = {"0A": [0], "0C": [0], "0G": [0], "0T": [1],
                            "1A": [0], "1C": [0], "1G": [0], "1T": [1],
                            "2A": [0], "2C": [0], "2G": [0], "2T": [1],
                            "3A": [1], "3C": [0], "3G": [0], "3T": [0],
                            "4A": [0], "4C": [0], "4G": [0], "4T": [1],
                            "5A": [0], "5C": [0], "5G": [0], "5T": [1],
                            "6A": [0], "6C": [0], "6G": [0], "6T": [1],
                            "7A": [0], "7C": [1], "7G": [0], "7T": [0],
                            "8A": [0], "8C": [1], "8G": [0], "8T": [0],
                            "9A": [1], "9C": [0], "9G": [0], "9T": [0],
                            "10A": [0], "10C": [1], "10G": [0], "10T": [0],
                            "11A": [1], "11C": [0], "11G": [0], "11T": [0],
                            "12A": [1], "12C": [0], "12G": [0], "12T": [0],
                            "13A": [1], "13C": [0], "13G": [0], "13T": [0],
                            "14A": [0], "14C": [0], "14G": [0], "14T": [1]}
    df_line_stop_codon_context = pd.DataFrame(data=dict_stop_codon_context)
    pdt.assert_frame_equal(df_line_stop_codon_context, res)

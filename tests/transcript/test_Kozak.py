"""
Test get_Kozak_seq, get_consensus_Kozak_seq and get_line_Kozak_matrix
"""

from tests.transcript.config import genome, transcript2

import pandas as pd
import pandas.util.testing as pdt


def test_get_Kosak_seq(transcript2):
    res = transcript2.get_Kozak_seq()
    assert res == "AAGCAGATGGTGGCT"


def test_get_consensus_Kozak_seq_numbers(genome):
    res = genome.get_consensus_Kozak_seq()
    assert res == "222211032212212"


def test_get_consensus_Kozak_seq(genome):
    res = genome.get_consensus_Kozak_seq(True)
    assert res == "GGGGCCATGGCGGCG"


def test_get_line_Kozak_matrix(transcript2):
    res = transcript2.get_Kozak_seq_as_df()
    dict_line_Kozak_matrix = {"0A": [1], "0C": [0], "0G": [0], "0T": [0],
                            "1A": [1], "1C": [0], "1G": [0], "1T": [0],
                            "2A": [0], "2C": [0], "2G": [1], "2T": [0],
                            "3A": [0], "3C": [1], "3G": [0], "3T": [0],
                            "4A": [1], "4C": [0], "4G": [0], "4T": [0],
                            "5A": [0], "5C": [0], "5G": [1], "5T": [0],
                            "6A": [1], "6C": [0], "6G": [0], "6T": [0],
                            "7A": [0], "7C": [0], "7G": [0], "7T": [1],
                            "8A": [0], "8C": [0], "8G": [1], "8T": [0],
                            "9A": [0], "9C": [0], "9G": [1], "9T": [0],
                            "10A": [0], "10C": [0], "10G": [0], "10T": [1],
                            "11A": [0], "11C": [0], "11G": [1], "11T": [0],
                            "12A": [0], "12C": [0], "12G": [1], "12T": [0],
                            "13A": [0], "13C": [1], "13G": [0], "13T": [0],
                            "14A": [0], "14C": [0], "14G": [0], "14T": [1]}
    df_line_Kozak_matrix = pd.DataFrame(data=dict_line_Kozak_matrix)
    pdt.assert_frame_equal(df_line_Kozak_matrix, res)

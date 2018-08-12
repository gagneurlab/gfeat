"""
Test get_codon_pairs_frequency
"""

from tests.transcript.config import transcript2
#import pandas as pd

def test_get_codon_pairs_frequency(transcript2):
    res = transcript2.get_codon_pairs_frequency()
    res_elem1 = res["AAAAAA"]
    test_elem1 = {0: 0.0018315018315018315}
    res_elem2 = res["TTTTAT"]
    test_elem2 = {0: 0.001221001221001221}
    res_elem3 = res["GCTAAG"]
    test_elem3 = {0: 0.001221001221001221}
    res_elem4 = res["AGCAAC"]
    test_elem4 = {0: 0.0018315018315018315}
    assert res_elem1[0] == test_elem1[0] and res_elem2[0] == test_elem2[0] \
           and res_elem3[0] == test_elem3[0] and res_elem4[0] == test_elem4[0] and len(res.columns) == 3600

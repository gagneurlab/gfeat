"""
Test upstreamATG.predict_on_sample
"""

# import pytest
from gfeat.upstreamAUG import UpstreamAUG
# from tests.uATG.config_uATG import uATG_object
import numpy as np


def test_predict_on_sample():
    uATG_object = UpstreamAUG()
    res = uATG_object.predict_on_sample("GTCCTAGAGCAGCTCCAGCAGGATGGCGGCTCCAGCGTCTCTAAGGCCTGCAGGGGGTCCAGCCCCATGGGGGGCGCC"
                                        "CTAGGCCTCCGACAGCTCCCCATCTGTGCTCCTGCCTGCCGGCCATCCTCAGGCCACTCGCC")
    assert np.array_equal(res, np.array([0, 0]))


def test_predict_on_sample_ORF():
    uATG_object = UpstreamAUG(True, True)
    res = uATG_object.predict_on_sample("GTCCTAGAGCAGCTCCAGCAGGATGGCGGCTCCAGCGTCTCTAAGGCCTGCAGGGGGTCCAGCCCCATGGGGGGCGCC"
                                        "CTAGGCCTCCGACAGCTCCCCATCTGTGCTCCTGCCTGCCGGCCATCCTCAGGCCACTCGCC")
    assert np.array_equal(res["frame"], np.array([0, 0])) and np.array_equal(res["uORF"], np.array([1, 0]))

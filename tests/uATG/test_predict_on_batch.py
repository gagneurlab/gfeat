"""
Test upstreamATG.predict_on_batch
"""

# import pytest
from gfeat.upstreamAUG import UpstreamAUG
from tests.uATG.config_uATG import seq_list
import numpy as np


def test_predict_on_batch(seq_list):
    uATG_object = UpstreamAUG()
    res = uATG_object.predict_on_batch(seq_list)
    assert np.array_equal(res[0], np.array([0, 0, 1])) and np.array_equal(res[1], np.array([1])) and \
           np.array_equal(res[2], np.array([0, 1, 1])) and np.array_equal(res[3], np.array([]))


def test_predict_on_batch_ORF(seq_list):
    uATG_object = UpstreamAUG(True, True)
    res = uATG_object.predict_on_batch(seq_list)
    assert np.array_equal((res[0])["frame"], np.array([0, 0, 1])) and np.array_equal((res[1])["frame"], np.array([1])) and \
           np.array_equal((res[2])["frame"], np.array([0, 1, 1])) and np.array_equal((res[3])["frame"], np.array([])) and \
           np.array_equal((res[0])["uORF"], np.array([1, 1, 1])) and np.array_equal((res[1])["uORF"], np.array([1])) and \
           np.array_equal((res[2])["uORF"], np.array([1, 0, 0])) and np.array_equal((res[3])["uORF"], np.array([]))

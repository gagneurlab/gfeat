"""
Test PCA_with_scaling
"""

from gfeat.utils import PCA_with_standard_sample_deviation_scaling
import numpy as np
import pandas as pd
import numpy.testing as npt


def test_get_codon_pairs_count():
    matrix = np.array([[62, 68, 5], [24, 94, 93], [96, 223, 23], [5, 9, 44]])
    df = pd.DataFrame(matrix)
    res = PCA_with_standard_sample_deviation_scaling(df)
    array_test = np.array([[0.47956485, -0.95669923],
                           [-1.0122183, 1.06750924],
                           [1.84616296, 0.44112221],
                           [-1.31350951, -0.55193222]])
    npt.assert_allclose(res, array_test)

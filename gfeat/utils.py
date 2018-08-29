"""
utils module
"""
from sklearn.decomposition import TruncatedSVD
import numpy as np


def reverse_complement(dna):
    """
    Reverse-complement a DNA sequence

    :param dna: string, DNA sequence
    :return: reverse-complement of a DNA sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])


def PCA_with_scaling(df, n_comp=2):
    """
    PCA with scaling. Scaling is done using sample standard deviation

    :param df: pandas.DataFrame, input matrix
    :param n_comp: number of principle components
    :return: pandas.DataFrame, truncated principal component decomposition
    """
    svd = TruncatedSVD(n_components=n_comp)
    # sklearn.preprocessing.StandardScaler uses population standard deviation (where the sum of squared deviations is
    # divided by the number of observations). In contrast, here sample standard deviations (denominator equal to number
    # of observations - 1) is used.
    manually_scaled_df = (df.values - np.mean(df.values, axis=0)) / np.std(df.values, axis=0, ddof=1)
    svd.fit(manually_scaled_df)
    return manually_scaled_df.dot(svd.components_.transpose())

"""
Test coding_sequence_codon_count
"""

from tests.gene.config_1 import gene


def test_coding_sequence_codon_count(gene):
    res = gene.coding_sequence_codon_count("AAC")
    assert res == 25

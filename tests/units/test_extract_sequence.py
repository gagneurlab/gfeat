"""
Test extract_sequence
"""

from tests.transcript.config import interval, fasta, vcf
from gfeat.units import extract_sequence

# def test_extract_sequence(interval, fasta):
#     res = test_extract_sequence(interval, fasta)
#     assert res == [("CCCTGAGTCATCCTTGC")]

def test_extract_sequence(interval, fasta, vcf):
    res = extract_sequence(interval, fasta, vcf)
    assert res == [("CCCTGAGGCGTCCTTGC", (12, 14)), ("CCCCAAGGCGTCCTTGC", (12, 14, 8)),
                   ("CCCTGAAGCGTCCTTGC", (12, 14, 11)), ("CCCTGAGGCGTCTTTGC", (12, 14, 17)),
                   ("CCCCAAAGCGTCCTTGC", (12, 14, 8, 11)), ("CCCCAAGGCGTCTTTGC", (12, 14, 8, 17)),
                   ("CCCTGAAGCGTCTTTGC", (12, 14, 11, 17)), ("CCCCAAAGCGTCTTTGC", (12, 14, 8, 11, 17))]

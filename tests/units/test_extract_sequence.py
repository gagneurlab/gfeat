"""
Test extract_sequence
"""

from config import interval, fasta, vcf

def test_extract_sequence(interval, fasta):
    res = test_extract_sequence(interval, fasta)
    assert res == [("CCCTGAGTCATCCTTGC")]

def test_extract_sequence(interval, fasta, vcf):
    res = test_extract_sequence(interval, fasta, vcf)
    assert res == [("CCCTGAGGCGTCCTTGC", (12, 14)), ("CCCCAAGGCGTCCTTGC", (8, 12, 14)),
                   ("CCCCAAAGCGTCCTTGC", (8, 11, 12, 14)), ("CCCCAAGGCGTCTTTGC", (8, 12, 14, 17)),
                   ("CCCCAAAGCGTCTTTGC", (8, 11, 12, 14, 17)), ("CCCTGAAGCGTCCTTGC", (11, 12, 14)),
                   ("CCCTGAAGCGTCTTTGC", (11, 12, 14, 17)),("CCCTGAGGCGTCTTTGC", (12, 14, 17))]

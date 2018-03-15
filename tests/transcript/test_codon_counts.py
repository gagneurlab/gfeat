"""
Test codon_counts
"""

from tests.transcript.config import transcript
# import os
# os.system("pyensembl install --reference-name hg38_test --annotation-name hg38_chr22_test --gtf \"./tests/data/gencode.v24.annotation_chr22.gtf\" --transcript-fasta \"./tests/data/hg38_chr22.fa\"")

def test_codon_counts(transcript):
    res = transcript.codon_counts()
    assert res == 2

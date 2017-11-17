"""
Test codon_counts
"""

import pytest
import pyensembl

from gfeat.transcript import GFTranscript


@pytest.fixture()
def transcript():
    data = pyensembl.ensembl_release.EnsemblRelease(75)
    test_transcript = GFTranscript("ENST00000369985", "MYO6-001", "6", 76458926, 76629253, "+", "protein_coding",
                                "ENSG00000196586", data)
    return test_transcript

def test_codon_counts(transcript):
    res = transcript.codon_counts()
    assert res == 1263

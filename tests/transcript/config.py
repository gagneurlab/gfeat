import pytest
import pyensembl
from pysam import FastaFile
from pybedtools import BedTool

from gfeat.transcript import GFTranscript


@pytest.fixture()
def transcript():
    data = pyensembl.ensembl_release.EnsemblRelease(75)
    test_transcript = GFTranscript("ENST00000369985", "MYO6-001", "6", 76458926, 76629253, "+", "protein_coding",
                                "ENSG00000196586", data)
    return test_transcript


@pytest.fixture()
def interval_plus():
    gtf = BedTool("./tests/data/chr22_testing.gtf")
    test_interval = gtf[0]  # get the first element in the bed file, returns an Interval object
    test_interval.strand = '+'
    return test_interval


@pytest.fixture()
def interval_minus():
    gtf = BedTool("./tests/data/chr22_testing.gtf")
    test_interval = gtf[0]  # get the first element in the bed file, returns an Interval object
    return test_interval


@pytest.fixture()
def fasta():
    test_fasta = FastaFile("./tests/data/chr22_testing.fa")
    return test_fasta


@pytest.fixture()
def vcf():
    test_vcf = "./tests/data/49470G_chr22_testing.vcf.gz"
    return test_vcf

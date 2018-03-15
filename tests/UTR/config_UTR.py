import pytest
import pyensembl
from pysam import FastaFile
from pybedtools import BedTool
from gfeat.genome import GFGenome
from gfeat.transcript import GFTranscript
from gfeat.UTR import FivePrimeUTRSeq


@pytest.fixture()
def FivePrimeUTRSeq_object():
    data = GFGenome(reference_name='hg38_test_FivePrimeUTRSeq',
                     annotation_name='hg38_chr22_test_FivePrimeUTRSeq',
                     gtf_path_or_url="./tests/data/gencode.v24.annotation_chr22_FivePrimeUTRSeq_testing.gtf",
                     transcript_fasta_paths_or_urls="./tests/data/hg38_chr22.fa",
                     )
    UTR_object = FivePrimeUTRSeq(data)
    return UTR_object


@pytest.fixture()
def FivePrimeUTRSeq_object_minus_strand():
    data_minus = GFGenome(reference_name='hg38_test_FivePrimeUTRSeq_minus',
                     annotation_name='hg38_chr22_test_FivePrimeUTRSeq_minus',
                     gtf_path_or_url="./tests/data/gencode.v24.annotation_chr22_FivePrimeUTRSeq_testing_minus-strand.gtf",
                     transcript_fasta_paths_or_urls="./tests/data/hg38_chr22.fa",
                     )
    UTR_object = FivePrimeUTRSeq(data_minus)
    return UTR_object

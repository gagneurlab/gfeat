import pytest
import pyensembl
from pysam import FastaFile
from pybedtools import BedTool
from gfeat.genome import GFGenome
from gfeat.transcript import GFTranscript


@pytest.fixture()
def genome():
    data = GFGenome(reference_name='hg38_test',
                     annotation_name='hg38_chr22_test',
                     gtf_path_or_url="./tests/data/gencode.v24.annotation_chr22.gtf",
                     transcript_fasta_paths_or_urls="./tests/data/hg38_chr22.fa",
                     )
    return data


@pytest.fixture()
def transcript():
    data = GFGenome(reference_name='hg38_test',
                     annotation_name='hg38_chr22_test',
                     gtf_path_or_url="./tests/data/gencode.v24.annotation_chr22.gtf",
                     transcript_fasta_paths_or_urls="./tests/data/hg38_chr22.fa",
                     )
    test_transcript = GFTranscript('ENST00000624155.1', 'BAGE5-201', '22', 11066501, 11068089, '-', 'None', 'ENSG00000279973.1', data)
    return test_transcript


@pytest.fixture()
def transcript2():
    data = GFGenome(reference_name='hg38_test',
                     annotation_name='hg38_chr22_test',
                     gtf_path_or_url="./tests/data/gencode.v24.annotation_chr22.gtf",
                     transcript_fasta_paths_or_urls="./tests/data/hg38_chr22.fa",
                     )
    test_transcript = GFTranscript('ENST00000343518.10', 'POTEH-001', '22', 15690026, 15721522, '+', 'None', 'ENSG00000198062.14', data)
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


@pytest.fixture()
def vcf_none():
    test_vcf = None
    return test_vcf


@pytest.fixture()
def ref_check():
    ref_check_test = True
    return ref_check_test


@pytest.fixture()
def seq():
    s = "CCCTGAGTCATCCTTGC"
    return s

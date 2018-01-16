import pytest
import pyensembl
from pyensembl import Genome
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
def interval():
    gtf = BedTool("/Users/veronikakotova/gfeat/tests/data/chr22_testing.gtf")
    test_interval = gtf[0]  # get the first element in the bed file, returns an Interval object
    return test_interval

@pytest.fixture()
def fasta():
    test_fasta = FastaFile("/Users/veronikakotova/gfeat/tests/data/chr22_testing.fa")
    return test_fasta

@pytest.fixture()
def vcf():
    test_vcf = "/Users/veronikakotova/gfeat/tests/data/49470G_chr22_testing.vcf.gz"
    return test_vcf

# data = Genome(reference_name='GRCh38',
#     annotation_name='my_genome_features',
#     gtf_path_or_url='//Users/veronikakotova/gfeat/tests/data/gencode.v24.annotation_chr22.gtf')
# # parse GTF and construct database of genomic features
# data.index()
# gene_names = data.gene_names_at_locus(contig=6, position=29945884)

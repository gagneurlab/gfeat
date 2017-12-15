import pytest
import pyensembl
from gfeat.genome import GFGenome
from gfeat.gene import GFGene

from gfeat.transcript import GFTranscript

@pytest.fixture()
def gene():
    data1 = GFGenome(reference_name='hg38_test',
                     annotation_name='hg38_chr22_test',
                     gtf_path_or_url='/Users/veronikakotova/Desktop/checking_pyensembl/gencode.v24.annotation_chr22_test.gtf',
                     transcript_fasta_paths_or_urls='/Users/veronikakotova/Desktop/checking_pyensembl/hg38_chr22_test.fa',
                     )
    test_gene = GFGene('ENSG00000008735.13', 'MAPK8IP2', '22', 50600685, 50613981, '+', None, data1)
    return test_gene

import pytest
import pyensembl
from pyensembl import Genome

from gfeat.transcript import GFTranscript

@pytest.fixture()
def transcript():
    data = pyensembl.ensembl_release.EnsemblRelease(75)
    test_transcript = GFTranscript("ENST00000369985", "MYO6-001", "6", 76458926, 76629253, "+", "protein_coding",
                                "ENSG00000196586", data)
    return test_transcript

# data = Genome(reference_name='GRCh38',
#     annotation_name='my_genome_features',
#     gtf_path_or_url='//Users/veronikakotova/gfeat/tests/data/gencode.v24.annotation_chr22.gtf')
# # parse GTF and construct database of genomic features
# data.index()
# gene_names = data.gene_names_at_locus(contig=6, position=29945884)

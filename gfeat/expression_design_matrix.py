from pyensembl import EnsemblRelease
from gfeat.genome import GFGenome
from cyvcf2 import VCF


def create_design_matrix(vcf, save_to=None, sample_id=None, ensembl_version=None, gtf=None, fasta=None):
    """

    :param vcf:
    :param save_to:
    :param sample_id:
    :param ensembl_version:
    :param gtf:
    :param fasta:
    :return:
    """

    if ensembl_version is not None:
        genome = EnsemblRelease(ensembl_version)
        GFdata = GFGenome(copy_genome=genome)
    else:
        GFdata = GFGenome(reference_name='NAME',
                          annotation_name='NAME',
                          gtf_path_or_url=gtf,
                          transcript_fasta_paths_or_urls=fasta,
                          )

    vcf_fields = ["G5", "GENEINFO", "GMAF", "GNO", "KGPilot123", "RSPOS", "SAO", "SLO", "SSR", "VC", "VLD", "WGT",
                  "dbSNPBuildID", "HD", "PH2", "G5A", "PM", "PMC", "ASP", "RV", "S3D", "GCF", "CLN", "LSD", "NOV",
                  "OM", "CFL", "NOV", "TPA", "MTP"]

    vcf_file = VCF(vcf)
    for field in vcf_fields:
        vcf_file.add_info_to_header({"ID": field, "Number": 1, "Type": "String", "Description": "dummy"})

    if save_to is not None:
        if save_to.endswith('.csv'):
            file = save_to
        else:
            file = save_to + "/" + sample_id + "_design_matrix.csv"
        df_design_matrix = GFdata.get_nucleobase_mutation_table(vcf)
        if sample_id is not None:
            df_design_matrix["sample_id"] = sample_id
        df_design_matrix.to_csv(file)
        pass
    else:
        return GFdata.get_nucleobase_mutation_table(vcf)

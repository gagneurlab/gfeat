from pyensembl import EnsemblRelease
from gfeat.genome import GFGenome
from cyvcf2 import VCF


def create_design_matrix(vcf, save_to=None, sample_id=None, ensembl_version=None, gtf=None, fasta=None):
    """
    1 possible model built using gfeat
    Creates and saves to a csv file a design matrix consisting of columns which correspond to Kozak sequence and stop
    codon context. If the first letter of the column names is "K", then the column corresponds to the Kozak
    sequence; "S" – to the stop codon context. Example: "K_3" means 3rd position in Kozak sequence.
    0 means that there is no variant, 1 – heterozygous variant, 2 – homozygous
    Please note that only either ensembl_version or gtf and fasta has to be provided. In case if all three variables
    are provided, only ensembl_version is going to be used
    :param vcf: string, path to the vcf.gz
    :param save_to: string, path to the folder or the file name where to save the pd.DataFrame output table
    :param sample_id: string, id of the sample
    :param ensembl_version: int, Ensembl version
           MIN_ENSEMBL_RELEASE = 54
           MAX_ENSEMBL_RELEASE = 85
    :param gtf: string, path to the gtf file; as an alternative to the Ensembl version, user gtf and Fasta files
                can be provided
    :param fasta: string, path to the Fasta file, string; as an alternative to the Ensembl version, user gtf
                  and Fasta files can be provided
    :return: pd.DataFrame with the following columns: K_i – where i shows position in Kozak sequence;
                                                      S_i – where i shows position in stop codon context,
                                                      gene_id, transctipt_id, sample_id
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

    # temporarily adding fields to the csv file in order to suppress the output and speed up the program
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

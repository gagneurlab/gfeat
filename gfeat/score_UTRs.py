from gfeat.upstreamAUG import UpstreamAUG
from gfeat.UTR import FivePrimeUTRSeq
from gfeat.units import mutate_sequence_Interval_vcf
from gfeat.genome import GFGenome
from gfeat.common_methods import reverse_complement
from pyensembl import EnsemblRelease
from cyvcf2 import VCF
import pandas as pd
import numpy as np


def score_utrs(vcf, save_to_csv, path, ensembl_version=None, gtf=None, fasta=None):
    """
    :param vcf: string, path to the vcf.gz
    :param save_to_csv: bool, whether to save the pd.DataFrame output table to a csv file or not
    :param path: string, path to the folder where to save the pd.DataFrame output table
    :param ensembl_versio:, int, Ensembl version
    :param grf: string, path to the gtf file, as an alternative to the Ensembl version, user gtf and Fasta files
                can be provided
    :param fasta: string, path to the Fasta file, string, as an alternative to the Ensembl version, user gtf
           and Fasta files can be provided
    :return pd.DataFrame with 6 columns: Gene, Transcript, in-frame_no_uORF, in-frame_uORF, not_in-frame_no_uORF,
            not_in-frame_uORF
    """
    model = UpstreamAUG(True, True)  # True: in frame, True: ORF

    if ensembl_version != None:
        data = EnsemblRelease(ensembl_version)
    else:
        data = GFGenome(reference_name='NAME',
                        annotation_name='NAME',
                        gtf_path_or_url=gtf,
                        transcript_fasta_paths_or_urls=fasta,
                        )

    contigs = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
               '20', '21', '22', 'X', 'Y']

    vcf_fields = ["G5", "GENEINFO", "GMAF", "GNO", "KGPilot123", "RSPOS", "SAO", "SLO", "SSR", "VC", "VLD", "WGT",
                  "dbSNPBuildID", "HD", "PH2", "G5A", "PM", "PMC", "ASP", "RV", "S3D", "GCF", "CLN", "LSD", "NOV",
                  "OM", "CFL", "NOV", "TPA", "MTP"]

    vcf_file = VCF(vcf)
    for field in vcf_fields:
        vcf_file.add_info_to_header({"ID": field, "Number": 1, "Type": "String", "Description": "dummy"})

    if save_to_csv:
        dictionary = {'Gene': [], 'Transcript': [], 'not_in-frame_no_uORF': [], 'not_in-frame_uORF': [],
                      'in-frame_no_uORF': [], 'in-frame_uORF': []}
        error_dictionary = {'More_than_14_mutations': []}
        sum_dictionary = {'Contig': [], 'Total_transcript_num': [], 'Mutated_transcript_num': [],
                          'not_in-frame_no_uORF_num': [],
                          'not_in-frame_uORF': [], 'in-frame_no_uORF_num': [], 'in-frame_uORF_num': []}
        mutated_dictionary = {'Gene': [], 'Transcript': [], 'Type': [], 'Position': []}

        df0 = pd.DataFrame(data=dictionary)
        df1 = pd.DataFrame(data=error_dictionary)
        df2 = pd.DataFrame(data=sum_dictionary)

        df0.to_csv(path + "/AUGs.csv")
        df1.to_csv(path + "/too_many_mutated.csv")
        df2.to_csv(path + "/summary.csv")

    for contig in contigs:

        ds = FivePrimeUTRSeq(data, False, contig, '+')  # Todo: check

        dictionary = {'Gene': [], 'Transcript': [], 'not_in-frame_no_uORF': [], 'not_in-frame_uORF': [],
                      'in-frame_no_uORF': [], 'in-frame_uORF': []}
        sum_dictionary = {'Contig': [contig], 'Total_transcript_num': [len(ds) - 2], 'Mutated_transcript_num': [0],
                          'not_in-frame_no_uORF_num': [0],
                          'not_in-frame_uORF_num': [0], 'in-frame_no_uORF_num': [0], 'in-frame_uORF_num': [0]}

        list_many_mut = []

        for i in range(2, len(ds)):  # ds[0] reserved for " Pyensembl error"
            sample = ds[i]
            ref_seq = ""
            if (list(sample.keys())[0] != "transcripts") and (list(sample.keys())[0] != "exons"):
                ref_seq = list(sample.keys())[0]

            list_exon_mut_seq = []

            mutated = 0

            for exon in sample["exons"]:
                e, m = mutate_sequence_Interval_vcf(exon[1], exon[0], vcf_file)
                mutated = mutated + m
                list_exon_mut_seq.append(e)

            if 0 < mutated < 14:

                sum_dictionary['Mutated_transcript_num'][0] = sum_dictionary['Mutated_transcript_num'][0] + mutated

                count_of_combinations = 1
                list_mut_seq = []
                for x in range(0, len(list_exon_mut_seq)):
                    count_of_combinations = count_of_combinations * len(list_exon_mut_seq[x])
                indexes = [0] * len(list_exon_mut_seq)
                counter = 1
                sequence = ""
                while counter <= count_of_combinations:
                    for x in range(0, len(list_exon_mut_seq)):
                        if indexes[x] == len(list_exon_mut_seq[x]):
                            indexes[x] = 0
                            indexes[x + 1] = indexes[x + 1] + 1
                    for x in range(0, len(list_exon_mut_seq)):
                        sequence = sequence + (list_exon_mut_seq[x])[indexes[x]][0]
                    list_mut_seq.append(sequence)
                    sequence = ""
                    indexes[0] = indexes[0] + 1
                    counter = counter + 1

                dictionary["Transcript"].append(sample["transcripts"])
                genes = []
                temp = ''
                for transcript in sample["transcripts"]:
                    gene = data.transcript_by_id(transcript).gene.name
                    if gene != temp:
                        genes.append(gene)
                        temp = gene
                dictionary["Gene"].append(genes)

                if sample[ref_seq].strand == "-":
                    model.predict_on_sample_with_pos_pandas(reverse_complement(ref_seq), dictionary, '-',
                                                            sample[ref_seq].start)
                    cur_0_0 = dictionary['not_in-frame_no_uORF'][len(dictionary['not_in-frame_no_uORF']) - 1]
                    cur_0_1 = dictionary['not_in-frame_uORF'][len(dictionary['not_in-frame_uORF']) - 1]
                    cur_1_0 = dictionary['in-frame_no_uORF'][len(dictionary['in-frame_no_uORF']) - 1]
                    cur_1_1 = dictionary['in-frame_uORF'][len(dictionary['in-frame_uORF']) - 1]
                else:
                    model.predict_on_sample_with_pos_pandas(ref_seq, dictionary,'+', sample[ref_seq].start)
                    cur_0_0 = dictionary['not_in-frame_no_uORF'][len(dictionary['not_in-frame_no_uORF']) - 1]
                    cur_0_1 = dictionary['not_in-frame_uORF'][len(dictionary['not_in-frame_uORF']) - 1]
                    cur_1_0 = dictionary['in-frame_no_uORF'][len(dictionary['in-frame_no_uORF']) - 1]
                    cur_1_1 = dictionary['in-frame_uORF'][len(dictionary['in-frame_uORF']) - 1]

                for seq in list_mut_seq:
                    dictionary["Transcript"].append(sample["transcripts"])
                    dictionary["Gene"].append(genes)
                    if sample[ref_seq].strand == "-":
                        model.predict_on_sample_with_pos_pandas(seq, dictionary, '-', sample[ref_seq].start)
                    else:
                        model.predict_on_sample_with_pos_pandas(seq, dictionary, '+', sample[ref_seq].start)
                    mut_0_0 = dictionary['not_in-frame_no_uORF'][len(dictionary['not_in-frame_no_uORF']) - 1]
                    if not np.array_equal(mut_0_0, cur_0_0):
                        sum_dictionary['not_in-frame_no_uORF_num'][0] = sum_dictionary['not_in-frame_no_uORF_num'][
                                                                            0] + 1
                        mutated_dictionary['Gene'].append(genes)
                        mutated_dictionary['Transcript'].append(sample["transcripts"])
                        mutated_dictionary['Type'].append("not_in-frame_no_uORF")
                        mutated_dictionary['Position'].append(
                            np.setdiff1d(np.union1d(mut_0_0, cur_0_0), np.intersect1d(mut_0_0, cur_0_0)))
                    mut_0_1 = dictionary['not_in-frame_uORF'][len(dictionary['not_in-frame_uORF']) - 1]
                    if not np.array_equal(mut_0_1, cur_0_1):
                        sum_dictionary['not_in-frame_uORF_num'][0] = sum_dictionary['not_in-frame_uORF_num'][0] + 1
                        mutated_dictionary['Gene'].append(genes)
                        mutated_dictionary['Transcript'].append(sample["transcripts"])
                        mutated_dictionary['Type'].append("not_in-frame_uORF")
                        mutated_dictionary['Position'].append(
                            np.setdiff1d(np.union1d(mut_0_1, cur_0_1), np.intersect1d(mut_0_1, cur_0_1)))
                    mut_1_0 = dictionary['in-frame_no_uORF'][len(dictionary['in-frame_no_uORF']) - 1]
                    if not np.array_equal(mut_1_0, cur_1_0):
                        sum_dictionary['in-frame_no_uORF_num'][0] = sum_dictionary['in-frame_no_uORF_num'][0] + 1
                        mutated_dictionary['Gene'].append(genes)
                        mutated_dictionary['Transcript'].append(sample["transcripts"])
                        mutated_dictionary['Type'].append("in-frame_no_uORF")
                        mutated_dictionary['Position'].append(
                            np.setdiff1d(np.union1d(mut_1_0, cur_1_0), np.intersect1d(mut_1_0, cur_1_0)))
                    mut_1_1 = dictionary['in-frame_uORF'][len(dictionary['in-frame_uORF']) - 1]
                    if not np.array_equal(mut_1_1, cur_1_1):
                        sum_dictionary['in-frame_uORF_num'][0] = sum_dictionary['in-frame_uORF_num'][0] + 1
                        mutated_dictionary['Gene'].append(genes)
                        mutated_dictionary['Transcript'].append(sample["transcripts"])
                        mutated_dictionary['Type'].append("in-frame_uORF")
                        mutated_dictionary['Position'].append(
                            np.setdiff1d(np.union1d(mut_1_1, cur_1_1), np.intersect1d(mut_1_1, cur_1_1)))
            else:
                if mutated > 14:
                    list_many_mut.append((ds[i])["transcripts"][0])

        if save_to_csv:
            df0 = pd.DataFrame(data=dictionary)
            error_list = list_many_mut
            error_dictionary = {"More_than_14_mutations": error_list}
            df1 = pd.DataFrame(data=error_dictionary)

            with open(path + "/AUGs.csv", 'a') as f:
                df0.to_csv(f, header=False)

            with open(path + "/too_many_mutated.csv", 'a') as f:
                df1.to_csv(f, header=False)

        ds = FivePrimeUTRSeq(data, False, contig, '-')

        dictionary = {'Gene': [], 'Transcript': [], 'not_in-frame_no_uORF': [], 'not_in-frame_uORF': [],
                      'in-frame_no_uORF': [], 'in-frame_uORF': []}

        sum_dictionary['Total_transcript_num'][0] = sum_dictionary['Total_transcript_num'][0] + len(ds) - 2

        list_many_mut = []

        for i in range(2, len(ds)):  # ds[1]
            sample = ds[i]
            ref_seq = ""
            if (list(sample.keys())[0] != "transcripts") and (list(sample.keys())[0] != "exons"):
                ref_seq = list(sample.keys())[0]

            list_exon_mut_seq = []

            mutated = 0

            for exon in sample["exons"]:
                e, m = mutate_sequence_Interval_vcf(exon[1], exon[0], vcf_file)
                mutated = mutated + m
                list_exon_mut_seq.append(e)

            if 0 < mutated < 14:

                sum_dictionary['Mutated_transcript_num'][0] = sum_dictionary['Mutated_transcript_num'][0] + mutated

                count_of_combinations = 1
                list_mut_seq = []
                for x in range(0, len(list_exon_mut_seq)):
                    count_of_combinations = count_of_combinations * len(list_exon_mut_seq[x])
                indexes = [0] * len(list_exon_mut_seq)
                counter = 1
                sequence = ""
                while counter <= count_of_combinations:
                    for x in range(0, len(list_exon_mut_seq)):
                        if indexes[x] == len(list_exon_mut_seq[x]):
                            indexes[x] = 0
                            indexes[x + 1] = indexes[x + 1] + 1
                    for x in range(0, len(list_exon_mut_seq)):
                        sequence = sequence + (list_exon_mut_seq[x])[indexes[x]][0]
                    list_mut_seq.append(sequence)
                    sequence = ""
                    indexes[0] = indexes[0] + 1
                    counter = counter + 1
                dictionary["Transcript"].append(sample["transcripts"])
                genes = []
                temp = ''
                for transcript in sample["transcripts"]:
                    gene = data.transcript_by_id(transcript).gene.name
                    if gene != temp:
                        genes.append(gene)
                        temp = gene
                dictionary["Gene"].append(genes)

                if sample[ref_seq].strand == "-":
                    model.predict_on_sample_with_pos_pandas(reverse_complement(ref_seq), dictionary, '-',
                                                            sample[ref_seq].start)
                    cur_0_0 = dictionary['not_in-frame_no_uORF'][len(dictionary['not_in-frame_no_uORF']) - 1]
                    cur_0_1 = dictionary['not_in-frame_uORF'][len(dictionary['not_in-frame_uORF']) - 1]
                    cur_1_0 = dictionary['in-frame_no_uORF'][len(dictionary['in-frame_no_uORF']) - 1]
                    cur_1_1 = dictionary['in-frame_uORF'][len(dictionary['in-frame_uORF']) - 1]
                else:
                    model.predict_on_sample_with_pos_pandas(ref_seq, dictionary, '+', sample[ref_seq].start)
                    cur_0_0 = dictionary['not_in-frame_no_uORF'][len(dictionary['not_in-frame_no_uORF']) - 1]
                    cur_0_1 = dictionary['not_in-frame_uORF'][len(dictionary['not_in-frame_uORF']) - 1]
                    cur_1_0 = dictionary['in-frame_no_uORF'][len(dictionary['in-frame_no_uORF']) - 1]
                    cur_1_1 = dictionary['in-frame_uORF'][len(dictionary['in-frame_uORF']) - 1]

                for seq in list_mut_seq:
                    dictionary["Transcript"].append(sample["transcripts"])
                    dictionary["Gene"].append(genes)
                    if sample[ref_seq].strand == "-":
                        model.predict_on_sample_with_pos_pandas(seq, dictionary, '-', sample[ref_seq].start)
                    else:
                        model.predict_on_sample_with_pos_pandas(seq, dictionary, '+', sample[ref_seq].start)
                    mut_0_0 = dictionary['not_in-frame_no_uORF'][len(dictionary['not_in-frame_no_uORF']) - 1]
                    if not np.array_equal(mut_0_0, cur_0_0):
                        sum_dictionary['not_in-frame_no_uORF_num'][0] = sum_dictionary['not_in-frame_no_uORF_num'][
                                                                            0] + 1
                        mutated_dictionary['Gene'].append(genes)
                        mutated_dictionary['Transcript'].append(sample["transcripts"])
                        mutated_dictionary['Type'].append("not_in-frame_no_uORF")
                        mutated_dictionary['Position'].append(
                            np.setdiff1d(np.union1d(mut_0_0, cur_0_0), np.intersect1d(mut_0_0, cur_0_0)))
                    mut_0_1 = dictionary['not_in-frame_uORF'][len(dictionary['not_in-frame_uORF']) - 1]
                    if not np.array_equal(mut_0_1, cur_0_1):
                        sum_dictionary['not_in-frame_uORF_num'][0] = sum_dictionary['not_in-frame_uORF_num'][0] + 1
                        mutated_dictionary['Gene'].append(genes)
                        mutated_dictionary['Transcript'].append(sample["transcripts"])
                        mutated_dictionary['Type'].append("not_in-frame_uORF")
                        mutated_dictionary['Position'].append(
                            np.setdiff1d(np.union1d(mut_0_1, cur_0_1), np.intersect1d(mut_0_1, cur_0_1)))
                    mut_1_0 = dictionary['in-frame_no_uORF'][len(dictionary['in-frame_no_uORF']) - 1]
                    if not np.array_equal(mut_1_0, cur_1_0):
                        sum_dictionary['in-frame_no_uORF_num'][0] = sum_dictionary['in-frame_no_uORF_num'][0] + 1
                        mutated_dictionary['Gene'].append(genes)
                        mutated_dictionary['Transcript'].append(sample["transcripts"])
                        mutated_dictionary['Type'].append("in-frame_no_uORF")
                        mutated_dictionary['Position'].append(
                            np.setdiff1d(np.union1d(mut_1_0, cur_1_0), np.intersect1d(mut_1_0, cur_1_0)))
                    mut_1_1 = dictionary['in-frame_uORF'][len(dictionary['in-frame_uORF']) - 1]
                    if not np.array_equal(mut_1_1, cur_1_1):
                        sum_dictionary['in-frame_uORF_num'][0] = sum_dictionary['in-frame_uORF_num'][0] + 1
                        mutated_dictionary['Gene'].append(genes)
                        mutated_dictionary['Transcript'].append(sample["transcripts"])
                        mutated_dictionary['Type'].append("in-frame_uORF")
                        mutated_dictionary['Position'].append(
                            np.setdiff1d(np.union1d(mut_1_1, cur_1_1), np.intersect1d(mut_1_1, cur_1_1)))
            else:
                if mutated > 14:
                    list_many_mut.append((ds[i])["transcripts"][0])

        if save_to_csv:
            df0 = pd.DataFrame(data=dictionary)
            error_list = list_many_mut
            error_dictionary = {"More_than_14_mutations": error_list}
            df2 = pd.DataFrame(data=sum_dictionary)
            df1 = pd.DataFrame(data=error_dictionary)

            with open(path + "/AUGs.csv", 'a') as f:
                df0.to_csv(f, header=False)

            with open(path + "/too_many_mutated.csv", 'a') as f:
                df1.to_csv(f, header=False)

            with open(path + "/summary.csv", 'a') as f:
                df2.to_csv(f, header=False)

    df3 = pd.DataFrame(data=mutated_dictionary)
    df3.to_csv(path + "/mutated.csv")

    return df3

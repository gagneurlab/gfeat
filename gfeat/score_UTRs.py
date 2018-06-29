from gfeat.upstreamAUG import UpstreamAUG
from gfeat.UTR import FivePrimeUTRSeq
from gfeat.units import VCFMutator
from gfeat.genome import GFGenome
from gfeat.utils import reverse_complement
from pyensembl import EnsemblRelease
from cyvcf2 import VCF
import pandas as pd
import numpy as np


def score_utrs(vcf, save_extra_to_csv, save_to, patient_id, ensembl_version=None, gtf=None, fasta=None):
    """
    1 possible model built using gfeat
    Looks for creation and deletion of AUGs which are located upstream from the corresponding canonical start.
    Writes down to the file "mutated.csv" geneId, position of an AUG loss/creation, transcriptID,
    AUGType (in-frame_no_uORF, in-frame_uORF, not_in-frame_no_uORF, not_in-frame_uORF),
    alternative value of the nucleobase, chromosome, reference value of the nucleobase, sampleID
    :param vcf: string, path to the vcf.gz
    :param save_extra_to_csv: bool, whether to save intermediate information (all AUGs, transcripts with too many
                              mutations and some very basic statistics) to a csv file or not
           currently work only if save_to is the name of the folder
    :param save_to: string, path to the folder or the file name where to save the pd.DataFrame output table
    :param ensembl_version:, int, Ensembl version
           MIN_ENSEMBL_RELEASE = 54
           MAX_ENSEMBL_RELEASE = 85
    :param gtf: string, path to the gtf file, as an alternative to the Ensembl version, user gtf and Fasta files
                can be provided
    :param fasta: string, path to the Fasta file, string, as an alternative to the Ensembl version, user gtf
                  and Fasta files can be provided
    :return pd.DataFrame with 6 columns: geneId, position, transcriptID, AUGType, alt, chr, ref, sampleID
    """
    model = UpstreamAUG(True, True)

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


    dictionary = {'Gene': [], 'Transcript': [], 'not_in-frame_no_uORF': [], 'not_in-frame_uORF': [],
                      'in-frame_no_uORF': [], 'in-frame_uORF': []}
    error_dictionary = {'More_than_14_mutations': []}
    sum_dictionary = {'Contig': [], 'Total_transcript_num': [], 'Mutated_transcript_num': [],
                          'not_in-frame_no_uORF_num': [],
                          'not_in-frame_uORF': [], 'in-frame_no_uORF_num': [], 'in-frame_uORF_num': []}
    mutated_dictionary = {'Gene': [], 'Transcript': [], 'Type': [], 'Position': [], 'ref': [],
                              'alt': [], 'chr': [], 'sample': []}
    if save_extra_to_csv:

        pass

        # df0 = pd.DataFrame(data=dictionary)
        # df1 = pd.DataFrame(data=error_dictionary)
        # df2 = pd.DataFrame(data=sum_dictionary)
        #
        # df0.to_csv(save_to + "/" + patient_id + "_AUGs.csv")
        # df1.to_csv(save_to + "/" + patient_id + "_too_many_mutated.csv")
        # df2.to_csv(save_to + "/" + patient_id + "_summary.csv")

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
            elif (list(sample.keys())[1] != "transcripts") and (list(sample.keys())[1] != "exons"):
                ref_seq = list(sample.keys())[1]
            else:
                ref_seq = list(sample.keys())[2]

            list_exon_mut_seq = []

            mutated = 0

            mutator = VCFMutator(False, True, vcf_file, True)

            for exon in sample["exons"]:
                e, m = mutator.mutate_sequence(exon[1], False, exon[0])
                mutated = m < 14
                list_exon_mut_seq.append(e)

            if mutated:

                sum_dictionary['Mutated_transcript_num'][0] = sum_dictionary['Mutated_transcript_num'][0] + mutated

                count_of_combinations = 1
                list_of_mut_seq_inf = []
                for x in range(0, len(list_exon_mut_seq)):
                    count_of_combinations = count_of_combinations * len(list_exon_mut_seq[x])
                indexes = [0] * len(list_exon_mut_seq)
                counter = 1
                sequence = ""
                tupl_ref = ()
                tupl_alt = ()
                tupl_pos = ()
                while counter <= count_of_combinations:
                    for x in range(0, len(list_exon_mut_seq)):
                        if indexes[x] == len(list_exon_mut_seq[x]):
                            indexes[x] = 0
                            indexes[x + 1] = indexes[x + 1] + 1
                    for x in range(0, len(list_exon_mut_seq)):
                        sequence = sequence + (list_exon_mut_seq[x])[indexes[x]][0]
                        tupl_pos = tupl_pos + (list_exon_mut_seq[x])[indexes[x]][1]
                        tupl_ref = tupl_ref + (list_exon_mut_seq[x])[indexes[x]][2]
                        tupl_alt = tupl_alt + (list_exon_mut_seq[x])[indexes[x]][3]
                    list_of_mut_seq_inf.append({'seq': sequence, 'pos': tupl_pos, 'ref': tupl_ref, 'alt': tupl_alt})
                    sequence = ""
                    tupl_ref = ()
                    tupl_alt = ()
                    tupl_pos = ()
                    indexes[0] = indexes[0] + 1
                    counter = counter + 1

                dictionary["Transcript"].append(sample["transcripts"])
                genes = []
                temp = ''
                for transcript in sample["transcripts"]:
                    gene = data.transcript_by_id(transcript).gene.id
                    if gene != temp:
                        genes.append(gene)
                        temp = gene
                dictionary["Gene"].append(genes)

                if sample[ref_seq].strand == "-":
                    model.predict_on_sample_with_stop_pandas(reverse_complement(ref_seq), dictionary, '-',
                                                             sample[ref_seq].start)
                    cur_0_0 = dictionary['not_in-frame_no_uORF'][len(dictionary['not_in-frame_no_uORF']) - 1]
                    cur_0_1 = dictionary['not_in-frame_uORF'][len(dictionary['not_in-frame_uORF']) - 1]
                    cur_1_0 = dictionary['in-frame_no_uORF'][len(dictionary['in-frame_no_uORF']) - 1]
                    cur_1_1 = dictionary['in-frame_uORF'][len(dictionary['in-frame_uORF']) - 1]
                else:
                    model.predict_on_sample_with_stop_pandas(ref_seq, dictionary, '+', sample[ref_seq].start)
                    cur_0_0 = dictionary['not_in-frame_no_uORF'][len(dictionary['not_in-frame_no_uORF']) - 1]
                    cur_0_1 = dictionary['not_in-frame_uORF'][len(dictionary['not_in-frame_uORF']) - 1]
                    cur_1_0 = dictionary['in-frame_no_uORF'][len(dictionary['in-frame_no_uORF']) - 1]
                    cur_1_1 = dictionary['in-frame_uORF'][len(dictionary['in-frame_uORF']) - 1]

                for info in list_of_mut_seq_inf:
                    dictionary["Transcript"].append(sample["transcripts"])
                    dictionary["Gene"].append(genes)
                    if sample[ref_seq].strand == "-":
                        model.predict_on_sample_with_stop_pandas(info['seq'], dictionary, '-', sample[ref_seq].start)
                    else:
                        model.predict_on_sample_with_stop_pandas(info['seq'], dictionary, '+', sample[ref_seq].start)
                    mut_0_0 = dictionary['not_in-frame_no_uORF'][len(dictionary['not_in-frame_no_uORF']) - 1]
                    if not np.array_equal(mut_0_0, cur_0_0):
                        for pos in np.setdiff1d(np.union1d(mut_0_0, cur_0_0), np.intersect1d(mut_0_0, cur_0_0)):
                            for i in range(0, len(info['pos'])):
                                if (info['pos'][i] >= pos) and (info['pos'][i] < pos + 3):
                                    sum_dictionary['not_in-frame_no_uORF_num'][0] = \
                                    sum_dictionary['not_in-frame_no_uORF_num'][
                                        0] + 1
                                    for transcript_id in sample["transcripts"]:
                                        mutated_dictionary['Gene'].append(genes)
                                        mutated_dictionary['Transcript'].append(transcript_id)
                                        mutated_dictionary['Type'].append("not_in-frame_no_uORF")
                                        mutated_dictionary['Position'].append(info['pos'][i])
                                        mutated_dictionary['ref'].append(info['ref'][i])
                                        mutated_dictionary['alt'].append(info['alt'][i])
                                        mutated_dictionary['chr'].append(contig)
                                        mutated_dictionary['sample'].append(patient_id)
                    mut_0_1 = dictionary['not_in-frame_uORF'][len(dictionary['not_in-frame_uORF']) - 1]
                    if not np.array_equal(mut_0_1, cur_0_1):
                        for pos in np.setdiff1d(np.union1d(mut_0_1, cur_0_1), np.intersect1d(mut_0_1, cur_0_1)):
                            for i in range(0, len(info['pos'])):
                                if (info['pos'][i] >= pos) and (info['pos'][i] < pos + 3):
                                    sum_dictionary['not_in-frame_uORF_num'][0] = \
                                    sum_dictionary['not_in-frame_uORF_num'][
                                        0] + 1
                                    for transcript_id in sample["transcripts"]:
                                        mutated_dictionary['Gene'].append(genes)
                                        mutated_dictionary['Transcript'].append(transcript_id)
                                        mutated_dictionary['Type'].append("not_in-frame_uORF")
                                        mutated_dictionary['Position'].append(info['pos'][i])
                                        mutated_dictionary['ref'].append(info['ref'][i])
                                        mutated_dictionary['alt'].append(info['alt'][i])
                                        mutated_dictionary['chr'].append(contig)
                                        mutated_dictionary['sample'].append(patient_id)
                    mut_1_0 = dictionary['in-frame_no_uORF'][len(dictionary['in-frame_no_uORF']) - 1]
                    if not np.array_equal(mut_1_0, cur_1_0):
                        for pos in np.setdiff1d(np.union1d(mut_1_0, cur_1_0), np.intersect1d(mut_1_0, cur_1_0)):
                            for i in range(0, len(info['pos'])):
                                if (info['pos'][i] >= pos) and (info['pos'][i] < pos + 3):
                                    sum_dictionary['in-frame_no_uORF_num'][0] = sum_dictionary['in-frame_no_uORF_num'][
                                                                                    0] + 1
                                    for transcript_id in sample["transcripts"]:
                                        mutated_dictionary['Gene'].append(genes)
                                        mutated_dictionary['Transcript'].append(transcript_id)
                                        mutated_dictionary['Type'].append("in-frame_no_uORF")
                                        mutated_dictionary['Position'].append(info['pos'][i])
                                        mutated_dictionary['ref'].append(info['ref'][i])
                                        mutated_dictionary['alt'].append(info['alt'][i])
                                        mutated_dictionary['chr'].append(contig)
                                        mutated_dictionary['sample'].append(patient_id)
                    mut_1_1 = dictionary['in-frame_uORF'][len(dictionary['in-frame_uORF']) - 1]
                    if not np.array_equal(mut_1_1, cur_1_1):
                        for pos in np.setdiff1d(np.union1d(mut_1_1, cur_1_1), np.intersect1d(mut_1_1, cur_1_1)):
                            for i in range(0, len(info['pos'])):
                                if (info['pos'][i] >= pos) and (info['pos'][i] < pos + 3):
                                    sum_dictionary['in-frame_uORF_num'][0] = sum_dictionary['in-frame_uORF_num'][
                                                                                 0] + 1
                                    for transcript_id in sample["transcripts"]:
                                        mutated_dictionary['Gene'].append(genes)
                                        mutated_dictionary['Transcript'].append(transcript_id)
                                        mutated_dictionary['Type'].append("in-frame_uORF")
                                        mutated_dictionary['Position'].append(info['pos'][i])
                                        mutated_dictionary['ref'].append(info['ref'][i])
                                        mutated_dictionary['alt'].append(info['alt'][i])
                                        mutated_dictionary['chr'].append(contig)
                                        mutated_dictionary['sample'].append(patient_id)
            else:
                list_many_mut.append((ds[i])["transcripts"][0])

        if save_extra_to_csv:

            pass

            # df0 = pd.DataFrame(data=dictionary)
            # error_list = list_many_mut
            # error_dictionary = {"More_than_14_mutations": error_list}
            # df1 = pd.DataFrame(data=error_dictionary)
            #
            # with open(save_to + "/" + patient_id + "_AUGs.csv", 'a') as f:
            #     df0.to_csv(f, header=False)
            #
            # with open(save_to + "/" + patient_id + "_too_many_mutated.csv", 'a') as f:
            #     df1.to_csv(f, header=False)

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
            elif (list(sample.keys())[1] != "transcripts") and (list(sample.keys())[1] != "exons"):
                ref_seq = list(sample.keys())[1]
            else:
                ref_seq = list(sample.keys())[2]

            list_exon_mut_seq = []

            mutated = 0

            mutator = VCFMutator(False, True, vcf_file, True)

            for exon in sample["exons"]:
                e, m = mutator.mutate_sequence(exon[1], False, exon[0])
                mutated = m < 14
                list_exon_mut_seq.append(e)

            if mutated:

                sum_dictionary['Mutated_transcript_num'][0] = sum_dictionary['Mutated_transcript_num'][0] + mutated

                count_of_combinations = 1
                list_of_mut_seq_inf = []
                for x in range(0, len(list_exon_mut_seq)):
                    count_of_combinations = count_of_combinations * len(list_exon_mut_seq[x])
                indexes = [0] * len(list_exon_mut_seq)
                counter = 1
                sequence = ""
                tupl_ref = ()
                tupl_alt = ()
                tupl_pos = ()
                while counter <= count_of_combinations:
                    for x in range(0, len(list_exon_mut_seq)):
                        if indexes[x] == len(list_exon_mut_seq[x]):
                            indexes[x] = 0
                            indexes[x + 1] = indexes[x + 1] + 1
                    for x in range(0, len(list_exon_mut_seq)):
                        sequence = sequence + (list_exon_mut_seq[x])[indexes[x]][0]
                        tupl_pos = tupl_pos + (list_exon_mut_seq[x])[indexes[x]][1]
                        tupl_ref = tupl_ref + (list_exon_mut_seq[x])[indexes[x]][2]
                        tupl_alt = tupl_alt + (list_exon_mut_seq[x])[indexes[x]][3]
                    list_of_mut_seq_inf.append({'seq': sequence, 'pos': tupl_pos, 'ref': tupl_ref, 'alt': tupl_alt})
                    sequence = ""
                    tupl_ref = ()
                    tupl_alt = ()
                    tupl_pos = ()
                    indexes[0] = indexes[0] + 1
                    counter = counter + 1

                dictionary["Transcript"].append(sample["transcripts"])
                genes = []
                temp = ''
                for transcript in sample["transcripts"]:
                    gene = data.transcript_by_id(transcript).gene.id
                    if gene != temp:
                        genes.append(gene)
                        temp = gene
                dictionary["Gene"].append(genes)

                if sample[ref_seq].strand == "-":
                    model.predict_on_sample_with_stop_pandas(reverse_complement(ref_seq), dictionary, '-',
                                                             sample[ref_seq].start)
                    cur_0_0 = dictionary['not_in-frame_no_uORF'][len(dictionary['not_in-frame_no_uORF']) - 1]
                    cur_0_1 = dictionary['not_in-frame_uORF'][len(dictionary['not_in-frame_uORF']) - 1]
                    cur_1_0 = dictionary['in-frame_no_uORF'][len(dictionary['in-frame_no_uORF']) - 1]
                    cur_1_1 = dictionary['in-frame_uORF'][len(dictionary['in-frame_uORF']) - 1]
                else:
                    model.predict_on_sample_with_stop_pandas(ref_seq, dictionary, '+', sample[ref_seq].start)
                    cur_0_0 = dictionary['not_in-frame_no_uORF'][len(dictionary['not_in-frame_no_uORF']) - 1]
                    cur_0_1 = dictionary['not_in-frame_uORF'][len(dictionary['not_in-frame_uORF']) - 1]
                    cur_1_0 = dictionary['in-frame_no_uORF'][len(dictionary['in-frame_no_uORF']) - 1]
                    cur_1_1 = dictionary['in-frame_uORF'][len(dictionary['in-frame_uORF']) - 1]

                for info in list_of_mut_seq_inf:
                    dictionary["Transcript"].append(sample["transcripts"])
                    dictionary["Gene"].append(genes)
                    if sample[ref_seq].strand == "-":
                        model.predict_on_sample_with_stop_pandas(info['seq'], dictionary, '-', sample[ref_seq].start)
                    else:
                        model.predict_on_sample_with_stop_pandas(info['seq'], dictionary, '+', sample[ref_seq].start)
                    mut_0_0 = dictionary['not_in-frame_no_uORF'][len(dictionary['not_in-frame_no_uORF']) - 1]
                    if not np.array_equal(mut_0_0, cur_0_0):
                        for pos in np.setdiff1d(np.union1d(mut_0_0, cur_0_0), np.intersect1d(mut_0_0, cur_0_0)):
                            for i in range(0, len(info['pos'])):
                                if (info['pos'][i] <= pos) and (info['pos'][i] > pos - 3):
                                    sum_dictionary['not_in-frame_no_uORF_num'][0] = \
                                    sum_dictionary['not_in-frame_no_uORF_num'][
                                        0] + 1
                                    for transcript_id in sample["transcripts"]:
                                        mutated_dictionary['Gene'].append(genes)
                                        mutated_dictionary['Transcript'].append(transcript_id)
                                        mutated_dictionary['Type'].append("not_in-frame_no_uORF")
                                        mutated_dictionary['Position'].append(info['pos'][i])
                                        mutated_dictionary['ref'].append(info['ref'][i])
                                        mutated_dictionary['alt'].append(info['alt'][i])
                                        mutated_dictionary['chr'].append(contig)
                                        mutated_dictionary['sample'].append(patient_id)
                    mut_0_1 = dictionary['not_in-frame_uORF'][len(dictionary['not_in-frame_uORF']) - 1]
                    if not np.array_equal(mut_0_1, cur_0_1):
                        for pos in np.setdiff1d(np.union1d(mut_0_1, cur_0_1), np.intersect1d(mut_0_1, cur_0_1)):
                            for i in range(0, len(info['pos'])):
                                if (info['pos'][i] <= pos) and (info['pos'][i] > pos - 3):
                                    sum_dictionary['not_in-frame_uORF_num'][0] = \
                                    sum_dictionary['not_in-frame_uORF_num'][
                                        0] + 1
                                    for transcript_id in sample["transcripts"]:
                                        mutated_dictionary['Gene'].append(genes)
                                        mutated_dictionary['Transcript'].append(transcript_id)
                                        mutated_dictionary['Type'].append("not_in-frame_uORF")
                                        mutated_dictionary['Position'].append(info['pos'][i])
                                        mutated_dictionary['ref'].append(info['ref'][i])
                                        mutated_dictionary['alt'].append(info['alt'][i])
                                        mutated_dictionary['chr'].append(contig)
                                        mutated_dictionary['sample'].append(patient_id)
                    mut_1_0 = dictionary['in-frame_no_uORF'][len(dictionary['in-frame_no_uORF']) - 1]
                    if not np.array_equal(mut_1_0, cur_1_0):
                        for pos in np.setdiff1d(np.union1d(mut_1_0, cur_1_0), np.intersect1d(mut_1_0, cur_1_0)):
                            for i in range(0, len(info['pos'])):
                                if (info['pos'][i] <= pos) and (info['pos'][i] > pos - 3):
                                    sum_dictionary['in-frame_no_uORF_num'][0] = sum_dictionary['in-frame_no_uORF_num'][
                                                                                    0] + 1
                                    for transcript_id in sample["transcripts"]:
                                        mutated_dictionary['Gene'].append(genes)
                                        mutated_dictionary['Transcript'].append(transcript_id)
                                        mutated_dictionary['Type'].append("in-frame_no_uORF")
                                        mutated_dictionary['Position'].append(info['pos'][i])
                                        mutated_dictionary['ref'].append(info['ref'][i])
                                        mutated_dictionary['alt'].append(info['alt'][i])
                                        mutated_dictionary['chr'].append(contig)
                                        mutated_dictionary['sample'].append(patient_id)
                    mut_1_1 = dictionary['in-frame_uORF'][len(dictionary['in-frame_uORF']) - 1]
                    if not np.array_equal(mut_1_1, cur_1_1):
                        for pos in np.setdiff1d(np.union1d(mut_1_1, cur_1_1), np.intersect1d(mut_1_1, cur_1_1)):
                            for i in range(0, len(info['pos'])):
                                if (info['pos'][i] <= pos) and (info['pos'][i] > pos - 3):
                                    sum_dictionary['in-frame_uORF_num'][0] = sum_dictionary['in-frame_uORF_num'][
                                                                                 0] + 1
                                    for transcript_id in sample["transcripts"]:
                                        mutated_dictionary['Gene'].append(genes)
                                        mutated_dictionary['Transcript'].append(transcript_id)
                                        mutated_dictionary['Type'].append("in-frame_uORF")
                                        mutated_dictionary['Position'].append(info['pos'][i])
                                        mutated_dictionary['ref'].append(info['ref'][i])
                                        mutated_dictionary['alt'].append(info['alt'][i])
                                        mutated_dictionary['chr'].append(contig)
                                        mutated_dictionary['sample'].append(patient_id)
            else:
                list_many_mut.append((ds[i])["transcripts"][0])

        if save_extra_to_csv:

            pass

            # df0 = pd.DataFrame(data=dictionary)
            # error_list = list_many_mut
            # error_dictionary = {"More_than_14_mutations": error_list}
            # df2 = pd.DataFrame(data=sum_dictionary)
            # df1 = pd.DataFrame(data=error_dictionary)
            #
            # with open(save_to + "/" + patient_id + "_AUGs.csv", 'a') as f:
            #     df0.to_csv(f, header=False)
            #
            # with open(save_to + "/" + patient_id + "_too_many_mutated.csv", 'a') as f:
            #     df1.to_csv(f, header=False)
            #
            # with open(save_to + "/" + patient_id + "_summary.csv", 'a') as f:
            #     df2.to_csv(f, header=False)

    df3 = pd.DataFrame(data=mutated_dictionary)
    df3 = df3.iloc[df3.astype(str).drop_duplicates().index]

    if save_to.endswith('.csv'):
        file = save_to
    else:
        file = save_to + "/mutated.csv"

    with open(file, 'a') as f:
        df3.to_csv(f, header=False)
    pass

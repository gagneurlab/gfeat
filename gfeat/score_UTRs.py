from gfeat.upstreamATG import UpstreamATG
from gfeat.UTR import FivePrimeUTRSeq
from gfeat.units import mutate_sequence_Interval
from gfeat.genome import GFGenome
from gfeat.common_methods import reverse_complement
import pandas as pd


def score_utrs(vcf, gtf, fasta):
    model = UpstreamATG(True, True)  # True: in frame, True: ORF

    data = GFGenome(reference_name='NAME',
                    annotation_name='NAME',
                    gtf_path_or_url=gtf,
                    transcript_fasta_paths_or_urls=fasta,
                    )

    ds = FivePrimeUTRSeq(data, False)  # Todo: check

    dictionary = {'Transcript': [], '0-0': [], '0-1': [], '1-0': [], '1-1': []}
    ref_dictionary = {'Transcript': [], '0-0': [], '0-1': [], '1-0': [], '1-1': []}

    for i in range(1, len(ds)):  # ds[1]
        sample = ds[i]
        ref_seq = ""
        if (list(sample.keys())[0] != "transcripts") and (list(sample.keys())[0] != "exons"):
            ref_seq = list(sample.keys())[0]
        else:
            if (list(sample.keys())[1] != "transcripts") and (list(sample.keys())[1] != "exons"):
                ref_seq = list(sample.keys())[1]
            else:
                ref_seq = list(sample.keys())[2]

        list_exon_mut_seq = []

        for exon in sample["exons"]:
            list_exon_mut_seq.append(mutate_sequence_Interval(exon[1], exon[0], vcf))

        CountOfConbinations = 1
        list_mut_seq = []
        for x in range(0, len(list_exon_mut_seq)):
            CountOfConbinations = CountOfConbinations * len(list_exon_mut_seq[x])
        indexes = [0] * len(list_exon_mut_seq)
        counter = 1
        sequence = ""
        while counter <= CountOfConbinations:
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

        ref_dictionary["Transcript"].append(sample["transcripts"])
        if sample[ref_seq].strand == "-":
            model.predict_on_sample_with_pos_pandas(reverse_complement(ref_seq), ref_dictionary)
        else:
            model.predict_on_sample_with_pos_pandas(ref_seq, ref_dictionary)

        for seq in list_mut_seq:
            dictionary["Transcript"].append(sample["transcripts"])
            model.predict_on_sample_with_pos_pandas(seq, dictionary)


    df1 = pd.DataFrame(data = ref_dictionary)
    df2 = pd.DataFrame(data = dictionary)

    writer = pd.ExcelWriter('output.xlsx')
    df1.to_excel(writer, 'reference')
    df2.to_excel(writer, 'mutated')
    writer.save()

        # mut_seq_list = mutate_sequence(sample[ref_seq_exons], ref_seq_exons, vcf)
        #
        # mut_exon_seq_list = []
        #
        # # print(mut_seq_list)
        #
        # for tuple in mut_seq_list:
        #     mut_exon_seq = ""
        #     for exon in sample["exons"]:
        #         mut_exon_seq = mut_exon_seq + (tuple[0])[(exon[1]).start: (exon[1]).end]
        #         # print(tuple)
        #     mut_exon_seq_list.append((mut_exon_seq, tuple[1]))  # atm: positions of all mutations, not only relevant
        #
        # ref_pred = model.predict_on_sample(ref_seq_exons)
        #
        # for tuple in mut_exon_seq_list:
        #     alt_pred = model.predict_on_sample(tuple[0])
        #     output.append((sample["transcripts"], alt_pred["frame"], alt_pred["uORF"], tuple[1]))

    pass
    ## Append to a pandas DataFrame

from gfeat.upstreamATG import UpstreamATG  # Todo
from gfeat.UTR import FivePrimeUTRSeq
from gfeat.units import mutate_sequence
from gfeat.genome import GFGenome


def score_utrs(vcf, gtf, fasta):
    model = UpstreamATG(True, True)  # True: in frame, True: ORF

    data = GFGenome(reference_name='NAME',
                    annotation_name='NAME',
                    gtf_path_or_url=gtf,
                    transcript_fasta_paths_or_urls=fasta,
                    )

    ds = FivePrimeUTRSeq(data, False)  # Todo: check

    output = []

    for i in range(len(ds)):
        sample = ds[i]
        ref_seq = ""
        if (list(sample.keys())[0] != "transcripts") and (list(sample.keys())[0] != "exons"):
            ref_seq = list(sample.keys())[0]
        mut_seq_list = mutate_sequence(sample[ref_seq], ref_seq, vcf)

        mut_exon_seq_list = []

        # print(mut_seq_list)

        for tuple in mut_seq_list:
            mut_exon_seq = ""
            for exon in sample["exons"]:
                mut_exon_seq = mut_exon_seq + (tuple[0])[(exon[1]).start: (exon[1]).end]
                # print(tuple)
            mut_exon_seq_list.append((mut_exon_seq, tuple[1]))  # atm: positions of all mutations, not only relevant

        ref_pred = model.predict_on_sample(ref_seq)

        for tuple in mut_exon_seq_list:
            alt_pred = model.predict_on_sample(tuple[0])
            output.append((sample["transcripts"], alt_pred["frame"], alt_pred["uORF"], tuple[1]))

    return output
    ## Append to a pandas DataFrame

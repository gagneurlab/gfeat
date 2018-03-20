from gfeat.upstreamATG import UpstreamATG
from gfeat.UTR import FivePrimeUTRSeq
from gfeat.units import mutate_sequence
from gfeat.genome import GFGenome

def score_utrs(vcf, gtf, fasta):
    model = UpstreamATG()

    data = GFGenome(reference_name='NAME',
                     annotation_name='NAME',
                     gtf_path_or_url=gtf,
                     transcript_fasta_paths_or_urls=fasta,
                     )
    ds = FivePrimeUTRSeq(data, False) #Todo: check what's needed - False or TRUE

    for i in range(len(ds)):
        sample = ds[i]
        ref_seq = sample['inputs']
        mut_seq = mutate_sequence(sample["metadata"]["interval"], ref_seq, vcf)
        ref_pred = model.predict_on_sample(ref_seq)
        alt_pred = model.predict_on_sample(alt_seq)
        ## Append to a pandas DataFrame

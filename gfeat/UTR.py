from pybedtools import Interval
from gfeat.genome import GFGenome

class FivePrimeUTRSeq:

    # def __init__(self, gtf, fasta = None): #NOTE: added additional parameter - fasta
    #     """
    #     constructor
    #     :param gtf: path to a gtf file
    #     :return initialises attribute five_prime_utrs – list of strings (5'UTR sequences)
    #     """
    #
    #     gtf_interval = BedTool(gtf)
    #
    #     self.five_prime_utrs = []
    #     if fasta != None:
    #         for interval in gtf_interval.intervals:
    #             if interval.fields[2] == "UTR":
    #                 self.five_prime_utrs.append(interval)
    #
    #     return self.five_prime_utrs

    def __init__(self, data):
        """
        constructor
        :param gtf: path to a gtf file # CHANGE IT IF IT WORKS
        :return initialises attribute five_prime_utrs – list of strings (5'UTR sequences)
        """

        self.five_prime_utrs = []
        self.intervals = []
        self.transcripts = {} # dictionary key – transcript id, value – index of the corresponding 5' UTR

        for contig in data.contigs():
            for transcript in data.transcripts(contig, "+"): # Todo: add "-" strand
                if transcript.contains_start_codon:
                    if transcript.five_prime_utr_sequence != self.five_prime_utrs:
                        self.five_prime_utrs.append(transcript.five_prime_utr_sequence)
                        # BedTool_string = "chr" + contig + " " + transcript.start + " " + \
                        #     + str(transcript.start + transcript.first_start_codon_spliced_offset)
                        # self.intrvals.append(BedTool(BedTool_string, from_string=True))
                        self.intervals.append(Interval(transcript.contig, transcript.start, \
                                                      int(transcript.start + transcript.first_start_codon_spliced_offset), \
                                                      "5' UTR", 0))
                    self.transcripts[transcript.id] = self.five_prime_utrs.index(transcript.five_prime_utr_sequence)
            for transcript in data.transcripts(contig, "-"): # Todo: add "-" strand
                if transcript.contains_start_codon:
                    if transcript.five_prime_utr_sequence != self.five_prime_utrs:
                        self.five_prime_utrs.append(transcript.five_prime_utr_sequence)
                        # BedTool_string = "chr" + contig + " " + transcript.start + " " + \
                        #     + str(transcript.start + transcript.first_start_codon_spliced_offset)
                        # self.intrvals.append(BedTool(BedTool_string, from_string=True))
                        self.intervals.append(interval[transcript.start, transcript.start + \
                                                  transcript.first_start_codon_spliced_offset])
                    self.transcripts[transcript.id] = self.five_prime_utrs.index(transcript.five_prime_utr_sequence)

    def __len__(self):
        """
        :return: the length of the list of 5'UTR sequences – five_prime_utrs
        """
        return len(self.five_prime_utrs)

    def __getitem__(self, idx):
        """
        :param idx: integer index of an element of five_prime_utrs
        :return: dictionary: first entry – utr's sequence,
        second entry – dictionary: first entry – sequence's interval,
                                    second entry – respective transcript's id
        """

        metadata = {}

        for transcript in self.transcripts:
            if self.transcripts[transcript] == idx:
                metadata[transcript] = self.intervals[idx]

        return {
            "inputs": self.five_prime_utrs[idx],
            "metadata" : metadata}
#
# data = GFGenome(reference_name='hg38_test_FivePrimeUTRSeq',
#                 annotation_name='hg38_chr22_test_FivePrimeUTRSeq',
#                 gtf_path_or_url="/Users/veronikakotova/gfeat/tests/data/gencode.v24.annotation_chr22_FivePrimeUTRSeq_testing.gtf",
#                 transcript_fasta_paths_or_urls="/Users/veronikakotova/gfeat/tests/data/hg38_chr22.fa",
#                 )
# FivePrimeUTRSeq_test = FivePrimeUTRSeq(data)
# print(FivePrimeUTRSeq_test[4])

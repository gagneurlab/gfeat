from pybedtools import Interval
from gfeat.genome import GFGenome # Todo
from gfeat.common_methods import reverse_complement # Todo

class FivePrimeUTRSeq:

    def __init__(self, data, reverse_complement_bool):
        """
        constructor
        :param data: GFGenome object
               reverse_complement: bool, True if the reverse_complement of 5'UTR sequence for "-" strand is required
        :return initialises the following attributes:
                    seq[] - list of 5' UTR sequences with exons and introns
                    seq_exons[] - list of 5' UTR sequences with only exons,
                            NOTE it gives exons corresponding to the 5' UTR, therefore the last corresponding exon gets
                             cropped at the start_codon_positions[0]
                    intervals[] - list of 5' UTR intervals
                    transcripts{} - dictionary: key – transcript id, value – index of the corresponding 5' UTR
                    exons{} - dictionary: key - 5' UTR index, value - list of tuples (exon sequence, exon interval)
        """

        # NOTE that pyensembl gives another sequence
        # NOTE 2 5' UTR are considered to be equal only if both their seq and seq_exons are equal
        self.seq = []
        self.seq_exons = []
        self.intervals = []
        self.transcripts = {}
        self.exons = {}

        count = 0

        for contig in data.contigs():
            for transcript in data.transcripts(contig, '+'):
                if transcript.contains_start_codon:

                    temp_seq_exons = ""
                    temp_exon_list = []
                    start = transcript.start_codon_positions[0]
                    for exon in transcript.exons:
                        if exon.start <= transcript.start_codon_positions[0] and start <= exon.start:
                            if exon.end > transcript.start_codon_positions[0]:
                                temp_seq_exons = temp_seq_exons + transcript.sequence[exon.start - transcript.start:
                                                              transcript.start_codon_positions[0] - transcript.start]
                                temp_exon_list.append((transcript.sequence[exon.start - transcript.start:
                                                         transcript.start_codon_positions[0] - transcript.start],
                                                  Interval(transcript.contig, exon.start, transcript.start_codon_positions[0] - 1,
                                                           exon.id, 0, "+")))  # Todo dummy score
                            else:
                                temp_exon_list.append((transcript.sequence[exon.start - transcript.start: exon.end - transcript.start + 1], \
                                                  Interval(transcript.contig, exon.start, exon.end, exon.id, 0, "+")))  # Todo dummy score
                                temp_seq_exons = temp_seq_exons + transcript.sequence[exon.start - transcript.start: \
                                                                                      exon.end - transcript.start + 1]
                            start = exon.start

                    if transcript.sequence[:transcript.start_codon_positions[0] - transcript.start] not in self.seq \
                        and temp_seq_exons not in self.seq_exons:
                        self.seq.append(transcript.sequence[:transcript.start_codon_positions[0] - transcript.start])
                        self.intervals.append(Interval(transcript.contig, transcript.start, int(transcript.start_codon_positions[0] - 1), "5' UTR", 0, "+")) # Todo dummy score

                        self.seq_exons.append(temp_seq_exons)
                        self.exons[count] = temp_exon_list
                        count = count + 1

                    self.transcripts[transcript.id] = self.seq.index(transcript.sequence[:transcript.start_codon_positions[0] - transcript.start])

            for transcript in data.transcripts(contig, '-'):
                if transcript.contains_start_codon:

                    # temp_seq_exons = ""
                    # temp_exon_list = []
                    # end = transcript.end
                    # for exon in transcript.exons:
                    #     if exon.end >= transcript.start_codon_positions[0] and end >= exon.end:
                    #         if exon.start < transcript.start_codon_positions[0]:
                    #             reverse_complemet_seq = reverse_complement(transcript.sequence[transcript.start_codon_positions[2]
                    #                                      - transcript.start + 1:exon.end - transcript.start + 1])
                    #             temp_seq_exons = temp_seq_exons + reverse_complemet_seq
                    #             temp_exon_list.append((reverse_complemet_seq, Interval(transcript.contig,
                    #                   transcript.start_codon_positions[2] + 1, exon.end, exon.id, 0, "-")))  # Todo dummy score
                    #         else:
                    #             reverse_complemet_seq = reverse_complement(transcript.sequence[exon.start - transcript.start: exon.end - transcript.start + 1])
                    #             temp_exon_list.append((reverse_complemet_seq, Interval(transcript.contig, exon.start, exon.end, exon.id, 0, "-")))  # Todo dummy score
                    #             temp_seq_exons = temp_seq_exons + reverse_complemet_seq
                    #         end = exon.end
                    #
                    # curent_transcript_seq = reverse_complement(transcript.sequence[transcript.start_codon_positions[2] - transcript.start + 1:])
                    # if curent_transcript_seq not in self.seq \
                    #     and temp_seq_exons not in self.seq_exons:
                    #     self.seq.append(curent_transcript_seq)
                    #     self.intervals.append(Interval(transcript.contig, transcript.start_codon_positions[2] + 1, transcript.end, "5' UTR", 0, "-")) # Todo dummy score
                    #
                    #     self.seq_exons.append(temp_seq_exons)
                    #     self.exons[count] = temp_exon_list
                    #     count = count + 1
                    #
                    # self.transcripts[transcript.id] = self.seq.index(curent_transcript_seq)

                    temp_seq_exons = ""
                    temp_exon_list = []
                    end = transcript.end
                    for exon in transcript.exons:
                        if exon.end >= transcript.start_codon_positions[0] and end >= exon.end:
                            if exon.start < transcript.start_codon_positions[0]:
                                temp_seq = transcript.sequence[transcript.start_codon_positions[2]
                                                         - transcript.start + 1:exon.end - transcript.start + 1]
                                temp_seq_exons = temp_seq + temp_seq_exons
                                temp_exon_list.append((temp_seq, Interval(transcript.contig,
                                      transcript.start_codon_positions[2] + 1, exon.end, exon.id, 0, "-")))  # Todo dummy score
                            else:
                                temp_seq = transcript.sequence[exon.start - transcript.start: exon.end - transcript.start + 1]
                                temp_exon_list.append((temp_seq, Interval(transcript.contig, exon.start, exon.end, exon.id, 0, "-")))  # Todo dummy score
                                temp_seq_exons = temp_seq + temp_seq_exons
                            end = exon.end

                    if reverse_complement_bool:
                        temp_exon_list_reverse = []
                        temp_seq_exons = reverse_complement(temp_seq_exons)
                        for temp_exon in temp_exon_list:
                            temp_exon_list_reverse.append((reverse_complement(temp_exon[0]), temp_exon[1]))
                        temp_exon_list = temp_exon_list_reverse

                    curent_transcript_seq = reverse_complement(transcript.sequence[transcript.start_codon_positions[2] - transcript.start + 1:])
                    if curent_transcript_seq not in self.seq \
                        and temp_seq_exons not in self.seq_exons:
                        self.seq.append(curent_transcript_seq)
                        self.intervals.append(Interval(transcript.contig, transcript.start_codon_positions[2] + 1, transcript.end, "5' UTR", 0, "-")) # Todo dummy score

                        self.seq_exons.append(temp_seq_exons)
                        self.exons[count] = temp_exon_list
                        count = count + 1

                    self.transcripts[transcript.id] = self.seq.index(curent_transcript_seq)

    def __len__(self):
        """
        :return: the length of the list of 5'UTR sequences – seq
        """
        return len(self.seq)

    def __getitem__(self, idx):
        """
        :param idx: integer index of an element of five_prime_utrs
        :return: dictionary: first entry – dictionary: key - 5' UTR's sequence (with both exons and introns),
                                                       value - interval for 5'UTR
                             second entry – dictionary: key - "transcripts",
                                                        value – list of respective transcripts' ids
                             third entry - dictionary: key - "exons",
                                                        value - list of tuples (exon sequence, interval)
        """

        transcripts = []

        for transcript in self.transcripts:
            if self.transcripts[transcript] == idx:
                transcripts.append(transcript)

        return {
                self.seq[idx]: self.intervals[idx],
                "transcripts": transcripts,
                "exons": self.exons[idx]
            }
#
#
# data_minus = GFGenome(reference_name='hg38_test_FivePrimeUTRSeq_minus',
#                      annotation_name='hg38_chr22_test_FivePrimeUTRSeq_minus',
#                      gtf_path_or_url="/Users/veronikakotova/gfeat/tests/data/gencode.v24.annotation_chr22_FivePrimeUTRSeq_testing_minus-strand.gtf",
#                      transcript_fasta_paths_or_urls="/Users/veronikakotova/gfeat/tests/data/hg38_chr22.fa",
#                      )
# FivePrimeUTRSeq_test_minus = FivePrimeUTRSeq(data_minus, True)
# print(FivePrimeUTRSeq_test_minus.seq)
# print(FivePrimeUTRSeq_test_minus[0])

from pybedtools import Interval
from gfeat.genome import GFGenome

class FivePrimeUTRSeq:

    def __init__(self, data):
        """
        constructor
        :param data: GFGenome object
        :return initialises the following attributes:
                    seq[] - list of 5' UTR sequences with exons and introns
                    seq_exons[] - list of 5' UTR sequences with only exons,
                            NOTE it gives exons corresponding to the 5' UTR, therefore the last corresponding exon gets
                             cropped at the start_codon_positions[0]
                    intervals[] - list of 5' UTR intervals
                    transcripts{} - dictionary: key – transcript id, value – index of the corresponding 5' UTR
                    exons{} - dictionary: key - 5' UTR index, value - list of tuples (exon sequence, exon interval)
        """

        # NOTE that pyensembles gives another sequence
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
                    for exon in transcript.exons:
                        if exon.start <= transcript.start_codon_positions[0]:
                            if exon.end > transcript.start_codon_positions[0]:
                                temp_seq_exons = temp_seq_exons + transcript.sequence[exon.start - transcript.start: \
                                                              transcript.start_codon_positions[0] - transcript.start]
                                temp_exon_list.append((transcript.sequence[exon.start - transcript.start: \
                                                         transcript.start_codon_positions[0] - transcript.start], \
                                                  Interval(transcript.contig, exon.start, transcript.start_codon_positions[0] - 1,\
                                                           exon.id, 0)))  # Todo dummy score
                            else:
                                temp_exon_list.append((transcript.sequence[exon.start - transcript.start: exon.end - transcript.start + 1], \
                                                  Interval(transcript.contig, exon.start, exon.end, exon.id, 0)))  # Todo dummy score
                                temp_seq_exons = temp_seq_exons + transcript.sequence[exon.start - transcript.start: \
                                                                                      exon.end - transcript.start + 1]

                    if transcript.sequence[:transcript.start_codon_positions[0] - transcript.start] not in self.seq \
                        and temp_seq_exons not in self.seq_exons:
                        self.seq.append(transcript.sequence[:transcript.start_codon_positions[0] - transcript.start])
                        self.intervals.append(Interval(transcript.contig, transcript.start, int(transcript.start_codon_positions[0] - 1), "5' UTR", 0)) # Todo dummy score

                        self.seq_exons.append(temp_seq_exons)
                        self.exons[count] = temp_exon_list
                        count = count + 1

                    self.transcripts[transcript.id] = self.seq.index(transcript.sequence[:transcript.start_codon_positions[0] - transcript.start])

            for transcript in data.transcripts(contig, '-'):
                if transcript.contains_start_codon:

                    temp_seq_exons = ""
                    temp_exon_list = []
                    for exon in transcript.exons:
                        if exon.start <= transcript.start_codon_positions[0]:
                            if exon.end > transcript.start_codon_positions[0]:
                                temp_seq_exons = temp_seq_exons + transcript.sequence[exon.start - transcript.start: \
                                                              transcript.start_codon_positions[0] - transcript.start]
                                temp_exon_list.append((transcript.sequence[exon.start - transcript.start: \
                                                         transcript.start_codon_positions[0] - transcript.start], \
                                                  Interval(transcript.contig, exon.start, transcript.start_codon_positions[0] - 1,\
                                                           exon.id, 0)))  # Todo dummy score
                            else:
                                temp_exon_list.append((transcript.sequence[exon.start - transcript.start: exon.end - transcript.start + 1], \
                                                  Interval(transcript.contig, exon.start, exon.end, exon.id, 0)))  # Todo dummy score
                                temp_seq_exons = temp_seq_exons + transcript.sequence[exon.start - transcript.start: \
                                                                                      exon.end - transcript.start + 1]

                    if transcript.sequence[:transcript.start_codon_positions[0] - transcript.start] not in self.seq \
                        and temp_seq_exons not in self.seq_exons:
                        self.seq.append(transcript.sequence[:transcript.start_codon_positions[0] - transcript.start])
                        self.intervals.append(Interval(transcript.contig, transcript.start, int(transcript.start_codon_positions[0] - 1), "5' UTR", 0)) # Todo dummy score

                        self.seq_exons.append(temp_seq_exons)
                        self.exons[count] = temp_exon_list
                        count = count + 1

                    self.transcripts[transcript.id] = self.seq.index(transcript.sequence[:transcript.start_codon_positions[0] - transcript.start])

    def __len__(self):
        """
        :return: the length of the list of 5'UTR sequences – seq
        """
        return len(self.seq)

    def __getitem__(self, idx):
        """
        :param idx: integer index of an element of five_prime_utrs
        :return: dictionary: first entry – 5' UTR's sequence (with both exons and introns),
        second entry – dictionary: first entry – transcripts dictionary:
                                        key - "transcripts",
                                    second entry – respective transcript's id
        """

        metadata = {}

        transcripts = []

        for transcript in self.transcripts:
            if self.transcripts[transcript] == idx:
                metadata[transcript] = self.intervals[idx]
                transcripts.append(transcript)

        return {
            self.seq[idx] : {
                "transcripts" : transcripts,
                "exons" : self.exons[idx]
            }}

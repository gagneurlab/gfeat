from pybedtools import Interval
from gfeat.common_methods import reverse_complement


class FivePrimeUTRSeq:

    def __init__(self, data, reverse_complement_bool, contig=None, strand=None):
        """
        constructor
        :param data: GFGenome or pyensembl.Genome object
               reverse_complement: bool, True if the reverse_complement of 5'UTR sequence for "-" strand is required
               contig: string, number of the chromosome without 'chr'
               strand: character '+' or '-', can be absent
        :return initialises the following attributes:
                    seq[] - list of 5' UTR sequences with exons and introns
                    seq_exons[] - list of 5' UTR sequences with only exons,
                            NOTE it gives exons corresponding to the 5' UTR, therefore the last corresponding exon gets
                             cropped at the start_codon_positions[0]
                    intervals[] - list

                    of 5' UTR intervals
                    transcripts{} - dictionary: key – transcript id, value – index of the corresponding 5' UTR
                    exons{} - dictionary: key - 5' UTR index, value - list of tuples (exon sequence, exon interval)
        """

        # NOTE that sometimes pyensembl 1.1.0 gives a wrong 5'UTR sequence, e.g. ENST00000373388 or ENST00000295591
        # NOTE 2 5' UTR are considered to be equal only if both their seq and seq_exons are equal
        self.seq_exons = []
        self.intervals = []
        self.transcripts = {}
        self.exons = {}

        count = 1

        self.seq_exons.append(" Pyensembl error")
        self.intervals.append(Interval(contig, 0, 0, "5' UTR", 0, "+"))  # Todo dummy score
        self.exons[0] = []

        if strand is None:
            for transcript in data.transcripts(contig, '+'):
                if transcript.contains_start_codon:

                    temp_exon_list = []
                    start = transcript.start
                    start_pos = 0
                    for exon in transcript.exons:
                        if transcript.start_codon_positions[0] >= exon.start >= start:
                            if exon.end > transcript.start_codon_positions[0]:
                                temp_exon_list.append((transcript.five_prime_utr_sequence[start_pos:
                                                                                          (start_pos +
                                                                                           transcript.start_codon_positions[
                                                                                               0] - exon.start)],
                                                       Interval(transcript.contig, exon.start,
                                                                transcript.start_codon_positions[0] - 1,
                                                                exon.id, 0, "+")))  # Todo dummy score
                            else:
                                temp_exon_list.append((transcript.five_prime_utr_sequence[start_pos:
                                                                                          (
                                                                                              start_pos + exon.end - exon.start + 1)],
                                                       Interval(transcript.contig, exon.start, exon.end, exon.id, 0,
                                                                "+")))  # Todo dummy score
                                start_pos = start_pos + exon.end - exon.start + 1
                            start = exon.start

                    # apparently 2 5'UTRs can have the same exonic sequences but different exonic+intronic sequences
                        if transcript.five_prime_utr_sequence not in self.seq_exons:
                            self.seq_exons.append(transcript.five_prime_utr_sequence)
                            self.intervals.append(
                                Interval(transcript.contig, transcript.exons[0].start,
                                         transcript.exons[len(transcript.exons) - 1].end,
                                         "5' UTR", 0, "+"))  # Todo dummy score
                            self.exons[count] = temp_exon_list
                            self.transcripts[transcript.id] = count
                            count = count + 1
                        else:
                            pos = self.seq_exons.index(transcript.five_prime_utr_sequence)
                            if (self.intervals[pos].strand == "+") and \
                                ((self.intervals[pos]).start != (transcript.exons[0]).start):
                                self.seq_exons.append(transcript.five_prime_utr_sequence)
                                self.intervals.append(
                                    Interval(transcript.contig, transcript.exons[0].start,
                                             transcript.exons[len(transcript.exons) - 1].end,
                                             "5' UTR", 0, "+"))  # Todo dummy score
                                self.exons[count] = temp_exon_list
                                self.transcripts[transcript.id] = count
                                count = count + 1
                            else:
                                self.transcripts[transcript.id] = self.seq_exons.index(
                                    transcript.five_prime_utr_sequence)

            for transcript in data.transcripts(contig, '-'):
                if transcript.contains_start_codon:

                    temp_exon_list = []
                    end = transcript.end
                    temp_reverse_seq = reverse_complement(transcript.five_prime_utr_sequence)
                    start_pos = len(temp_reverse_seq)
                    for exon in transcript.exons:
                        if transcript.start_codon_positions[2] <= exon.end <= end:
                            if exon.start < transcript.start_codon_positions[2]:
                                temp_exon_list.append((temp_reverse_seq[:start_pos],
                                                       Interval(transcript.contig,
                                                                transcript.start_codon_positions[2] + 1, exon.end,
                                                                exon.id, 0,
                                                                "-")))  # Todo dummy score
                            else:
                                temp_exon_list.append(
                                    (temp_reverse_seq[start_pos - (exon.end - exon.start + 1):start_pos],
                                     Interval(transcript.contig, exon.start, exon.end, exon.id, 0,
                                              "-")))  # Todo dummy score
                                start_pos = start_pos - (exon.end - exon.start) - 1
                            end = exon.end

                    if reverse_complement_bool:
                        temp_exon_list_reverse = []
                        for temp_exon in temp_exon_list:
                            temp_exon_list_reverse.append((reverse_complement(temp_exon[0]), temp_exon[1]))
                        temp_exon_list = temp_exon_list_reverse
                        curent_transcript_seq = transcript.five_prime_utr_sequence
                    else:
                        curent_transcript_seq = reverse_complement(transcript.five_prime_utr_sequence)

                    if curent_transcript_seq not in self.seq_exons:
                        self.seq_exons.append(curent_transcript_seq)
                        self.intervals.append(Interval(transcript.contig, transcript.exons[0].start,
                                                       transcript.exons[len(transcript.exons) - 1].end, "5' UTR", 0,
                                                       "-"))  # Todo dummy score
                        self.exons[count] = temp_exon_list
                        self.transcripts[transcript.id] = count
                        count = count + 1
                    else:
                        pos = self.seq_exons.index(curent_transcript_seq)
                        if (self.intervals[pos].strand == "-") and \
                            (self.intervals[pos].start != transcript.exons[len(transcript.exons) - 1].end):
                            self.seq_exons.append(curent_transcript_seq)
                            self.intervals.append(
                                Interval(transcript.contig, transcript.exons[0].start,
                                         transcript.exons[len(transcript.exons) - 1].end, "5' UTR", 0,
                                         "-"))  # Todo dummy score
                            self.exons[count] = temp_exon_list
                            self.transcripts[transcript.id] = count
                            count = count + 1
                        else:
                            self.transcripts[transcript.id] = self.seq_exons.index(curent_transcript_seq)
        else:
            if strand == "+":
                for transcript in data.transcripts(contig, '+'):
                    if transcript.contains_start_codon:

                        temp_exon_list = []
                        start = transcript.start
                        start_pos = 0
                        for exon in transcript.exons:
                            if transcript.start_codon_positions[0] >= exon.start >= start:
                                if exon.end > transcript.start_codon_positions[0]:
                                    temp_exon_list.append((transcript.five_prime_utr_sequence[start_pos:
                                                                                              (start_pos +
                                                                                               transcript.start_codon_positions[
                                                                                                   0] - exon.start)],
                                                           Interval(transcript.contig, exon.start,
                                                                    transcript.start_codon_positions[0] - 1,
                                                                    exon.id, 0, "+")))  # Todo dummy score
                                else:
                                    temp_exon_list.append((transcript.five_prime_utr_sequence[start_pos:
                                                                                              (
                                                                                                  start_pos + exon.end - exon.start + 1)],
                                                           Interval(transcript.contig, exon.start, exon.end, exon.id, 0,
                                                                    "+")))  # Todo dummy score
                                    start_pos = start_pos + exon.end - exon.start + 1
                                start = exon.start

                    # apparently 2 5'UTRs can have the same exonic sequences but different exonic+intronic sequences
                        if transcript.five_prime_utr_sequence not in self.seq_exons:
                            self.seq_exons.append(transcript.five_prime_utr_sequence)
                            self.intervals.append(
                                Interval(transcript.contig, transcript.exons[0].start,
                                         transcript.exons[len(transcript.exons) - 1].end,
                                         "5' UTR", 0, "+"))  # Todo dummy score
                            self.exons[count] = temp_exon_list
                            self.transcripts[transcript.id] = count
                            count = count + 1
                        else:
                            pos = self.seq_exons.index(transcript.five_prime_utr_sequence)
                            if (self.intervals[pos].strand == "+") and \
                                ((self.intervals[pos]).start != (transcript.exons[0]).start):
                                self.seq_exons.append(transcript.five_prime_utr_sequence)
                                self.intervals.append(
                                    Interval(transcript.contig, transcript.exons[0].start,
                                             transcript.exons[len(transcript.exons) - 1].end,
                                             "5' UTR", 0, "+"))  # Todo dummy score
                                self.exons[count] = temp_exon_list
                                self.transcripts[transcript.id] = count
                                count = count + 1
                            else:
                                self.transcripts[transcript.id] = self.seq_exons.index(
                                    transcript.five_prime_utr_sequence)

            else:
                for transcript in data.transcripts(contig, '-'):
                    if transcript.contains_start_codon:

                        # if transcript.id == 'ENST00000547691':
                        #     print("herer")

                        temp_exon_list = []
                        end = transcript.end
                        temp_reverse_seq = reverse_complement(transcript.five_prime_utr_sequence)
                        start_pos = len(temp_reverse_seq)
                        for exon in transcript.exons:
                            if transcript.start_codon_positions[2] <= exon.end <= end:
                                if exon.start < transcript.start_codon_positions[2]:
                                    temp_exon_list.append((temp_reverse_seq[:start_pos],
                                                           Interval(transcript.contig,
                                                                    transcript.start_codon_positions[2] + 1, exon.end,
                                                                    exon.id, 0,
                                                                    "-")))  # Todo dummy score
                                else:
                                    temp_exon_list.append(
                                        (temp_reverse_seq[start_pos - (exon.end - exon.start + 1):start_pos],
                                         Interval(transcript.contig, exon.start, exon.end, exon.id, 0,
                                                  "-")))  # Todo dummy score
                                    start_pos = start_pos - (exon.end - exon.start) - 1
                                end = exon.end

                        if reverse_complement_bool:
                            temp_exon_list_reverse = []
                            for temp_exon in temp_exon_list:
                                temp_exon_list_reverse.append((reverse_complement(temp_exon[0]), temp_exon[1]))
                            temp_exon_list = temp_exon_list_reverse
                            curent_transcript_seq = transcript.five_prime_utr_sequence
                        else:
                            curent_transcript_seq = reverse_complement(transcript.five_prime_utr_sequence)

                        if curent_transcript_seq not in self.seq_exons:
                            self.seq_exons.append(curent_transcript_seq)
                            self.intervals.append(Interval(transcript.contig, transcript.start_codon_positions[2] + 1,
                                                           transcript.end, "5' UTR", 0, "-"))  # Todo dummy score
                            self.exons[count] = temp_exon_list
                            self.transcripts[transcript.id] = count
                            count = count + 1
                        else:
                            pos = self.seq_exons.index(curent_transcript_seq)
                            if (self.intervals[pos].strand == "-") and \
                                (self.intervals[pos].start != transcript.exons[len(transcript.exons) - 1].end):
                                self.seq_exons.append(curent_transcript_seq)
                                self.intervals.append(
                                    Interval(transcript.contig, transcript.start_codon_positions[2] + 1, transcript.end,
                                             "5' UTR", 0, "-"))  # Todo dummy score
                                self.exons[count] = temp_exon_list
                                self.transcripts[transcript.id] = count
                                count = count + 1
                            else:
                                self.transcripts[transcript.id] = self.seq_exons.index(curent_transcript_seq)

    def __len__(self):
        """
        :return: the length of the list of 5'UTR sequences – seq
        """
        return len(self.seq_exons)

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
            self.seq_exons[idx]: self.intervals[idx],
            "transcripts": transcripts,
            "exons": self.exons[idx]
        }

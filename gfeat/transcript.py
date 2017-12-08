import pyensembl
#from .units import IndexUnit

# Boyer Moore string search algorithm
def boyer_moore_search(text, pattern):
    occurrences = dict()
    for letter in {"A", "B", "C", "D"}:
        occurrences[letter] = pattern.rfind(letter)
    m = len(pattern)
    n = len(text)
    itext = m - 1
    ipattern = m - 1
    count = 0
    while itext < n:
        if text[itext] == pattern[ipattern]:
            if ipattern == 0:
                count += 1
                itext += m - 1
            else:
                itext -= 1
                ipattern -= 1
        else:
            l = occurrences[letter]
            itext = itext + m - min(ipattern, 1+l)
            ipattern = m - 1
    return -1

#class GFTranscript(IndexUnit, pyensembl.Transcript):
# Transcript class inherited from pyensembl's Transcript class
class GFTranscript(pyensembl.Transcript):

    # Init already implemented by `pyensembl.Transcript`
    # TODO  implement the feature extractors

    # Counts how many codons are in the transcript's sequence
    def codon_count(self):
        # Removing 5' UTR and 3' UTR sequences as they don't have any codons
        sequence = self.sequence.replace(self.five_prime_utr_sequence, "").replace(self.three_prime_utr_sequence, "")
        return len(sequence)/3

    # Counts how many times a given motif is presented in the 3' UTR sequence
    def utr3_motif_counts(self, pattern):
        return boyer_moore_search(self.three_prime_utr_sequence.upper(), pattern.upper())

    # Counts how many times a given motif is presented in the 5' UTR sequence
    def utr5_motif_counts(self, pattern):
        return boyer_moore_search(self.five_prime_utr_sequence.upper(), pattern.upper())

    @classmethod
    def from_pyensembl(cls, obj):
        pass

    @classmethod
    def iter_all(cls, genome):
        pass


    # TODO - implement
    #         - codon_counts
    #         - utr3_motif_counts()
    #            - unit-test each of them
    #               - deposit a small fasta and gtf (chr22) to the test files

# Loading testing data


import pyensembl
import unicodedata
#from .units import IndexUnit

# Boyer Moore string search algorithm
def boyer_moore_search(text, pattern):
    # alphabet = set(text)
    occurrences = dict()
    for letter in {'A', 'C', 'G', 'T'}:
        occurrences[letter] = pattern.rfind(letter)
    #last = last_occurrence(pattern, alphabet)
    m = len(pattern)
    n = len(text)
    i = m - 1  # text index
    j = m - 1  # pattern index
    count = 0;
    while i < n:
        if text[i] == pattern[j]:
            if j == 0:
                count += 1
                #print text[i-1]
                print i
                print "---" + str(count)
                i += 2*m - 1
                j = m - 1
                #return i
            else:
                print i
                i -= 1
                j -= 1
        else:
            l = occurrences[text[i]]
            #print i
            #print l
            #print text[i]
            print i
            i = i + m - min(j, 1 + l)
            #print i
            #print "-----------------------"
            j = m - 1
    # print count
    return count

#class GFTranscript(IndexUnit, pyensembl.Transcript):
# Transcript class inherited from pyensembl's Transcript class
class GFTranscript(pyensembl.Transcript):

    # Init already implemented by `pyensembl.Transcript`
    # TODO  implement the feature extractors

    # Counts how many codons are in the transcript's sequence
    def codon_counts(self):
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
data = pyensembl.ensembl_release.EnsemblRelease(75)
mytranscript = GFTranscript("ENST00000369985", "MYO6-001", "6", 76458926, 76629253, "+", "protein_coding", "ENSG00000196586", data)
#print mytranscript.codon_counts()
print mytranscript.three_prime_utr_sequence
print mytranscript.utr3_motif_counts("ACCA")

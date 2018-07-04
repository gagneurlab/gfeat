import pyensembl
from pyensembl import Genome


from six import string_types

from pyensembl import memory_cache
from pyensembl import download_cache
from pyensembl.common import dump_pickle
from gfeat.transcript import GFTranscript
import pandas as pd


class GFGenome(pyensembl.Genome):

    def __init__(
            self,
            reference_name: object,
            annotation_name: object,
            annotation_version: object = None,
            gtf_path_or_url: object = None,
            transcript_fasta_paths_or_urls: object = None,
            protein_fasta_paths_or_urls: object = None,
            decompress_on_download: object = False,
            copy_local_files_to_cache: object = False,
            require_ensembl_ids: object = True,
            cache_directory_path: object = None):

        super(GFGenome, self).__init__(reference_name,
            annotation_name,
            annotation_version,
            gtf_path_or_url,
            transcript_fasta_paths_or_urls,
            protein_fasta_paths_or_urls,
            decompress_on_download,
            copy_local_files_to_cache,
            require_ensembl_ids,
            cache_directory_path)
        # first_dict_key = 0;
        key_values = []

        # parsing of ordinary ??? Fasta files
        if ((len(self.transcript_sequences.fasta_dictionary) <= 1) and (self.transcript_fasta_path is not None)):
            for key in self.transcript_sequences.fasta_dictionary.keys():
                key_values.append(key);
            for key in key_values:
                for trans in self.transcripts():
                    self.transcript_sequences.fasta_dictionary[trans.transcript_id] = self.transcript_sequences.fasta_dictionary[key][trans.start-1:trans.end];
        dump_pickle(self._transcript_sequences._fasta_dictionary, self._transcript_sequences.fasta_dictionary_pickle_path)

    def check_fasta_dictionary(self):

        pass

    def transcripts(self, contig=None, strand=None):
        """
        Construct Transcript object for every transcript entry in
        the database. Optionally restrict to a particular
        chromosome using the `contig` argument.
        """
        transcript_ids = self.transcript_ids(contig=contig, strand=strand)
        return [
            GFTranscript(copy_transcript = self.transcript_by_id(transcript_id))
            for transcript_id in transcript_ids
        ]

    def get_consensus_Kozak_seq(self, seq = False):
        """
        TODO

        BY DEFAULT RETURN NUMBERS

        :param pattern: string motif to be found in the 5' UTR sequence
        :return: how many times a given motif is presented in the 5' UTR sequence
        """

        nucleobase_count = {"0A": [0], "0C": [0], "0G": [0], "0T": [0],
                            "1A": [0], "1C": [0], "1G": [0], "1T": [0],
                            "2A": [0], "2C": [0], "2G": [0], "2T": [0],
                            "3A": [0], "3C": [0], "3G": [0], "3T": [0],
                            "4A": [0], "4C": [0], "4G": [0], "4T": [0],
                            "5A": [0], "5C": [0], "5G": [0], "5T": [0],
                            "6A": [0], "6C": [0], "6G": [0], "6T": [0],
                            "7A": [0], "7C": [0], "7G": [0], "7T": [0],
                            "8A": [0], "8C": [0], "8G": [0], "8T": [0],
                            "9A": [0], "9C": [0], "9G": [0], "9T": [0],
                            "10A": [0], "10C": [0], "10G": [0], "10T": [0],
                            "11A": [0], "11C": [0], "11G": [0], "11T": [0],
                            "12A": [0], "12C": [0], "12G": [0], "12T": [0],
                            "13A": [0], "13C": [0], "13G": [0], "13T": [0],
                            "14A": [0], "14C": [0], "14G": [0], "14T": [0]}

        for contig in self.contigs():
            for transcript in self.transcripts(contig, '+'):
                if transcript.contains_start_codon and transcript.contains_stop_codon:
                    i = 0
                    Kozak_seq = transcript.get_Kozak_seq()
                    for c in Kozak_seq:
                        nucleobase_count[str(i)+c][0] = nucleobase_count[str(i)+c][0] + 1
                        i = i + 1
            for transcript in self.transcripts(contig, '-'):
                if transcript.contains_start_codon and transcript.contains_stop_codon:
                    i = 0
                    Kozak_seq = transcript.get_Kozak_seq()
                    for c in Kozak_seq:
                        nucleobase_count[str(i)+c][0] = nucleobase_count[str(i)+c][0] + 1
                        i = i + 1

        nucliobase_int = {0: "A", 1: "C", 2: "G", 3: "T"}
        consensus_seq = ""

        if seq:
            for i in range(15):
                temp_list = [nucleobase_count[str(i)+"A"][0],
                             nucleobase_count[str(i) + "C"][0],
                             nucleobase_count[str(i) + "G"][0],
                             nucleobase_count[str(i) + "T"][0]]
                consensus_seq = consensus_seq + nucliobase_int[temp_list.index(max(temp_list))]
        else:
            for i in range(15):
                temp_list = [nucleobase_count[str(i)+"A"][0],
                             nucleobase_count[str(i) + "C"][0],
                             nucleobase_count[str(i) + "G"][0],
                             nucleobase_count[str(i) + "T"][0]]
                consensus_seq = consensus_seq + str(temp_list.index(max(temp_list)))

        return consensus_seq

    def get_Kozak_matrix(self):
        """
        TODO
        :param pattern: string motif to be found in the 5' UTR sequence
        :return: how many times a given motif is presented in the 5' UTR sequence
        """

        lines = []
        for contig in self.contigs():
            for transcript in self.transcripts(contig, '+'):
                if transcript.contains_start_codon and transcript.contains_stop_codon:
                    lines.append(transcript.get_line_Kozak_matrix())
            for transcript in self.transcripts(contig, '-'):
                if transcript.contains_start_codon and transcript.contains_stop_codon:
                    lines.append(transcript.get_line_Kozak_matrix())

        df_Kozak_matrix = pd.concat(lines)
        consensus_Kozak_seq = self.get_consensus_Kozak_seq(True)

        i = 0
        for c in consensus_Kozak_seq:
            df_Kozak_matrix = df_Kozak_matrix.drop(columns=[str(i)+c])
            i = i + 1

        return df_Kozak_matrix

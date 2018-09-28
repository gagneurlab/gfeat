import pyensembl
from pyensembl.common import dump_pickle
from gfeat.transcript import GFTranscript
from gfeat.utils import PCA_with_standard_sample_deviation_scaling, reverse_complement
from gfeat.units import VCFMutator
import pandas as pd
from pybedtools import Interval


class GFGenome(pyensembl.Genome):
    """
    GFGenome class
    """

    def __init__(
        self,
        reference_name: object = None,
        annotation_name: object = None,
        annotation_version: object = None,
        gtf_path_or_url: object = None,
        transcript_fasta_paths_or_urls: object = None,
        protein_fasta_paths_or_urls: object = None,
        decompress_on_download: object = False,
        copy_local_files_to_cache: object = False,
        require_ensembl_ids: object = True,
        cache_directory_path: object = None,
        copy_genome=None):

        if copy_genome is not None:
            reference_name = copy_genome.reference_name
            annotation_name = copy_genome.annotation_name
            annotation_version = copy_genome.annotation_version
            gtf_path_or_url = copy_genome._gtf_path_or_url
            transcript_fasta_paths_or_urls = copy_genome._transcript_fasta_path_or_url
            protein_fasta_paths_or_urls = copy_genome._protein_fasta_path_or_url
            decompress_on_download = copy_genome.decompress_on_download
            copy_local_files_to_cache = copy_genome.copy_local_files_to_cache
            require_ensembl_ids = copy_genome.require_ensembl_ids
            cache_directory_path = copy_genome.cache_directory_path
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

        if copy_genome is None:
            key_values = []
            if (len(self.transcript_sequences.fasta_dictionary) <= 1) and (self.transcript_fasta_path is not None):
                for key in self.transcript_sequences.fasta_dictionary.keys():
                    key_values.append(key);
                for key in key_values:
                    for trans in self.transcripts():
                        self.transcript_sequences.fasta_dictionary[trans.transcript_id] = \
                            self.transcript_sequences.fasta_dictionary[key][trans.start - 1:trans.end];
            dump_pickle(self._transcript_sequences._fasta_dictionary,
                        self._transcript_sequences.fasta_dictionary_pickle_path)

    def check_fasta_dictionary(self):

        pass

    def transcripts(self, contig=None, strand=None):
        """
        Get GFTranscript objects. Optionally restrict to a particular chromosome and strand

        :param contig: optional, chromosome transcripts for which are required
        :type contig: int
        :param strand: optional, chromosome strand transcripts for which are required
        :type strand: str
        :return: list of gfeat.transcript.GFTranscript
        """

        transcript_ids = self.transcript_ids(contig=contig, strand=strand)
        return [
            GFTranscript(copy_transcript=self.transcript_by_id(transcript_id))
            for transcript_id in transcript_ids
        ]

    def transcripts_by_name(self, transcript_name):
        """
        Get GFTranscript objects by the specified transcript name

        :param transcript_name: name of the transcripts
        :type transcript_name: str
        :return: gfeat.transcript.GFTranscript object
        """
        transcripts = self.transcripts_by_name(transcript_name)
        return [
            GFTranscript(copy_transcript=self.transcript_by_id(transcript))
            for transcript in transcripts
        ]

    def transcript_by_protein_id(self, protein_id):
        """
        Get GFTranscript objects by the specified protein id

        :param protein_id: protein id transcript for which is required
        :type protein_id: str
        :return: gfeat.transcript.GFTranscript object
        """
        return GFTranscript(copy_transcript=self.transcript_by_id(protein_id))

    def gftranscript_by_id(self, transcript_id):
        """
        Get GFTranscript object by the transcript id

        :param transcript_id: id of the transcript
        :type transcript_id: str
        :return: gfeat.transcript.GFTranscript object
        """
        return GFTranscript(copy_transcript=self.transcript_by_id(transcript_id))

    def get_consensus_Kozak_seq(self, seq=False):
        """
        Get the consensus Kozak sequence (Kozak sequence which is present more times than all other for this being)

        :param seq: True – to return the letter representation of the sequence

            False – to return the digit representation of the sequence: 0 – A, 1 – C, 2 – G, 3 – T

        :type seq: bool
        :return: string, consensus Kozak sequence constituting of either digits or letters depending on the specified
		seq parameter
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
                        nucleobase_count[str(i) + c][0] = nucleobase_count[str(i) + c][0] + 1
                        i = i + 1
            for transcript in self.transcripts(contig, '-'):
                if transcript.contains_start_codon and transcript.contains_stop_codon:
                    i = 0
                    Kozak_seq = transcript.get_Kozak_seq()
                    for c in Kozak_seq:
                        nucleobase_count[str(i) + c][0] = nucleobase_count[str(i) + c][0] + 1
                        i = i + 1

        nucleobase_int = {0: "A", 1: "C", 2: "G", 3: "T"}
        consensus_seq = ""

        if seq:
            for i in range(15):
                temp_list = [nucleobase_count[str(i) + "A"][0],
                             nucleobase_count[str(i) + "C"][0],
                             nucleobase_count[str(i) + "G"][0],
                             nucleobase_count[str(i) + "T"][0]]
                consensus_seq = consensus_seq + nucleobase_int[temp_list.index(max(temp_list))]
        else:
            for i in range(15):
                temp_list = [nucleobase_count[str(i) + "A"][0],
                             nucleobase_count[str(i) + "C"][0],
                             nucleobase_count[str(i) + "G"][0],
                             nucleobase_count[str(i) + "T"][0]]
                consensus_seq = consensus_seq + str(temp_list.index(max(temp_list)))

        return consensus_seq

    def get_Kozak_matrix(self):
        """
        Get Kozak matrix with all transcripts for this being. Consensus columns are removed

        :return: pandas.DataFrame,

                column names – first number corresponds to the base position in the Kozak sequence,
                e.g. 0A means base A 5 bases upstream from the start codon

                rows – 1 if it has the corresponding base, 0 otherwise

        NOTE: columns corresponding to the most frequent bases at each position are removed

        """

        lines = []
        for contig in self.contigs():
            for transcript in self.transcripts(contig, '+'):
                if transcript.contains_start_codon and transcript.contains_stop_codon:
                    lines.append(transcript.get_Kozak_seq_as_df())
            for transcript in self.transcripts(contig, '-'):
                if transcript.contains_start_codon and transcript.contains_stop_codon:
                    lines.append(transcript.get_Kozak_seq_as_df())

        df_Kozak_matrix = pd.concat(lines)
        consensus_Kozak_seq = self.get_consensus_Kozak_seq(True)

        i = 0
        for c in consensus_Kozak_seq:
            df_Kozak_matrix = df_Kozak_matrix.drop(columns=[str(i) + c])
            i = i + 1

        return df_Kozak_matrix

    def get_consensus_stop_codon_context(self, seq=False):
        """
        Get the consensus stop codon context sequence (stop codon context sequence which is present more times than
        all other for this being)

        :param seq: True – to return the letter representation of the sequence

                          False – to return the digit representation of the sequence: 0 – A, 1 – C, 2 – G, 3 – T
        :type seq: bool
        :return: string, consensus stop codon context sequence constituting of either digits or letters
		depending on the specified seq parameter
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
                    Kozak_seq = transcript.get_stop_codon_context()
                    for c in Kozak_seq:
                        nucleobase_count[str(i) + c][0] = nucleobase_count[str(i) + c][0] + 1
                        i = i + 1
            for transcript in self.transcripts(contig, '-'):
                if transcript.contains_start_codon and transcript.contains_stop_codon:
                    i = 0
                    Kozak_seq = transcript.get_stop_codon_context()
                    for c in Kozak_seq:
                        nucleobase_count[str(i) + c][0] = nucleobase_count[str(i) + c][0] + 1
                        i = i + 1

        nucleobase_int = {0: "A", 1: "C", 2: "G", 3: "T"}
        consensus_seq = ""

        if seq:
            for i in range(15):
                temp_list = [nucleobase_count[str(i) + "A"][0],
                             nucleobase_count[str(i) + "C"][0],
                             nucleobase_count[str(i) + "G"][0],
                             nucleobase_count[str(i) + "T"][0]]
                consensus_seq = consensus_seq + nucleobase_int[temp_list.index(max(temp_list))]
        else:
            for i in range(15):
                temp_list = [nucleobase_count[str(i) + "A"][0],
                             nucleobase_count[str(i) + "C"][0],
                             nucleobase_count[str(i) + "G"][0],
                             nucleobase_count[str(i) + "T"][0]]
                consensus_seq = consensus_seq + str(temp_list.index(max(temp_list)))

        return consensus_seq

    def get_stop_codon_context_matrix(self):
        """
        Get stop codon context matrix with all transcripts for this being. Consensus columns are removed

        :return: pandas.DataFrame,

                column names – first number corresponds to the base position in the stop codon
                context sequence, e.g. 0A means base A 5 bases upstream from the start codon

                rows – 1 if it has the corresponding base, 0 otherwise

        NOTE: columns corresponding to the most frequent bases at each position are removed
        """

        lines = []
        for contig in self.contigs():
            for transcript in self.transcripts(contig, '+'):
                if transcript.contains_start_codon and transcript.contains_stop_codon:
                    lines.append(transcript.get_stop_codon_context_as_df())
            for transcript in self.transcripts(contig, '-'):
                if transcript.contains_start_codon and transcript.contains_stop_codon:
                    lines.append(transcript.get_stop_codon_context_as_df())

        df_stop_codon_context_matrix = pd.concat(lines)
        consensus_stop_codon_context = self.get_consensus_stop_codon_context(True)

        i = 0
        for c in consensus_stop_codon_context:
            df_stop_codon_context_matrix = df_stop_codon_context_matrix.drop(columns=[str(i) + c])
            i = i + 1

        return df_stop_codon_context_matrix

    def get_codon_pair_bias(self):
        """
        Get truncated PCA of the codon frequency matrix with 2 principal components

        :return: pandas.DataFrame, truncated PCA of the codon frequency matrix with 2 principal components
        """

        lines = []

        for contig in self.contigs():
            for transcript in self.transcripts(contig, '+'):
                if transcript.contains_start_codon and transcript.contains_stop_codon:
                    lines.append(transcript.get_codon_pairs_frequency())
            for transcript in self.transcripts(contig, '-'):
                if transcript.contains_start_codon and transcript.contains_stop_codon:
                    lines.append(transcript.get_codon_pairs_frequency())

        df_codon_pairs_frequency = pd.concat(lines)
        PCA_decomposition = PCA_with_scaling(df_codon_pairs_frequency)
        return PCA_decomposition

    def get_nucleobase_mutation_table(self, vcf):
        """
        Get a table which shows whether a certain nucleobase in Kozak sequence or stop codon context was mutated or not.

        :param vcf: path to the vcf.gz or file opened using cyvcf2
        :type vcf: string or an "opened" file
        :return: pd.DataFrame,

                column names – K_i – where i shows position in Kozak sequence;
                S_i – where i shows position in stop codon context; gene_id; transcript_id

                rows – NaN – no variant, 1 – heterozygous variant, 2 – homozygous variant
        """

        # only for Kozak sequence and stop codon context + transcript_id column
        columns = ["K_0", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "K_11", "K_12", "K_13",
                   "K_14",
                   "S_0", "S_1", "S_2", "S_3", "S_4", "S_5", "S_6", "S_7", "S_8", "S_9", "S_10", "S_11", "S_12", "S_13",
                   "S_14",
                   "transcript_id", "gene_id", "name"]
        df_nucleobases = pd.DataFrame(columns=columns)
        nucleobases_lines = []

        mutator = VCFMutator(False, True, vcf, True)

        contigs = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18',
                   '19', '20', '21', '22', 'X', 'Y']

        for contig in contigs:
            for transcript in self.transcripts(contig, '+'):
                if transcript.contains_stop_codon and transcript.contains_start_codon:
                    Kozak_seq = transcript.get_Kozak_seq()
                    Interval_Kozak = Interval("chr" + transcript.contig, transcript.start_codon_positions[0] - 6,
                                              transcript.start_codon_positions[0] + 9, "NA", 0, "+")
                    stop_codon_context = transcript.get_stop_codon_context()
                    Interval_stop = Interval("chr" + transcript.contig, transcript.stop_codon_positions[0] - 6,
                                             transcript.stop_codon_positions[0] + 9, "NA", 0, "+")
                    df_nucleobases_line = mutator.mutate_codon_context([Interval_Kozak, Interval_stop],
                                                                       [Kozak_seq, stop_codon_context], ["K_", "S_"])
                    if len(Kozak_seq) < 15:
                        new_columns = []
                        for column in df_nucleobases_line:
                            if column.find("K_") != -1:
                                new_columns.append("K_" + str(int(column[2:]) + (15 - len(Kozak_seq))))
                            else:
                                new_columns.append(column)
                        df_nucleobases_line.columns = new_columns

                    df_nucleobases_line["transcript_id"] = transcript.id
                    df_nucleobases_line["gene_id"] = transcript.gene_id
                    nucleobases_lines.append(df_nucleobases_line)
            for transcript in self.transcripts(contig, '-'):
                if transcript.contains_stop_codon and transcript.contains_start_codon:
                    Kozak_seq = reverse_complement(transcript.get_Kozak_seq())
                    Interval_Kozak = Interval("chr" + transcript.contig, transcript.start_codon_positions[0] - 6,
                                              transcript.start_codon_positions[0] + 9, "NA", 0, "-")
                    stop_codon_context = reverse_complement(transcript.get_stop_codon_context())
                    Interval_stop = Interval("chr" + transcript.contig, transcript.stop_codon_positions[0] - 6,
                                             transcript.stop_codon_positions[0] + 9, "NA", 0, "-")
                    df_nucleobases_line = mutator.mutate_codon_context([Interval_Kozak, Interval_stop],
                                                                       [Kozak_seq, stop_codon_context], ["K_", "S_"])
                    if len(Kozak_seq) < 15:
                        new_columns = []
                        for column in df_nucleobases_line:
                            if column.find("K_") != -1:
                                new_columns.append("K_" + str(int(column[2:]) + (15 - len(Kozak_seq))))
                            else:
                                new_columns.append(column)
                        df_nucleobases_line.columns = new_columns

                    df_nucleobases_line["transcript_id"] = transcript.id
                    df_nucleobases_line["gene_id"] = transcript.gene_id
                    nucleobases_lines.append(df_nucleobases_line)

            df_nucleobases = pd.concat(nucleobases_lines, ignore_index=True)
            df_nucleobases = df_nucleobases.drop(['name'], axis=1)
        return df_nucleobases

import pyensembl
from pyensembl import Genome


from six import string_types

from pyensembl import memory_cache
from pyensembl import download_cache
from pyensembl.common import dump_pickle


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

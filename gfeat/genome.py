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

        # if isinstance(transcript_fasta_paths_or_urls, string_types):
        #     transcript_fasta_paths_or_urls = [transcript_fasta_paths_or_urls]
        #
        # if isinstance(protein_fasta_paths_or_urls, string_types):
        #     protein_fasta_paths_or_urls = [protein_fasta_paths_or_urls]
        #
        # self.reference_name = reference_name
        # self.annotation_name = annotation_name
        # self.annotation_version = annotation_version
        # self.decompress_on_download = decompress_on_download
        # self.copy_local_files_to_cache = copy_local_files_to_cache
        # self.require_ensembl_ids = require_ensembl_ids
        # self.cache_directory_path = cache_directory_path
        # self._gtf_path_or_url = gtf_path_or_url
        # self._transcript_fasta_paths_or_urls = transcript_fasta_paths_or_urls
        # self._protein_fasta_paths_or_urls = protein_fasta_paths_or_urls
        #
        # self.download_cache = pyensembl.DownloadCache(
        #     reference_name=self.reference_name,
        #     annotation_name=self.annotation_name,
        #     annotation_version=self.annotation_version,
        #     decompress_on_download=self.decompress_on_download,
        #     copy_local_files_to_cache=self.copy_local_files_to_cache,
        #     install_string_function=self.install_string,
        #     cache_directory_path=self.cache_directory_path)
        # self.cache_directory_path = self.download_cache.cache_directory_path
        #
        # self.has_gtf = self._gtf_path_or_url is not None
        # self.has_transcript_fasta = self._transcript_fasta_paths_or_urls is not None
        # self.has_protein_fasta = self._protein_fasta_paths_or_urls is not None
        # self.memory_cache = pyensembl.MemoryCache()
        #
        # self._init_lazy_fields()





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

        if ((len(self.transcript_sequences.fasta_dictionary) <= 1) and (self.transcript_fasta_path is not None)):
            for key in self.transcript_sequences.fasta_dictionary.keys():
                key_values.append(key);
            for key in key_values:
                for trans in self.transcripts():
                    self.transcript_sequences.fasta_dictionary[trans.transcript_id] = self.transcript_sequences.fasta_dictionary[key][trans.start-1:trans.end];
        dump_pickle(self._transcript_sequences._fasta_dictionary, self._transcript_sequences.fasta_dictionary_pickle_path)


    def check_fasta_dictionary(self):

        pass

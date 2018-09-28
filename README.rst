===============================
gfeat
===============================


.. image:: https://img.shields.io/pypi/v/gfeat.svg
        :target: https://pypi.python.org/pypi/gfeat

.. image:: https://img.shields.io/travis/avsecz/gfeat.svg
        :target: https://travis-ci.org/avsecz/gfeat

.. image:: https://readthedocs.org/projects/gfeat/badge/?version=latest
        :target: https://gfeat.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://pyup.io/repos/github/avsecz/gfeat/shield.svg
     :target: https://pyup.io/repos/github/avsecz/gfeat/
     :alt: Updates


Python genomic features extractor from raw files (Fasta, GTF and VCF).

*gfeat* is a convenient tool for extracting diverse DNA features for genomic modeling and analysis. It allows to get,
for example, transcript Kozak sequence and GC content, genome codon pair bias and etc. It is also possible to inject
SNPs using a VCF file, get position-type variant matrices and investigate sequences with different combinations
of heterozygous variants (when homozygous always stay injected in the sequence).

You can find several short usage examples in the Tutorial section of the documentation.

* Free software: MIT license
* Documentation: https://i12g-gagneurweb.in.tum.de/public/docs/gfeat


Features
--------

GFGenome:

* GFGenome.get_consensus_Kozak_seq(seq=False)
    Consensus Kozak sequence

* GFGenome.get_Kozak_matrix()
    Kozak sequence matrix (rows - transcripts)

* GFGenome.get_consensus_stop_codon_context(seq=False)
    Consensus stop codon context sequence

* GFGenome.get_stop_codon_context_matrix()
    Stop codon context matrix (rows - transcripts)

* GFGenome.get_codon_pair_bias()
    Codon pair bias

* GFGenome.get_nucleobase_mutation_table()
    Nucleobase mutation matrix (table with positions of variants and their type in Kozak sequence and stop codon context)

GFTranscript:

* GFTranscript.codon_counts()
    Coding sequence codon count

* GFTranscript.utr5_motif_counts(pattern)
    5'UTR motif count

* GFTranscript.utr3_motif_counts(pattern)
    3'UTR motif count

* GFTranscript.codon_usage()
    Coding sequence codon usage

* GFTranscript.gc_content(region)
    Coding sequence, 5'UTR or 3'UTR G and C content

* GFTranscript.get_Kozak_seq()
    Kozak sequence

* GFTranscript.get_stop_codon_context()
    Stop codon context

* GFTranscript.get_codon_pairs_frequency()
    Coding sequence codon pair frequency

Upstream AUG:

* UpstreamAUG.predict_on_sample(seq)
    Predict on sample (get all AUGs with the information whether they are in-frame or not and whether they have a stop codon or not)

* UpstreamAUG.predict_on_sample_with_pos(seq)
    Predict on sample with position (predict on sample plus positions of AUGs)

* UpstreamAUG.predict_on_sample_with_pos_pandas()
    Predict on sample with position and appending to the passed dictionary

* UpstreamAUG.predict_on_batch(seq_list)
    Predict on batch

5'UTR class:

* FivePrimeUTRSeq(data, reverse_complement_bool, contig=None, strand=None)
    An object contains list of all 5'UTR sequences, their positions, exons, positions of exons and corresponding transcripts

Auxiliary functions:

* VCFMutator.mutate_sequence(interval, fasta=None, seq_whole=None)
    Mutate sequence

* VCFMutator.mutate_codon_context(intervals, seqs, column_names)
    Mutate codon context

* reverse_complement(dna)
    Reverse complement

* PCA_with_standard_sample_deviation_scaling(df, n_comp=2)
    Pca with scaling

Credits
---------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template. *gfeat* is based
on pyensembl_ package and can be partially viewed as its extension.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _pyensembl: https://github.com/openvax/pyensembl

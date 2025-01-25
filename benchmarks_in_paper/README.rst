Contents
--------
Benchmark data for the TIPP3 paper.

Overview
~~~~~~~~
This directory contains all benchmarks obtained for the TIPP3 paper, including
the unprocessed abundance profiles produced by TIPP3, TIPP3-fast, Kraken2,
Bracken, Metabuli, MetaPhlAn4, and mOTUsv3.

Kraken2, Bracken, and Metabuli were
run with two types of inputs: all reads, or just the filtered reads by TIPP3.
In the paper, we reported that filtering generally helps with the profiling
accuracy (except for Bracken, which under hard conditions we observed worse
performance).

Fold Structure
~~~~~~~~~~~~~~
``PROFILES.tar.gz`` contains all the raw benchmark results of each method, for
each different experiment. A general structure of the abundance profile at the
7 taxonomic levels (species, genus, family, order, class, phylum, superkingdom)
is as follows::

   PROFILES
     |_<experiment_id>
       |_<method_id>
         |_<method_id>.abundance.species.tsv
         |_<method_id>.abundance.genus.tsv
         |_...

Experiment IDs
++++++++++++++
#. ``cami_long``: CAMI-II Marine PacBio reads, replicate 1.

#. ``cami_short``: CAMI-II Marine Illumina reads, replicate 1.

#. ``known_exp2``: mock microbial community with 50 known genomes,
   Illumina reads.

#. ``known_metaphlan_exp2``: mock microbial community with 25 known genomes,
   Illumina reads. Designed to compared TIPP3, TIPP3-fast, and MetaPhlAn4.

#. ``known_nanopore_exp2_2``: mock microbial community with 50 known genomes,
   Nanopore reads.

#. ``known_pacbio_exp2``: mock microbial community with 50 known genomes,
   PacBio reads.

#. ``mixed_exp2``: mock microbial community with 53 known genomes and 47
   novel genomes, Illumina reads.

#. ``mixed_metaphlan_exp2``: mock microbial community with 25 known genomes
   and 25 novel genomes, Illumina reads.
   Designed to compared TIPP3, TIPP3-fast, and MetaPhlAn4.

#. ``mixed_nanopore_exp2_2``: mock microbial community with 53 known genomes
   and 47 novel genomes, Nanopore reads.

#. ``mixed_pacbio_exp2``: mock microbial community with 53 known genomes and
   47 novel genomes, PacBio reads.

#. ``novel_exp2``: mock microbial community with 50 novel genomes,
   Illumina reads.

#. ``novel_metaphlan_exp2``: mock microbial community with 25 novel genomes,
   Illumina reads. Designed to compared TIPP3, TIPP3-fast, and MetaPhlAn4.

#. ``novel_nanopore_exp2_2``: mock microbial community with 50 novel genomes,
   Nanopore reads.

#. ``novel_pacbio_exp2``: mock microbial community with 50 novel genomes,
   PacBio reads.

Method IDs
++++++++++
#. ``tipp3``: TIPP3.

#. ``tipp3-fast``: TIPP3-fast.

#. ``kraken_all_reads``: Kraken2 with all reads as input.

#. ``kraken_filtered``: Kraken2 with filtered reads as input.

#. ``bracken_all_reads``: Bracken with all reads as input.

#. ``bracken_filtered``: Bracken with filtered reads as input.

#. ``metabuli_all_reads``: Metabuli with all reads as input.

#. ``metabuli_filtered``: Metabuli with filtered reads as input.

#. ``metaphlan_all_reads``: MetaPhlAn4 with all reads as input.

#. ``motus``: mOTUsv3 with all reads as input.

Computed Evaluation Metrics
~~~~~~~~~~~~~~~~~~~~~~~~~~~
All computed evaluation metrics, normalized Hellinger distance and Hellinger
distance, are included in the ``abundance_profile_combined.csv`` file. The
method code names are different from the ones in the ``PROFILES/``.

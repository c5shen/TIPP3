Taxonomic Identification and Phylogenetic Profiling (TIPP)
==========================================================
Latest version is TIPP2

TIPP(2) is a method for the following problems:

Taxonomic identification:
+ Input: A query sequence *q*
+ Output: The taxonomic lineage of *q*

Abundance profiling:
+ Input: A set *Q* of query sequences
+ Output: An abundance profile estimated on *Q*

TIPP is a modification of SEPP for classifying query sequences (i.e. reads) using phylogenetic placement. TIPP inserts each read into a taxonomic tree and uses the insertion location to identify the taxonomic lineage of the read. The novel idea behind TIPP is that rather than using the single best alignment and placement for taxonomic identification, we use a collection of alignments and placements and consider statistical support for each alignment and placement. Our study shows that TIPP provides improved classification accuracy on novel sequences and on sequences with evolutionarily divergent datasets. TIPP can also be used for abundance estimation by computing an abundance profile on the reads binned to marker genes in a reference dataset. TIPP2 provides an new reference dataset with 40 marker genes, assembled from the NCBI RefSeq database (learn more [here](https://github.com/shahnidhi/TIPP_reference_package)). In addition, TIPP2 updates how query sequences (i.e. reads) are mapped to marker genes. This repository corresponds to TIPP2, and henceforth we use the terms TIPP and TIPP2 interchangeably.

Developers of TIPP: TODO.

### Publications:
(TIPP3)

(TIPP2) Nguyen, Nam, Siavash Mirarab, Bo Liu, Mihai Pop, and Tandy Warnow, "TIPP: Taxonomic identification and phylogenetic profiling," *Bioinformatics*, 2014. [doi:10.1093/bioinformatics/btu721](http://bioinformatics.oxfordjournals.org/content/30/24/3548.full.pdf).

(TIPP) Shah, Nidhi, Erin K. Molloy, Mihai Pop, and Tandy Warnow, "TIPP2: metagenomic taxonomic profiling using phylogenetic markers," *Bioinformatics*, 2020. [doi:10.1093/bioinformatics/btab023](https://doi.org/10.1093/bioinformatics/btab023). Datasets used in the TIPP2 paper can be found [here](https://obj.umiacs.umd.edu/tipp/tipp-datasets.tar.gz). 

### Note and Acknowledgment: 
- TIPP uses the [Dendropy](http://pythonhosted.org/DendroPy/) package. 

-------------------------------------

Installing TIPP
===============
TODO

Requirements:
-------------
TODO

Installation Steps:
-------------------
TODO


---------------------------------------------

Running TIPP
============
TODO


---------------------------------------------

Bugs and Errors
===============
TODO

TIPP3 - Taxonomic Identification and Phylogenetic Profiling
==========================================================
|CHANGELOG|
  
:Developer:
    Chengze Shen

.. contents:: Table of contents
   :backlinks: top
   :local:

News
----
* Nothing new right now.

TODO list
---------
* 12.11.2024 - Working on a PyPI installation for TIPP3.

Method Overview
---------------
TIPP3 is a metagenomic profiling method that solves the following problems:

#. Taxonomic identification:
> Input: A query sequence *q*
> Output: The taxonomic lineage of *q* (if mapped)

#. Abundance profiling:
+ Input: A set *Q* of query sequences
+ Output: An abundance profile estimated on *Q*

TIPP3 continues the TIPP-family methods (prior methods: TIPP and TIPP2), which use a marker gene database to identify the taxonomic lineage of input reads (if the read comes from a marker gene). See the pipeline below for the TIPP3 workflow.

+-------------------------------------------+
| Publication                               |
+===========================================+
| (TIPP2) Nguyen, Nam, Siavash Mirarab,     |
| Bo Liu, Mihai Pop, and Tandy Warnow,      |
| "TIPP: Taxonomic identification and       |
| phylogenetic profiling."                  |
| Bioinformatics, 2014.                     |
| https://doi:10.1093/bioinformatics/btu721 |
+===========================================+
|                                           |
+-------------------------------------------+

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


.. |CHANGELOG| image:: https://img.shields.io/badge/CHANGELOG-gray?style=flat
   :alt: Static Badge
   :target: CHANGELOG

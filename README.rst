TIPP3 - Taxonomic Identification and Phylogenetic Profiling
===========================================================
|CHANGELOG|
  
:Developer:
    Chengze Shen

.. contents:: Table of contents
   :backlinks: top
   :local:

News
----
* 1.20.2025 - TIPP3 now is feature complete for abundance profiling, for both
  the more accurate TIPP3 mode or the fast TIPP3-fast mode. By default,
  TIPP3-fast is used.
* 12.17.2024 - TIPP3 now can handle both fasta (``.fasta`` or ``.fa``)
  and fastq (``.fastq`` or ``.fq``) files as inputs, and they can be in
  gzip format (e.g., ``queries.fasta.gz`` or ``queries.fq.gz``).

TODO list
---------
* 1.20.2025 - Add a script for downloading the latest reference package to a
  user specified directory.
* 12.11.2024 - Working on a PyPI installation for TIPP3.

Method Overview
---------------
TIPP3 is a metagenomic profiling method that solves the following problems:

Taxonomic identification
  - **Input**: A query read *q*
  - **Output**: The taxonomic lineage of *q* (if identified)

Abundance profiling
  - **Input**: A set *Q* of query reads
  - **Output**: An abundance profile estimated on *Q*

TIPP3 continues the TIPP-family methods (prior methods: TIPP and TIPP2),
which use a marker gene database to identify the taxonomic lineage of input
reads (if the read comes from a marker gene).
See the pipeline below for the TIPP3 workflow.

.. image:: https://chengzeshen.com/documents/tipp3/tipp3_overview.png
   :alt: TIPP3 pipeline
   :width: 100%
   :align: center

+------------------------------------------------+
| Publication(s)                                 |
+================================================+
| (TIPP3) Shen, Chengze, Eleanor Wedell,         |
| Mihai Pop, and Tandy Warnow, "TIPP3 and        |
| TIPP3-fast: improved abundance profiling in    |
| metagenomics." TBD.                            |
+------------------------------------------------+
| (TIPP2) Nguyen, Nam, Siavash Mirarab,          |
| Bo Liu, Mihai Pop, and Tandy Warnow,           |
| "TIPP: Taxonomic identification and            |
| phylogenetic profiling."                       |
| Bioinformatics, 2014.                          |
| https://doi.org/10.1093/bioinformatics/btu721  |
+------------------------------------------------+
| (TIPP) Shah, Nidhi, Erin K. Molloy, Mihai      |
| Pop, and Tandy Warnow,                         |
| "TIPP2: metagenomic taxonomic profiling        |
| using phylogenetic markers."                   |
| Bioinformatics, 2020.                          |
| https://doi.org/10.1093/bioinformatics/btab023 |
+------------------------------------------------+

Note and Acknowledgment 
~~~~~~~~~~~~~~~~~~~~~~~
TIPP3 includes and uses:

#. `pplacer <https://github.com/matsen/pplacer>`__ (v1.1.alpha19).
#. `Batch-SCAMPP <https://github.com/ewedell/BSCAMPP>`__ (v1.0.0).

External Requirements
---------------------
**BLAST** is a hard requirement to run TIPP3. The software will automatically
look for ``blastn`` in the ``$PATH`` environment variable.
If you have not installed BLAST, you can find the latest version from
`<https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`__. 

Installation
------------
TODO

Install with PyPI (``pip``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
TODO

Install from source files
~~~~~~~~~~~~~~~~~~~~~~~~~
TODO

Requirements
++++++++++++

::

   python>=3.7
   configparser>=5.0.0
   DendroPy>=4.5.2
   numpy>=1.21.6
   psutil>=5.0.0
   setuptools>=60.0.0
   treeswift>=1.1.28
   witch-msa>=1.0.4

Installation Steps
++++++++++++++++++

.. code:: bash

   # 1. Install via GitHub repo
   git clone https://github.com/c5shen/TIPP3.git

   # 2. Install all requirements
   pip3 install -r requirements.txt

   # 3. Execute tipp3.py executable for the first time with "-h" to see
   #    allowed commandline parameters and example usages
   #    Running TIPP3 for the first time will also create the main config
   #    file at "~/.tipp3/main.config", which stores the default behavior
   #    for running TIPP3 (including all binary executable paths)
   python3 tipp3.py [-h]    # OR tipp3-accurate [-h] OR tipp3-fast [-h]

``main.config``
~~~~~~~~~~~~~~~

``main.config`` file will be created the first time running TIPP3 at the user
root directory (``~/.tipp3/main.config``). This file stores the default
behavior for running TIPP3 and the paths to all binary executables that TIPP3
need to use.

user-specified config file
~~~~~~~~~~~~~~~~~~~~~~~~~~
In addition, a user can specify a customized config file with ``-c`` or
``--config-file`` parameter option when running TIPP3 (e.g.,
``tipp3.py -c user.config``). The ``user.config`` file will override settings
from ``main.config`` (if overlaps). Command-line arguments still have the
highest priority and will override both config files, if any parameters overlap.

Usage
-----
The general command to run TIPP3 is listed below. By default, TIPP3-fast is run,
which is significantly faster than the more accurate TIPP3 mode. See examples
below for how to customize the TIPP3 pipeline.

.. code:: bash

   # (Optional) change the logging level to DEBUG for more verbose logging
   export TIPP_LOGGING_LEVEL=debug

   python3 tipp3.py -r [reference package path] -i [query reads] -d [output directory]


Examples
~~~~~~~~

Some examples of TIPP3 usage can be found at the bottom of the help text
running:

.. code:: bash

   python3 tipp3.py -h


All of the following examples can be found in the **examples/run.sh** bash
script, with example data stored under **examples/data**. The default example
data used is a small set of Illumina short reads denoted as
``illumina.small.queries.fasta``.

Scenario 1
++++++++++
(TIPP3-fast) Use BLAST for query alignment, and Batch-SCAMPP (``bscampp``) for
query placement. Keep all temporary files during the run.

.. code:: bash

   python3 tipp3.py -i examples/illumina.small.queries.fasta \
      --reference-package [reference package dir] --outdir tipp3_scenario1 \
      --alignment-method blast --placement-method bscampp \
      -t 16 --keeptemp

Scenario 2
++++++++++
Use BLAST for query alignment, and pplacer with the taxtastic package for
query placement (``pplacer-taxtastic``). Keep all temporary files. 

.. code:: bash

   python3 tipp3.py -i examples/illumina.small.queries.fasta \
      --reference-package [reference package dir] --outdir tipp3_scenario1 \
      --alignment-method blast --placement-method pplacer-taxtastic \
      -t 16 --keeptemp

Scenario 3
++++++++++
(TIPP3) Use WITCH for query alignment, and ``pplacer-taxtastic`` for query
placement. Keep all temporary files.

.. code:: bash

   python3 tipp3.py -i examples/illumina.small.queries.fasta \
      --reference-package [reference package dir] --outdir tipp3_scenario1 \
      --alignment-method witch --placement-method pplacer-taxtastic \
      -t 16 --keeptemp


.. |CHANGELOG| image:: https://img.shields.io/badge/CHANGELOG-gray?style=flat
   :alt: Static Badge
   :target: CHANGELOG

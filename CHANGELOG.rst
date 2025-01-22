TIPP3 0.1b
----------
#. Included other minor bug fixes.
#. Changed the default file name from ``tipp3.py`` to ``run_tipp3.py`` to avoid
   conflict with versioning and installed packages.
#. Fixed installed binaries to make sure not conflicting with the actual
   ``tipp3`` packages. Now the installed binaries with PyPI are:
   ``run_tipp3.py`` (for customizing parameters),
   ``tipp3-accurate`` (for most accurate settings of TIPP3), and
   ``tipp3`` (for fastest settings of TIPP3).

TIPP3 0.1a
-----------
#. Working on an installation for PyPI, almost done.
#. Support ``.fasta, .fa, .fastq, .fq`` files as inputs. Also support them in gzipped format (e.g., ``.fasta.gz or .fasta.gzip``)
#. Lint-rolled all codes to fix unused variables and undefined variables.

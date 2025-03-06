TIPP3 v0.1
----------
#. Added the full pipeline to create a customized TIPP reference package. Please
   refer to the Wiki page and the Jupyter notebook in ``refpkg_scripts``.
#. Added subcommands: "abundance" and "download_refpkg".
   ``run_tipp3.py abundance`` has the original TIPP3 behavior for abundance
   profiling. ``run_tipp3.py download_refpkg`` can be used to download the
   latest TIPP3 reference package to a designated directory. For adjustment,
   The other two installed binaries ``tipp3`` and ``tipp3-fast`` now function
   only for the subcommand "abundance".
#. Fixed help text typos and adjusted exception logging to exit with
   return code 1.

TIPP3 v0.1b2
------------
#. Fixed a missing exit command for function ``tipp3_stop()``.
#. Minor help text fixes.

TIPP3 v0.1b1
------------
#. Fixed a bug in code that prevented BLASTN from reading in fasta/fa files
   correctly.
#. Added new output file ``query_classifications.tsv`` that aggregates all
   mapped query reads with their taxonomic identifications.
#. Included other minor bug fixes and code updates. 

TIPP3 v0.1b
----------
#. Included other minor bug fixes.
#. Changed the default file name from ``tipp3.py`` to ``run_tipp3.py`` to avoid
   conflict with versioning and installed packages.
#. Fixed installed binaries to make sure not conflicting with the actual
   ``tipp3`` packages. Now the installed binaries with PyPI are:
   ``run_tipp3.py`` (for customizing parameters),
   ``tipp3-accurate`` (for most accurate settings of TIPP3), and
   ``tipp3`` (for fastest settings of TIPP3).

TIPP3 v0.1a
-----------
#. Working on an installation for PyPI, almost done.
#. Support ``.fasta, .fa, .fastq, .fq`` files as inputs. Also support them in gzipped format (e.g., ``.fasta.gz or .fasta.gzip``)
#. Lint-rolled all codes to fix unused variables and undefined variables.

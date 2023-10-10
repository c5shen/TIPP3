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

Developers of TIPP: Nam Nguyen, Siavash Mirarab, Nidhi Shah, Erin Molloy, and Tandy Warnow.

### Publications:
Nguyen, Nam, Siavash Mirarab, Bo Liu, Mihai Pop, and Tandy Warnow, "TIPP: Taxonomic identification and phylogenetic profiling," *Bioinformatics*, 2014. [doi:10.1093/bioinformatics/btu721](http://bioinformatics.oxfordjournals.org/content/30/24/3548.full.pdf).

Shah, Nidhi, Erin K. Molloy, Mihai Pop, and Tandy Warnow, "TIPP2: metagenomic taxonomic profiling using phylogenetic markers," *Bioinformatics*, 2020. [doi:10.1093/bioinformatics/btab023](https://doi.org/10.1093/bioinformatics/btab023). Datasets used in the TIPP2 paper can be found [here](https://obj.umiacs.umd.edu/tipp/tipp-datasets.tar.gz). 

### Note and Acknowledgment: 
- TIPP uses the [Dendropy](http://pythonhosted.org/DendroPy/) package. 
- TIPP uses the [SEPP](https://github.com/smirarab/sepp/) package for alignment and placement steps.

-------------------------------------

Installing TIPP
===============
This section details steps for installing and running TIPP. We have run TIPP on Linux and MAC. If you experience difficulty installing or running the software, please contact one of us (Tandy Warnow or Siavash Mirarab).

Requirements:
-------------
Before installing the software you need to make sure the following programs are installed on your machine.

- Python: Version >= 3.7
- Java: Version >= 1.8
- Blast: Version > 2.10.1
- SEPP: Version >= 4.5.1

Installation Steps:
-------------------
TIPP requires SEPP to be installed. <br />
SEPP can be easily installed with conda - `conda install -c bioconda sepp` or you can install from source directly. Check instructions [here](https://github.com/smirarab/sepp/blob/master/README.SEPP.md). <br />
Once SEPP is installed, do the following


1. Download and decompress the reference dataset available at [https://obj.umiacs.umd.edu/tipp/tipp2-refpkg.tar.gz](https://obj.umiacs.umd.edu/tipp/tipp2-refpkg.tar.gz).
2. Set the environment variables (`REFERENCE` and `BLAST`) that will be used to create the configuration file. The `REFERENCE` environment variable should point to the location of the reference dataset (i.e. it should point to the `tipp2-refpkg` directory with its full path). The `BLAST` environment variable should to point to the location of the binary: `blastn`. Environment variables can be set using the following (shell-dependent) commands:
	- `export VARIABLE_NAME=/path/to/file` (bash shell)
	- `setenv VARIABLE_NAME /path/to/file` (tcsh shell)
3. Create the TIPP configuration file by running the command: `python setup.py tipp` or `python setup.py tipp -c` to avoid using the home directory. This  creates a `tipp.config` config file.
4. Install: run `python setup.py install`

Common Problems:
----------------
1. TIPP requires SEPP to be installed. If TIPP is not running, first check to see if SEPP was installed correctly by typing `run_sepp.py -h` to see it's user options.
2. If you are installing sepp with the `--user` flag, then you may need to set your `PYTHONPATH` environmental variable. To run TIPP, you may also need to include your python bin directory in your `PATH` environmental variable. If you are having issues installing or running sepp, try adding the following lines to your bash profile file `~/.bash_profile`:
```
export PYTHONPATH="$HOME/.local/lib/python3.7/site-packages:$PYTHONPATH"
export PATH="$HOME/.local/bin:$PATH"
```
where `3.7` is replaced with the version of python that you are running.

3. Several cluster systems that we have worked on have different commands to run different versions of python e.g. `python` runs python version 2.X and `python3` runs python version 3.Y. In this scenario, you need to be careful to be consistent with the python command (e.g. `python` or `python3`) when installing/configuring sepp and tipp.

4. TIPP relies on `blastn` for the binning of metagenomic reads, so BLAST needs to be downloaded and installed separately (learn more [here](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)). Then, point the `BLAST` environment variable to your installation of `blastn`. Alternatively, you can manually point TIPP to the `blastn` installation by modifying the `tipp.config` file. 

5. TIPP performs abundance profiling uses a set of 40 marker genes. This reference dataset needs to be downloaded separately from [here](https://obj.umiacs.umd.edu/tipp/tipp2-refpkg.tar.gz). Then, point the `REFERENCE` environment variable to the decompressed directory before installing TIPP. Alternatively, you can manually point TIPP to the reference dataset by modifying the `tipp.config` file. 

---------------------------------------------

Running TIPP
============
TIPP is a general pipeline for classifying reads belonging to a specific marker gene or complete set of marker genes.  We provide precomputed marker gene datasets for a collection of genes found in the tipp2-refpkg archive.  

The general command for running TIPP for a specific marker is:

`run_abundance.py -G <markers-v3> -f <input sequences> -d <output directory> --tempdir <intermediate results folder> --cpu 8 `

If you want to run abundance profile using only a select few genes, then provide gene names comma separated to the `-g` parameter. For example, 
`run_abundance.py -G <markers-v3> -g ArgS_COG0018,CysS_COG0215,Ffh_COG0541 -f <input sequences> -d <output directory> --tempdir <intermediate results folder> --cpu 8 `

To see options for running the script, use the command:

`run_abundance.py -h`

The main output of TIPP is a `_classification.txt` file that annotation for each read, and abundance.<taxonomic-level>.csv files with relative abundance at main taxonomic levels - phylum, class, order, family, genus, and species.

TIPP also outputs other information about the fragments, including the alignments (note that there could be multiple alignment files created, each corresponding to a different placement subset) as well as the phylogenetic placements.
Placements are given in the `.json` file, created according to *pplacer* format. Please refer to [pplacer website](http://matsen.github.com/pplacer/generated_rst/pplacer.html#json-format-specification) for more information on the form at of the `.json` file. Also note that *pplacer* package provides a program called *guppy* that can read `.json` files and  perform downstream steps such as visualization.

By setting `SEPP_DEBUG` environment variable to `True`, you can instruct SEPP to output more information that can be helpful for debugging.  

This [tutorial](tipp-tutorial.md) contains examples of running TIPP for read classification as well as abundance profiling using the script `run_abundance.py`.

---------------------------------------------

Bugs and Errors
===============
TIPP is under active research development by the Warnow Lab. Please report any errors to Tandy Warnow (warnow@illinois.edu), Siavash Mirarab (smirarab@ucsd.edu), Nidhi Shah (nidhi@umd.edu) and Erin Molloy (ekmolloy@cs.ucla.edu).



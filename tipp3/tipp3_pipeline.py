import gzip
import os
import sys
import tempfile
import re
#from sepp.alignment import MutableAlignment
#from sepp.alignment import _write_fasta
#from sepp.config import options
#from sepp import get_logger
from tipp3 import get_logger
from tipp3.read_binning import read_binning
from tipp3.read_alignment import read_alignment
from tipp3.read_placement import read_placement
from tipp3.utils import *
from tipp3.config import Config 

_LOG = get_logger(__name__)
global levels, character_map, taxon_map, level_map, key_map, refpkg
character_map = {'A': 'T', 'a': 't', 'C': 'G', 'c': 'g', 'T': 'A',
                 't': 'a', 'G': 'C', 'g': 'c', '-': '-'}
levels = ["species", "genus", "family", "order",
            "class", "phylum", "superkingdom"]


# Main pipeline of TIPP3
def tipp3_pipeline(*args, **kwargs):
    # (0) load refpkg


    # (1) read binning against the TIPP3 refpkg using BLAST


    # (2) read alignment to corresponding marker genes


    # (3) read placement to corresponding marker gene taxonomic trees


    # (4) collect results and abundance profile

'''
Load TIPP3 reference package in
'''
def load_reference_package():
    #TODO
    global refpkg

    refpkg = {}

    path = os.path.join(options().__getattribute__('reference').path,
                        options().genes)
    input = os.path.join(path, "file-map-for-tipp.txt")

    refpkg["genes"] = []
    with open(input) as f:
        for line in f.readlines():
            [key, val] = line.split('=')

            [key1, key2] = key.strip().split(':')
            val = os.path.join(path, val.strip())

            try:
                refpkg[key1][key2] = val
            except KeyError:
                refpkg[key1] = {}
                refpkg[key1][key2] = val

            if (key1 != "blast") and (key1 != "taxonomy"):
                refpkg["genes"].append(key1)

    refpkg["genes"] = set(refpkg["genes"])
    refpkg["genes"] = list(refpkg["genes"])

import time, os, sys, tempfile, re, gzip
from argparse import ArgumentParser, Namespace

from tipp3.configs import Configs, _root_dir, main_config_path, _read_config_file
from tipp3.configs import *
from tipp3 import get_logger, __version__
from tipp3.refpkg_loader import load_reference_package
from tipp3.read_binning import read_binning
from tipp3.read_alignment import read_alignment
from tipp3.read_placement import read_placement

from tipp3.helpers.general_tools import SmartHelpFormatter

_LOG = get_logger(__name__)
global levels, character_map, taxon_map, level_map, key_map
character_map = {'A': 'T', 'a': 't', 'C': 'G', 'c': 'g', 'T': 'A',
                 't': 'a', 'G': 'C', 'g': 'c', '-': '-'}
levels = ["species", "genus", "family", "order",
            "class", "phylum", "superkingdom"]

# Main pipeline of TIPP3
def tipp3_pipeline(*args, **kwargs):
    s1 = time.time()

    # (pre) parse arguments and build configurations
    parse_arguments()

    # (0) load refpkg
    refpkg = load_reference_package()

    # (1) read binning against the TIPP3 refpkg using BLAST
    read_binning(refpkg)

    # (2) read alignment to corresponding marker genes


    # (3) read placement to corresponding marker gene taxonomic trees


    # (4) collect results and abundance profile

    s2 = time.time()
    _LOG.INFO('TIPP3 completed in {} seconds...'.format(s2 - s1))
    print('TIPP3 completed in {} seconds...'.format(s2 - s1))

'''
parse argument and populate Configs
'''
def parse_arguments():
    global _root_dir, main_config_path
    parser = _init_parser()
    cmdline_args = sys.argv[1:]

    buildConfigs(parser, cmdline_args)
    getConfigs()

    _LOG.INFO('TIPP3 is running with: {}'.format(' '.join(cmdline_args)))
    if os.path.exists(main_config_path):
        _LOG.INFO('Main configuration loaded from {}'.format(
            main_config_path))

'''
initialize parser to read user inputs
'''
def _init_parser():
    parser = ArugmentParser(
            description=(
                "This program runs TIPP3, a taxonomic identification "
                "and abundance profiling tool for metagenomic reads. ",
            conflict_handler='resolve',
            formatter_class=SmartHelpFormatter)
    parser.add_argument('-v', '--version', action='version',
        version="%(prog)s " + __version__)

    # basic settings
    basic_group = parser.add_argument_group(
            "Basic parameters".upper(),
            ("These are basic fields for running TIPP3. "
             "Users need to provide the path to a TIPP3-compatible refpkg "
             "and the path to the query reads they wish to profile."))
    parser.groups = dict()
    parser.groups['basic_group'] = basic_group
#TODO

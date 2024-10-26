import time, os, sys, tempfile, re, gzip
from argparse import ArgumentParser, Namespace, RawDescriptionHelpFormatter

from tipp3.configs import Configs, _root_dir, main_config_path, _read_config_file
from tipp3.configs import *
from tipp3 import get_logger, __version__
from tipp3.refpkg_loader import loadReferencePackage
from tipp3.query_binning import queryBinning 
from tipp3.query_alignment import queryAlignment
from tipp3.query_placement import queryPlacement

#from tipp3.helpers.general_tools import SmartHelpFormatter

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
    parseArguments()

    # (0) load refpkg
    refpkg = loadReferencePackage(Configs.refpkg_path, Configs.refpkg_version)

    # (1) read binning against the TIPP3 refpkg using BLAST
    query_paths, query_alignment_paths = queryBinning(refpkg)

    # (2) read alignment to corresponding marker genes
    if Configs.alignment_method != 'blast':
        query_alignment_paths = queryAlignment(refpkg, query_paths)

    # (3) read placement to corresponding marker gene taxonomic trees


    # (4) collect results and abundance profile

    s2 = time.time()
    _LOG.info('TIPP3 completed in {} seconds...'.format(s2 - s1))
    print('TIPP3 completed in {} seconds...'.format(s2 - s1))

'''
parse argument and populate Configs
'''
def parseArguments():
    global _root_dir, main_config_path
    parser = _init_parser()
    cmdline_args = sys.argv[1:]

    buildConfigs(parser, cmdline_args)
    #_LOG = get_logger(__name__, log_path=Configs.log_path)

    _LOG.info('TIPP3 is running with: {}'.format(' '.join(cmdline_args)))
    #if os.path.exists(main_config_path):
    #    _LOG.info('Main configuration loaded from {}'.format(
    #        main_config_path))
    getConfigs()

'''
initialize parser to read user inputs
'''
def _init_parser():
    # example usages
    example_usages = '''Example usages:
> TIPP3 default behavior
    %(prog)s -r refpkg_dir/ -i queries.fasta
> only output read alignment to marker genes (then exit) 
    %(prog)s -r refpkg_dir/ -i queries.fasta --alignment-only
'''

    parser = ArgumentParser(
            description=(
                "This program runs TIPP3, a taxonomic identification "
                "and abundance profiling tool for metagenomic reads. "),
            conflict_handler='resolve',
            epilog=example_usages,
            formatter_class=RawDescriptionHelpFormatter)
            #formatter_class=SmartHelpFormatter)
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
    basic_group.add_argument('-i', '--query-path', type=str,
        help='Path to a set of unaligned query reads for classification',
        required=True)
    #TODO - not making this required and invoke download for TIPP3 refpkg
    # if no refpkg is supplied
    basic_group.add_argument('-r', '--refpkg-path', '--refpkg',
            '--reference-package',
        type=str, help='Path to a TIPP3-compatible refpkg.',
        required=False, default=None)
    basic_group.add_argument('--refpkg-version',
        type=str, help='Version of the refpkg. [default: markers-v4]',
        default='markers-v4', required=False)
    basic_group.add_argument('--alignment-method',
        type=str, choices=['witch', 'blast', 'hmm'], default='witch',
        help=' '.join(['Alignment method to use for aligning reads',
            'to marker genes. [default: witch]']),
        required=False)
    basic_group.add_argument('--placement-method',
        type=str, choices=['pplacer', 'bscampp'],
        help=' '.join(['Placement method to use for placing aligned reads',
            'to marker gene taxonomic trees. [default: pplacer']),
        default='pplacer', required=False)
    basic_group.add_argument('-c', '--config-file',
        type=str, help='Path to a user-defined config file. [default: None]',
        required=False, default=None)
    basic_group.add_argument('-d', '--outdir',
        type=str, help='Path to the desired output directory [default: ./tipp3_output]',
        required=False, default='./tipp3_output')
    basic_group.add_argument('-t', '--num-cpus',
        type=int, help='Number of CPUs for multi-processing. [default: -1 (all)]',
        required=False, default=-1)

    # miscellaneous group
    misc_group = parser.add_argument_group(
            "Miscellaneous options".upper(),
            ("Optional parameters for TIPP3 setup/config etc."))
    parser.groups['misc_group'] = misc_group
    misc_group.add_argument('--keeptemp', action='store_const', const=True,
            help='Keep temporary files in the running process.',
            default=False)
    misc_group.add_argument('-y', '--bypass-setup', action='store_const',
            const=True, default=True,
            help=' '.join(['(Optional) Include this argument to bypass',
                'the initial step when running TIPP3 to set up the',
                'configuration directory (will use ~/.tipp3).',
                'Note: By default this option is one.',
                'You also only need to run this option once.']))

    return parser

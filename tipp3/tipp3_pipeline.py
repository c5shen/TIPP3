import time, os, sys, shutil
from argparse import ArgumentParser, Namespace, RawDescriptionHelpFormatter

from tipp3 import get_logger, __version__
from tipp3.configs import Configs, _root_dir, main_config_path, _read_config_file
from tipp3.configs import *
from tipp3.refpkg_loader import loadReferencePackage
from tipp3.query_binning import queryBinning 
from tipp3.query_alignment import queryAlignment
from tipp3.query_placement import queryPlacement
from tipp3.query_abundance import queryAbundance 
from tipp3.helpers.general_tools import tqdm_styles 

from multiprocessing import Lock, Manager
from concurrent.futures import ProcessPoolExecutor

_LOG = get_logger(__name__)

global levels, character_map, taxon_map, level_map, key_map
character_map = {'A': 'T', 'a': 't', 'C': 'G', 'c': 'g', 'T': 'A',
                 't': 'a', 'G': 'C', 'g': 'c', '-': '-'}
levels = ["species", "genus", "family", "order",
            "class", "phylum", "superkingdom"]

# Main pipeline of TIPP3
def tipp3_pipeline(*args, **kwargs):
    s1 = time.time()
    m = Manager()
    lock = m.Lock()

    # (pre) parse arguments and build configurations
    parser, cmdline_args = parseArguments()

    # initialize ProcessPoolExecutor
    _LOG.warning('Initializing ProcessorPoolExecutor instance...')
    pool = ProcessPoolExecutor(Configs.num_cpus,
            initializer=initiate_pool, initargs=(parser, cmdline_args,))

    # (0) load refpkg
    refpkg = loadReferencePackage(Configs.refpkg_path, Configs.refpkg_version)

    # (1) read binning against the TIPP3 refpkg using BLAST
    query_paths, query_alignment_paths = queryBinning(refpkg)

    # (2) read alignment to corresponding marker genes
    if Configs.alignment_method != 'blast':
        query_alignment_paths = queryAlignment(refpkg, query_paths)

    # early stop --> alignment-only 
    if Configs.alignment_only:
        _LOG.warning("User specifies to output query alignment to marker "
                "genes only. Stopping TIPP3 now.")
        _LOG.warning("You can find the alignment files at: "
                f"{os.path.join(Configs.outdir, 'query_alignments')}")
        tipp3_stop(s1)

    # (3) read placement to corresponding marker gene taxonomic trees
    query_placement_paths = queryPlacement(refpkg, query_alignment_paths)

    # (4) collect results and abundance profile
    queryAbundance(refpkg, query_placement_paths, pool, lock)

    # close ProcessPoolExecutor
    _LOG.warning('Closing ProcessPoolExecutor instance...')
    pool.shutdown()
    _LOG.warning('ProcessPoolExecutor instance closed.')

    # cleaning up temporary files
    if not Configs.keeptemp:
        _LOG.info("Removing intermediate output files...")
        tipp3_clean_temp()

    # stop TIPP3
    tipp3_stop(s1)

def tipp3_clean_temp():
    temp_dirs = ['blast_output', 'query', 'query_alignments',
            'query_placements']#, 'query_classification']
    for temp in temp_dirs:
        shutil.rmtree(os.path.join(Configs.outdir, temp))
        _LOG.info(f"Removed {temp}")

def tipp3_stop(start_time):
    s2 = time.time()
    _LOG.info('TIPP3 completed in {} seconds...'.format(s2 - start_time))
    print('TIPP3 completed in {} seconds...'.format(s2 - start_time))

'''
Init function for a queue and get configurations for each worker
'''
def initiate_pool(parser, cmdline_args):
    buildConfigs(parser, cmdline_args, child_process=True)

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
    return parser, cmdline_args

'''
initialize parser to read user inputs
'''
def _init_parser():
    # example usages
    example_usages = '''Example usages:
> TIPP3 default behavior
    %(prog)s -r refpkg_dir/ -i queries.fasta
> Only output read alignment to marker genes (then exit) 
    %(prog)s -r refpkg_dir/ -i queries.fasta --alignment-only
> Running TIPP3-fast
    %(prog)s -r refpkg_dir/ -i queries.fasta --alignment-method blast --placement-method bscampp
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
        help=' '.join(['Path to a set of unaligned query reads',
            'for classification.', 'Accepted format:'
            '.fa/.fasta (can be compressed as a .gz file).']),
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
        type=str, choices=['pplacer-taxtastic', 'bscampp'],
        help=' '.join(['Placement method to use for placing aligned reads',
            'to marker gene taxonomic trees. [default: pplacer-taxtastic']),
        default='pplacer-taxtastic', required=False)
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
    #misc_group.add_argument('--verbose', type=int,
    #        help=' '.join(["Verbose level for logging.",
    #        "0: error, 1: info, >1: debug. [default: 1]"]),
    #        default=1, required=False)
    misc_group.add_argument('--alignment-only', action='store_const',
            const=True, default=False,
            help='Only obtain query alignments to marker genes and stop TIPP3.')
    misc_group.add_argument('--keeptemp', action='store_const', const=True,
            help='Keep temporary files in the running process.',
            default=False)
    misc_group.add_argument('-y', '--bypass-setup', action='store_const',
            const=True, default=True,
            help=' '.join(['(Optional) Include this argument to bypass',
                'the initial step when running TIPP3 to set up the',
                'configuration directory (will use ~/.tipp3).',
                'Note: By default this option is enabled.']))

    return parser

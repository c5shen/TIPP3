import time, os, sys, shutil
from argparse import ArgumentParser, Namespace, RawDescriptionHelpFormatter

from tipp3 import get_logger, __version__
from tipp3.configs import Configs, _root_dir, main_config_path, _read_config_file
from tipp3.configs import *
from tipp3.refpkg_loader import loadReferencePackage, downloadReferencePackage
from tipp3.query_binning import queryBinning 
from tipp3.query_alignment import queryAlignment
from tipp3.query_placement import queryPlacement
from tipp3.query_abundance import queryAbundance 
from tipp3.helpers.general_tools import *

from multiprocessing import Lock, Manager
from concurrent.futures import ProcessPoolExecutor

_LOG = get_logger(__name__)

global levels, character_map, taxon_map, level_map, key_map
character_map = {'A': 'T', 'a': 't', 'C': 'G', 'c': 'g', 'T': 'A',
                 't': 'a', 'G': 'C', 'g': 'c', '-': '-'}
levels = ["species", "genus", "family", "order",
            "class", "phylum", "superkingdom"]
# available subcommands. If nothing is provided, will default to 'abundance' 
subcommands = ['abundance', 'download_refpkg']

'''
Preset mode to run TIPP3 for abundance profiling
'''
def run_tipp3():
    tipp3_pipeline(mode='tipp3', subcommand='abundance')

'''
Preset mode to run TIPP3-fast for abundance profiling
'''
def run_tipp3_fast():
    tipp3_pipeline(mode='tipp3-fast', subcommand='abundance')

# Main pipeline of TIPP3
def tipp3_pipeline(*args, **kwargs):
    s1 = time.time()
    m = Manager()
    lock = m.Lock()

    # (pre) parse arguments and build configurations
    parser, cmdline_args = parseArguments(
            mode=kwargs.get('mode', None),
            subcommand=kwargs.get('subcommand', None))

    if Configs.command == 'download_refpkg':
        downloadReferencePackage(Configs.outdir, Configs.decompress)
    elif Configs.command == 'abundance':
        # initialize ProcessPoolExecutor
        _LOG.warning('Initializing ProcessorPoolExecutor instance...')
        pool = ProcessPoolExecutor(Configs.num_cpus,
                initializer=initiate_pool, initargs=(parser, cmdline_args,))
        
        # create output directory
        if not os.path.exists(Configs.outdir):
            os.makedirs(Configs.outdir)

        # (0) load refpkg
        refpkg = loadReferencePackage(Configs.refpkg_path, Configs.refpkg_version)

        # (1) read binning against the TIPP3 refpkg using BLAST
        query_paths, query_alignment_paths = queryBinning(refpkg)
        s2 = time.time()
        _LOG.info(f"Runtime for mapping reads to marker genes (seconds): {s2 - s1}") 

        # (2) read alignment to corresponding marker genes
        if Configs.alignment_method != 'blast':
            query_alignment_paths = queryAlignment(refpkg, query_paths)
        s3 = time.time()
        _LOG.info(f"Runtime for aligning reads to marker genes (seconds): {s3 - s2}") 

        # early stop --> alignment-only 
        if Configs.alignment_only:
            _LOG.warning("User specifies to output query alignment to marker "
                    "genes only. Stopping TIPP3 now.")
            _LOG.warning("You can find the alignment files at: "
                    f"{os.path.join(Configs.outdir, 'query_alignments')}")
            tipp3_stop(s1)

        # (3) read placement to corresponding marker gene taxonomic trees
        query_placement_paths = queryPlacement(refpkg, query_alignment_paths)
        s4 = time.time()
        _LOG.info(f"Runtime for placing reads to marker gene taxonomies (seconds): {s4 - s3}") 

        # (4) collect results and abundance profile
        queryAbundance(refpkg, query_placement_paths, pool, lock)
        s5 = time.time()
        _LOG.info(f"Runtime for obtaining abundance profile (seconds): {s5 - s4}") 

        # close ProcessPoolExecutor
        _LOG.warning('Closing ProcessPoolExecutor instance...')
        pool.shutdown()
        _LOG.warning('ProcessPoolExecutor instance closed.')

        # cleaning up temporary files
        if not Configs.keeptemp:
            _LOG.info("Removing intermediate output files...")
            tipp3_clean_temp()
            s6 = time.time()
            _LOG.info(f"Runtime for cleaning temporary files (seconds): {s6 - s5}") 

    # stop TIPP3
    tipp3_stop(s1)

def tipp3_clean_temp():
    temp_dirs = ['blast_output', 'query', 'query_alignments',
            'query_placements']#, 'query_classifications']
    for temp in temp_dirs:
        shutil.rmtree(os.path.join(Configs.outdir, temp))
        _LOG.info(f"Removed {temp}")

def tipp3_stop(start_time):
    send = time.time()
    _LOG.info('TIPP3 completed in {} seconds...'.format(send - start_time))
    #print('TIPP3 completed in {} seconds...'.format(s2 - start_time))
    exit(0)

'''
Init function for a queue and get configurations for each worker
'''
def initiate_pool(parser, cmdline_args):
    buildConfigs(parser, cmdline_args, child_process=True)

'''
parse argument and populate Configs
'''
def parseArguments(mode=None, subcommand=None):
    global _root_dir, main_config_path
    parser = _init_parser(mode)

    cmdline_args = sys.argv[1:]
    # check if subcommand is provided, if not try to append a subcommand
    # if the default is given (e.g., running the binary `tipp3` or
    # `tipp3-accurate`
    if cmdline_args[0] not in subcommands:
        # check if it is an argument or an invalid subcommand
        if cmdline_args[0].startswith('-') and subcommand is not None:
            cmdline_args = [subcommand] + cmdline_args

    buildConfigs(parser, cmdline_args)
    #_LOG = get_logger(__name__, log_path=Configs.log_path)

    getConfigs(arguments=cmdline_args)
    _LOG.info('TIPP3 is running with: {}'.format(' '.join(cmdline_args)))
    #if os.path.exists(main_config_path):
    #    _LOG.info('Main configuration loaded from {}'.format(
    #        main_config_path))
    return parser, cmdline_args

'''
initialize parser to read user inputs
'''
def _init_parser(mode=None):
    # example usages
    example_usages = '''Example usages:
> abundance: run default for profiling (TIPP3-fast)
    %(prog)s abundance -r refpkg_dir/ -i queries.fasta[.gz]
    %(prog)s abundance -r refpkg_dir/ -i queries.fq[.gz]
> abundance: Only output read alignment to marker genes (then exit) 
    %(prog)s abundance -r refpkg_dir/ -i queries.fasta[.gz] --alignment-only
> download_refpkg: Download the latest TIPP3 reference package to current directory and decompress
    %(prog)s download_refpkg -d ./ --decompress
'''

    # determine which mode we have by default (default to tipp3-fast)
    _mode = 'tipp3-fast'
    if mode is not None:
        _mode = mode

    parser = ArgumentParser(
            description=(
                "This program runs TIPP3, a taxonomic identification "
                "and abundance profiling tool for metagenomic reads."),
            conflict_handler='resolve',
            epilog=example_usages,
            #formatter_class=RawDescriptionHelpFormatter)
            formatter_class=SmartHelpFormatter)
    parser.add_argument('-v', '--version', action='version',
        version="%(prog)s " + __version__)

    # add sub-commands
    subparsers = parser.add_subparsers(dest='command', help=None)

    # (1) DEFAULT abundance -- abundance profiling
    subparser_abs = subparsers.add_parser('abundance',
            help="(Default) Abundance profiling on input reads.")

    # (2) download_refpkg -- download reference package to target directory
    subparser_refpkg = subparsers.add_parser('download_refpkg',
            help="Download the latest TIPP3 reference package.")

    # (3) TBD

########################## subcommand: abundance ##############################
    # basic settings
    basic_group = subparser_abs.add_argument_group(
            "Basic parameters".upper(),
            ("These are basic fields for running TIPP3. "
             "Users need to provide the path to a TIPP3-compatible refpkg "
             "and the path to the query reads they wish to profile."))
    subparser_abs.groups = dict()
    subparser_abs.groups['basic_group'] = basic_group
    basic_group.add_argument('-i', '--query-path', type=str,
        help=' '.join(['Path to a set of unaligned query reads',
            'for classification.', 'Accepted format:'
            '.fa/.fasta/.fq/.fastq (can be compressed as a .gz file).']),
        required=True)
    #TODO - not making this required and invoke download for TIPP3 refpkg
    # if no refpkg is supplied
    basic_group.add_argument('-r', '--refpkg-path', '--refpkg',
            '--reference-package',
        type=str, help=' '.join(['Path to a TIPP3-compatible refpkg.',
            'Use subcommand \'download_refpkg\' to download the latest',
            'TIPP3 reference package.']),
        required=False, default=None)
    basic_group.add_argument('--refpkg-version',
        type=str, help='Version of the refpkg. [default: markers-v4]',
        default='markers-v4', required=False)
    basic_group.add_argument('--mode',
        type=str, choices=['tipp3', 'tipp3-fast'], default=mode,
        help=' '.join(['Preset mode for running TIPP3.', f'[default: {_mode}]',
            '\n\"tipp3\": the most accurate setting, with WITCH alignment',
            'and pplacer placement.',
            '\n\"tipp3-fast\": the fastest setting, with BLAST alignment',
            'and Batch-SCAMPP placement.',
            '\nThe mode will be overridden by parameters --alignment-method',
            'and --placement-method.']), required=False)
    basic_group.add_argument('--alignment-method',
        type=str, choices=['witch', 'blast', 'hmm'], default=None,
        help=' '.join(['Alignment method to use for aligning reads',
            'to marker genes. [default: using --mode]']),
        required=False)
    basic_group.add_argument('--placement-method',
        type=str, choices=['pplacer-taxtastic', 'bscampp'], default=None,
        help=' '.join(['Placement method to use for placing aligned reads',
            'to marker gene taxonomic trees. [default: using --mode]']),
        required=False)
    basic_group.add_argument('-c', '--config-file',
        type=str, help=' '.join(['Path to a user-defined config file.',
            'see an example at\n',
            'https://github.com/c5shen/TIPP3/blob/main/custom.config',
            '[default: None]']),
        required=False, default=None)
    basic_group.add_argument('-d', '--outdir',
        type=str, help='Path to the desired output directory [default: ./tipp3_output]',
        required=False, default='./tipp3_output')
    basic_group.add_argument('-t', '--num-cpus',
        type=int, help='Number of CPUs for multi-processing. [default: -1 (all)]',
        required=False, default=-1)

    # miscellaneous group
    misc_group = subparser_abs.add_argument_group(
            "Miscellaneous options".upper(),
            ("Optional parameters for TIPP3 setup/config etc."))
    subparser_abs.groups['misc_group'] = misc_group
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

######################### subcommand: download_refpkg #########################
    refpkg_basic_group = subparser_refpkg.add_argument_group(
            "Basic parameters".upper(),
            ("Use this subcommand to download the latest TIPP3 reference "
             "package, and optionally decompress it for usage."))
    refpkg_basic_group.add_argument('-d', '--outdir',
            type=str, help=' '.join(['Path to put the downloaded refpkg file.',
                '[default: ./]']), default='./')
    refpkg_basic_group.add_argument('--decompress', action='store_const',
            const=True, help='After download, decompress the file for usage.',
            default=False)

    return parser

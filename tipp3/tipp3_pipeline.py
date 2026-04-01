"""
TIPP3 main pipeline: orchestrates binning, alignment, placement, and profiling.
"""

import time, os, sys, shutil
from argparse import ArgumentParser
from multiprocessing import Manager
from concurrent.futures import ProcessPoolExecutor

from tipp3 import get_logger, __version__
from tipp3.configs import (
    Configs, _root_dir, main_config_path, _read_config_file,
    buildConfigs, getConfigs,
)
from tipp3.refpkg_loader import loadReferencePackage, downloadReferencePackage
from tipp3.query_binning import queryBinning
from tipp3.query_alignment import queryAlignment
from tipp3.query_placement import queryPlacement
from tipp3.query_abundance import (
    getAllClassification, getAbundanceProfile, getSpeciesDetection,
)
from tipp3.helpers.general_tools import SmartHelpFormatter

_LOG = get_logger(__name__)

SUBCOMMANDS = ['abundance', 'download_refpkg', 'detection']

def run_tipp3():
    """Entry point for `tipp3-accurate` (WITCH + pplacer)."""
    tipp3_pipeline(mode='tipp3', subcommand='abundance')


def run_tipp3_fast():
    """Entry point for `tipp3` (BLAST + BSCAMPP)."""
    tipp3_pipeline(mode='tipp3-fast', subcommand='abundance')


def tipp3_pipeline(*args, **kwargs):
    """Main TIPP3 pipeline: binning -> alignment -> placement -> profiling."""
    s1 = time.time()

    parser, cmdline_args = parseArguments(
        mode=kwargs.get('mode', None),
        subcommand=kwargs.get('subcommand', None))

    if Configs.command == 'download_refpkg':
        downloadReferencePackage(Configs.outdir, Configs.decompress)
        _tipp3_finish(s1)
        return

    if Configs.command not in ('abundance', 'detection'):
        _LOG.error(f"Unknown command: {Configs.command}")
        sys.exit(1)

    m = Manager()
    lock = m.Lock()
    _LOG.info('Initializing ProcessPoolExecutor...')
    pool = ProcessPoolExecutor(
        Configs.num_cpus,
        initializer=initiate_pool, initargs=(parser, cmdline_args))

    try:
        os.makedirs(Configs.outdir, exist_ok=True)

        # (0) Load reference package
        refpkg = loadReferencePackage(Configs.refpkg_path, Configs.refpkg_version)

        # (1) Read binning via BLAST
        query_paths, query_alignment_paths = queryBinning(refpkg, Configs.query_path)
        s2 = time.time()
        _LOG.info(f"Runtime for query read binning (seconds): {s2 - s1}")

        # (2) Read alignment
        if Configs.alignment_method != 'blast':
            query_alignment_paths = queryAlignment(refpkg, query_paths)
        s3 = time.time()
        _LOG.info(f"Runtime for query read alignment (seconds): {s3 - s2}")

        # Early stop for alignment-only mode
        if Configs.alignment_only:
            _LOG.info("Alignment-only mode: stopping after alignment. "
                      f"Output: {os.path.join(Configs.outdir, 'query_alignments')}")
            _tipp3_finish(s1)
            return

        # (3) Read placement
        query_placement_paths = queryPlacement(refpkg, query_alignment_paths)
        s4 = time.time()
        _LOG.info(f"Runtime for query read placement (seconds): {s4 - s3}")

        # (4) Classification and profiling
        classification_paths, filtered_paths = getAllClassification(
            refpkg, query_placement_paths, pool, lock)

        if Configs.command == 'abundance':
            getAbundanceProfile(refpkg, filtered_paths)
            s5 = time.time()
            _LOG.info(f"Runtime for abundance profiling (seconds): {s5 - s4}")

        elif Configs.command == 'detection':
            detection_thresholds = {'conservative': 0.2, 'sensitive': 0.12}
            if Configs.detection_threshold is not None:
                dt = Configs.detection_threshold
                if 0.0 <= dt <= 1.0:
                    detection_thresholds['custom'] = dt
                else:
                    _LOG.warning(
                        f"Detection threshold {dt} out of [0, 1] range, ignored.")

            getSpeciesDetection(detection_thresholds, refpkg, classification_paths)
            s5 = time.time()
            _LOG.info(f"Runtime for species detection (seconds): {s5 - s4}")

    finally:
        _LOG.info('Shutting down ProcessPoolExecutor...')
        pool.shutdown(wait=True)

    if not Configs.keeptemp:
        _LOG.info("Cleaning up intermediate files...")
        _tipp3_clean_temp()

    _tipp3_finish(s1)


def _tipp3_clean_temp():
    """Remove intermediate output directories."""
    temp_dirs = ['blast_output', 'query', 'query_alignments', 'query_placements']
    for temp in temp_dirs:
        path = os.path.join(Configs.outdir, temp)
        if os.path.isdir(path):
            shutil.rmtree(path)
            _LOG.debug(f"Removed {path}")
        else:
            _LOG.debug(f"Skipping {path} (not found)")


def _tipp3_finish(start_time):
    """Log completion time."""
    elapsed = time.time() - start_time
    _LOG.info(f'TIPP3 completed in {elapsed} seconds.')

def initiate_pool(parser, cmdline_args):
    """Initialize each worker process with the current configuration."""
    buildConfigs(parser, cmdline_args, child_process=True)


def parseArguments(mode=None, subcommand=None):
    """Parse command-line arguments and build Configs."""
    parser = _init_parser(mode)
    cmdline_args = sys.argv[1:]

    if not cmdline_args:
        parser.print_help()
        sys.exit(1)

    if cmdline_args[0] not in SUBCOMMANDS:
        if cmdline_args[0].startswith('-') and subcommand is not None:
            cmdline_args = [subcommand] + cmdline_args

    buildConfigs(parser, cmdline_args)
    getConfigs(arguments=cmdline_args)
    _LOG.info(f'TIPP3 v{__version__} running: {" ".join(cmdline_args)}')
    return parser, cmdline_args

def _init_parser(mode=None):
    """Build the argument parser with all subcommands and options."""
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
            formatter_class=SmartHelpFormatter)
    parser.add_argument('-v', '--version', action='version',
        version="%(prog)s " + __version__)

    # add sub-commands
    subparsers = parser.add_subparsers(dest='command', help=None)

    # (1) DEFAULT abundance -- abundance profiling
    subparser_abs = subparsers.add_parser('abundance',
        help="(Default) Abundance profiling on input reads.",
        formatter_class=SmartHelpFormatter)

    # (2) download_refpkg -- download reference package to target directory
    subparser_refpkg = subparsers.add_parser('download_refpkg',
        help="Download the latest TIPP3 reference package.",
        formatter_class=SmartHelpFormatter)

    # (3) species detection -- detecting if species is present or not
    subparser_detection = subparsers.add_parser('detection',
        help="Species detection on input reads.",
        formatter_class=SmartHelpFormatter)

#################### subcommand: abundance/detection ##########################
    # abundance and detection share the same groups of arguments
    for parser_name, _subparser in {
            'abundance': subparser_abs,
            'detection': subparser_detection
            }.items():
        _subparser.groups = dict()
        # basic settings
        basic_group = _subparser.add_argument_group(
                "Basic parameters".upper(),
                ("These are basic fields for running TIPP3. "
                 "Users need to provide the path to a TIPP3-compatible refpkg "
                 "and the path to the query reads they wish to profile."))
        _subparser.groups['basic_group'] = basic_group

        basic_group.add_argument('-i', '--query-path', type=str,
            help=' '.join(['Path to a set of unaligned query reads',
                'for classification.', 'Accepted format:'
                '.fa/.fasta/.fq/.fastq (can be compressed as a .gz file).']),
            required=True)
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
            type=str, choices=['tipp3', 'tipp3-fast'], default=_mode,
            help=' '.join(['Preset mode for running TIPP3.', f'[default: {_mode}]',
                '\n\"tipp3\": the most accurate setting, with WITCH alignment',
                'and pplacer placement.',
                '\n\"tipp3-fast\": the fastest setting, with BLAST alignment',
                'and Batch-SCAMPP placement.',
                'The mode will be overridden by parameters --alignment-method',
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
                'see an example at',
                'https://github.com/c5shen/TIPP3/blob/main/custom.config',
                '[default: None]']),
            required=False, default=None)
        basic_group.add_argument('-d', '--outdir',
            type=str, help='Path to the desired output directory [default: ./tipp3_output]',
            required=False, default='./tipp3_output')
        basic_group.add_argument('-t', '--num-cpus',
            type=int, help='Number of CPUs for multi-processing. [default: -1 (all)]',
            required=False, default=-1)
        # detection specific parameter
        if parser_name == 'detection': 
            basic_group.add_argument('-B', '--detection-threshold',
                type=float, help=' '.join(['(Detection only) Customized threshold for',
                    'species detection, other than the conservative (0.2)',
                    'and sensitive (0.12) values, ranging from [0, 1].',
                    'Will write detected species to',
                    'a separate file named \"detected_species_custom.tsv\".',
                    '[default: None]']),
                required=False, default=None)


        # miscellaneous group
        misc_group = _subparser.add_argument_group(
                "Miscellaneous options".upper(),
                ("Additional parameters for TIPP3 setup/config etc."))
        _subparser.groups['misc_group'] = misc_group
        misc_group.add_argument('--bscampp-mode', type=str,
            choices=['pplacer', 'epa-ng'], default='pplacer',
            help=' '.join(['Base placement method to use in BSCAMPP,',
                'currently supporting pplacer and epa-ng.',
                'Has priority and will override settings in',
                'main.config or a customized config file.'
                '[default: pplacer]'])) 
        misc_group.add_argument('--alignment-only', action='store_const',
            const=True, default=False,
            help='Only obtain query alignments to marker genes and stop TIPP3.')
        misc_group.add_argument('--keeptemp', action='store_const', const=True,
            help='Keep temporary files in the running process.',
            default=False)
        misc_group.add_argument('-y', '--bypass-setup', action='store_const',
            const=True, default=True,
            help=' '.join(['(DEPRECATED) Include this argument to bypass',
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

import os, time
try:
    import configparser
except ImportError:
    from ConfigParser import configparser
from argparse import ArgumentParser, Namespace
from platform import platform
from tipp3.init_configs import init_config_file
from tipp3 import get_logger

# detect home.path location or creating one if it is missing
homepath = os.path.dirname(__file__) + '/home.path'
_root_dir, main_config_path = init_config_file(homepath)

# set of valid configparse section names
valid_config_sections = ['witch', 'bscampp', 'pplacer-taxtastic', 
        'blast', 'refpkg'] 
logging_levels = set(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'])

_LOG = get_logger(__name__)

# default tqdm style for progress bar
tqdm_styles = {
        'desc': '\tRunning...', 'ascii': False,
        'ncols': 80, 
        #'disable': True,
        #'colour': 'green',
        'mininterval': 0.5
        }

'''
Configurations defined by users
'''
class Configs:
    global _root_dir

    # subcommand and verbose
    command = None
    verbose = 'INFO'

    # Multiprocessing settings
    num_cpus = -1
    max_concurrent_jobs = None

    # basic input items
    query_path = None
    refpkg_path = None     # e.g., xxx/yyy/tipp3-refpkg
    outdir = None
    config_file = None     # Added @ 7.25.2024 

    # choices of parameters
    # default to TIPP3-fast
    mode = 'tipp3-fast'
    alignment_method = 'blast'  # or blast
    placement_method = 'bscampp'  # or other method

    # binary paths (from config file)
    # these are just the default ones. User can modularize the method
    # for each step and add the corresponding executable paths
    # in main.config (by default at ~/.tipp3/main.config)
    pplacer_path = None
    bscampp_path = None
    blastn_path = None
    witch_path = None
    tippjsonmerger_path = None

    # reference package dir path
    refpkg_version = 'markers-v4'


    # miscellaneous
    alignment_only = False
    keeptemp = False
    bypass_setup = True
    
    ########## configs specific for download_refpkg ###########
    decompress = False

# check for valid configurations and set them
def set_valid_configuration(name, conf):
    if not isinstance(conf, Namespace):
        _LOG.warning('Looking for Namespace object from \'{}\' but find {}'.format(
            name, type(conf)))
        return

    # backbone alignment settings
    if name == 'basic':
        for k in conf.__dict__.keys():
            attr = getattr(conf, k)
            if not attr:
                continue

            if k == 'alignment_method':
                assert str(attr).lower() in ['witch', 'blast', 'hmm'], \
                    'Alignment method {} not implemented'.format(attr)
            elif k == 'placement_method':
                pass
                #assert int(attr).lower() in ['pplacer', 'bscampp'], \
                #    'Placement method {} not implemented'.format(attr)
            #elif k == 'path':
            #    assert os.path.exists(os.path.realpath(str(attr))), \
            #        '{} does not exist'.format(os.path.realpath(str(attr)))
            setattr(Configs, k, attr)
    elif name in valid_config_sections:
        setattr(Configs, name, conf)
    else:
        # not reading any invalid sections
        pass

# valid attribute check
def valid_attribute(k, v):
    assert isinstance(k, str)
    if isinstance(v, staticmethod):
        return False
    if not k.startswith('_'):
        return True
    return False

# print a list of all configurations
def getConfigs(arguments=None):
    msg = '\n********** Configurations **********\n' + \
            '\thome.path: {}\n'.format(homepath) + \
            '\tmain.config: {}\n'.format(main_config_path) + \
            '\targuments: {}\n\n'.format(arguments)
    for k, v in Configs.__dict__.items():
        if valid_attribute(k, v):
            msg += '\tConfigs.{}: {}\n'.format(k, v)
    print(msg, flush=True)

'''
Read in from config file if it exists. Any cmd-line provided configs will
override the config file.

Original functionality comes from SEPP -> sepp/config.py
'''
def _read_config_file(filename, cparser, opts,
        child_process=False, expand=None):
    #if not child_process:
    #    _LOG.info('Reading config from {}'.format(filename))
    config_defaults = []

    with open(filename, 'r') as cfile:
        cparser.read_file(cfile)
        if cparser.has_section('commandline'):
            for k, v in cparser.items('commandline'):
                config_defaults.append('--{}'.format(k))
                config_defaults.append(v)

        for section in cparser.sections():
            if section == 'commandline':
                continue
            if getattr(opts, section, None):
                section_name_space = getattr(opts, section)
            else:
                section_name_space = Namespace()
            for k, v in cparser.items(section):
                if expand and k == 'path':
                    v = os.path.join(expand, v)
                section_name_space.__setattr__(k, v)
            opts.__setattr__(section, section_name_space)
    return config_defaults

'''
Validate that Configs is set up correctly, mainly checking if the binary paths
exist
'''
def validateConfigs():
    b_valid = True
    ret = []
    files_to_check = {'PPLACER': Configs.pplacer_path,
            'BSCAMPP': Configs.bscampp_path,
            'BLASTN': Configs.blastn_path,
            'WITCH': Configs.witch_path}
    for method, path in files_to_check.items():
        if (not path) or (path == '') or (not os.path.exists(path)):
            b_valid = False
            ret.append((method, path))
    return b_valid, ret

'''
Build configurations
'''
def buildConfigs(parser, cmdline_args, child_process=False, rerun=False):
    # config parser, which first reads in main.config and later overrides
    # with user.config (if specified)
    cparser = configparser.ConfigParser()
    cparser.optionxform = str

    # load cmdline args first and identify the output directory
    # (to quickly create the outdir and log file
    args = parser.parse_args(cmdline_args)
    Configs.outdir = os.path.realpath(args.outdir)
    #if not os.path.exists(Configs.outdir):
    #    os.makedirs(Configs.outdir)

    # load default_args from main.config
    default_args = Namespace()
    cmdline_default = _read_config_file(main_config_path,
            cparser, default_args, child_process=child_process)
    
    # load cmdline args first, then search for user.config if specified
    args = parser.parse_args(cmdline_args)
    cmdline_user = []
    if 'config_file' in args.__dict__ and args.config_file != None:
        # override default_args
        Configs.config_file = args.config_file
        cmdline_user = _read_config_file(Configs.config_file,
                cparser, default_args, child_process=child_process)

    # finally, re-parse cmdline args in the order:
    #   [cmdline_default, cmd_user, cmdline_args] 
    args = parser.parse_args(cmdline_default + cmdline_user + cmdline_args,
            namespace=default_args)

    # store sub-command given
    Configs.command = args.command

    ######### subcommand: abundance ##########
    if Configs.command == 'abundance':
        # Must have
        Configs.query_path = os.path.realpath(args.query_path)
        if args.refpkg_path:
            Configs.refpkg_path = os.path.realpath(args.refpkg_path)

        Configs.alignment_only = args.alignment_only
        Configs.keeptemp = args.keeptemp

        # set up preset mode, TIPP3 or TIPP3-fast
        Configs.mode = args.mode
        if args.mode == 'tipp3-fast':
            Configs.alignment_method = 'blast'
            Configs.placement_method = 'bscampp'
        elif args.mode == 'tipp3':
            Configs.alignment_method = 'witch'
            Configs.placement_method = 'pplacer-taxtastic'

        # alignment_method and placement_method, and refpkg version
        if args.alignment_method:
            Configs.alignment_method = args.alignment_method
        if args.placement_method:
            Configs.placement_method = args.placement_method
        Configs.refpkg_version = args.refpkg_version

        if args.num_cpus > 0:
            Configs.num_cpus = min(os.cpu_count(), args.num_cpus)
        else:
            Configs.num_cpus = os.cpu_count()

        # verbose level
        verbose = os.getenv('TIPP_LOGGING_LEVEL', 'info').upper()
        if verbose in logging_levels:
            Configs.verbose = verbose
    ############ subcommand: download_refpkg #############
    elif Configs.command == 'download_refpkg':
        Configs.outdir = args.outdir
        Configs.decompress = args.decompress
    else:
        raise NotImplementedError
    
    #if args.max_concurrent_jobs:
    #    Configs.max_concurrent_jobs = args.max_concurrent_jobs
    #else:
    #    Configs.max_concurrent_jobs = min(50, 10 * Configs.num_cpus)

    # add any additional arguments to Configs
    for k in args.__dict__.keys():
        if k not in Configs.__dict__:
            k_attr = getattr(args, k)

            # check whether the configuration is valid
            set_valid_configuration(k, k_attr)

    # try once for validating Configs being set up correctly
    # if not, try once for re-initializing the main config file
    b_valid, invalid_paths = validateConfigs()
    if not b_valid and not rerun:
        _LOG.warning('Some software required by TIPP3 do not have valid binaries!')

        # trying once for regenerating
        _LOG.warning('Re-initializing the config file for once...')
        init_config_file(homepath, rerun=True)
        buildConfigs(parser, cmdline_args, rerun=True)
    elif not b_valid and rerun:
        # if re-initializing does not fix the problem
        if not b_valid:
            errmsg = 'Failed to find valid binaries for some software:\n'
            for item in invalid_paths:
                errmsg += f'\t{item[0]}: {item[1]}\n'
                _LOG.error(errmsg)
                raise ValueError(errmsg)

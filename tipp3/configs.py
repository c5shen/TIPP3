import os, time
try:
    import configparser
except ImportError:
    from ConfigParser import configparser
from argparse import ArgumentParser, Namespace
from platform import platform
from tipp3.init_configs import init_config_file
from tipp3 import get_logger

_LOG = get_logger(__name__)

# detect home.path location or creating one if it is missing
homepath = os.path.dirname(__file__) + '/home.path'
_root_dir, main_config_path = init_config_file(homepath)

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

    # basic input items
    query_path = None
    refpkg_path = None      # e.g., xxx/yyy/tipp3-refpkg
    outdir = None
    config_file = None     # Added @ 7.25.2024 
    keeptemp = False

    # choices of parameters
    alignment_method = 'witch'  # or blast
    placement_method = 'pplacer-taxtastic'  # or other method

    # binary paths (from config file)
    pplacer_path = None
    bscampp_path = None
    blastn_path = None
    witch_path = None

    # reference package dir path
    refpkg_version = 'markers-v4'

    # Multiprocessing settings
    num_cpus = -1
    max_concurrent_jobs = None

    # logging output path
    log_path = None


# check for valid configurations and set them
def set_valid_configuration(name, conf):
    assert isinstance(conf, Namespace), \
            'Looking for Namespace object but find {}'.format(type(conf))

    # backbone alignment settings
    if name == 'Basic':
        for k in conf.__dict__.keys():
            attr = getattr(conf, k)
            if not attr:
                continue

            if k == 'alignment_method':
                assert str(attr).lower() in ['witch', 'blast', 'hmm'], \
                    'Alignment method {} not implemented'.format(attr)
            elif k == 'placement_method':
                assert int(attr).lower() in ['pplacer-taxtastic', 'bscampp'], \
                    'Placement method {} not implemented'.format(attr)
            #elif k == 'path':
            #    assert os.path.exists(os.path.realpath(str(attr))), \
            #        '{} does not exist'.format(os.path.realpath(str(attr)))
            setattr(Configs, k, attr)
    elif name == 'WITCH':
        setattr(Configs, name, conf)
    elif name == 'Refpkg':
        setattr(Configs, name, conf)
    else:
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
def getConfigs():
    print('\n********* Configurations **********')
    print('\thome.path: {}'.format(homepath))
    print('\tmain.config: {}\n'.format(main_config_path))
    for k, v in Configs.__dict__.items():
        if valid_attribute(k, v):
            print('\tConfigs.{}: {}'.format(k, v))

'''
Read in from config file if it exists. Any cmd-line provided configs will
override the config file.

Original functionality comes from SEPP -> sepp/config.py
'''
def _read_config_file(filename, cparser, opts, expand=None):
    Configs.debug('Reading config from {}'.format(filename))
    config_defaults = []
    #cparser = configparser.ConfigParser()
    #cparser.optionxform = str
    cparser.read_file(filename)

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
Build configurations
'''
def buildConfigs(parser, cmdline_args):
    # config parser, which first reads in main.config and later overrides
    # with user.config (if specified)
    cparser = configparser.ConfigParser()
    cparser.optionxform = str

    # load default_args from main.config
    default_args = Namespace()
    cmdline_default = []
    with open(main_config_path, 'r') as cfile:
        cmdline_default = _read_config_file(cfile, cparser, default_args)
    
    # load cmdline args first, then search for user.config if specified
    args = parser.parse_args(cmdline_args)
    cmdline_user = []
    if args.config_file != None:
        # override default_args
        Configs.config_file = args.config_file
        with open(Configs.config_file, 'r') as cfile:
            cmdline_user = _read_config_file(cfile, cparser, default_args)

    # finally, re-parse cmdline args in the order:
    #   [cmdline_default, cmd_user, cmdline_args] 
    args = parser.parse_args(cmdline_default + cmdline_user + cmdline_args,
            namespace=default_args)

    # Must have
    Configs.query_path = os.path.realpath(args.query_path)
    Configs.refpkg_path = os.path.realpath(args.refpkg_path)
    Configs.outdir = os.path.realpath(args.outdir)

    if not os.path.exists(Configs.outdir):
        os.makedirs(Configs.outdir)
    Configs.log_path = os.path.join(Configs.outdir, 'tipp3.log')

    Configs.keeptemp = args.keeptemp

    # alignment_method and placement_method, and refpkg version
    Configs.alignment_method = args.alignment_method
    Configs.placement_method = args.placement_method
    Configs.refpkg_version = args.refpkg_version

    if args.num_cpus > 0:
        Configs.num_cpus = min(os.cpu_count(), args.num_cpus)
    else:
        Configs.num_cpus = os.cpu_count()

    if args.max_concurrent_jobs:
        Configs.max_concurrent_jobs = args.max_concurrent_jobs
    else:
        Configs.max_concurrent_jobs = min(50, 10 * Configs.num_cpus)

    # add any additional arguments to Configs
    for k in args.__dict__.keys():
        if k not in Configs.__dict__:
            k_attr = getattr(args, k)

            # check whether the configuration is valid
            set_valid_configuration(k, k_attr)

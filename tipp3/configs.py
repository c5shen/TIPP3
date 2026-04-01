import os
import configparser
from argparse import Namespace
from tipp3.init_configs import init_config_file
from tipp3 import get_logger

homepath = os.path.join(os.path.dirname(__file__), 'home.path')
_root_dir, main_config_path = init_config_file(homepath)

VALID_CONFIG_SECTIONS = ['witch', 'bscampp', 'pplacer-taxtastic',
                         'blast', 'refpkg']
LOGGING_LEVELS = {'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'}

_LOG = get_logger(__name__)


class Configs:
    """Central configuration store for TIPP3 pipeline settings."""

    command = None
    verbose = 'INFO'

    num_cpus = -1
    max_concurrent_jobs = None

    query_path = None
    refpkg_path = None
    outdir = None
    config_file = None

    mode = 'tipp3-fast'
    alignment_method = 'blast'
    placement_method = 'bscampp'

    pplacer_path = None
    bscampp_path = None
    blastn_path = None
    witch_path = None

    refpkg_version = 'markers-v4'

    bscampp_mode = None
    alignment_only = False
    keeptemp = False
    bypass_setup = True

    ########## configs specific for detection ##################
    detection_threshold = None
    
    ########## configs specific for download_refpkg ###########
    decompress = False

def set_valid_configuration(name, conf):
    """Validate and apply a named configuration section."""
    if not isinstance(conf, Namespace):
        _LOG.warning(f"Expected Namespace for '{name}', got {type(conf)}")
        return

    # backbone alignment settings
    if name == 'basic':
        for k in conf.__dict__:
            attr = getattr(conf, k)
            if not attr:
                continue
            if k == 'alignment_method':
                valid_methods = ('witch', 'blast', 'hmm')
                if str(attr).lower() not in valid_methods:
                    raise ValueError(
                        f"Alignment method '{attr}' not supported. "
                        f"Choose from: {valid_methods}")
            if attr is not None:
                setattr(Configs, k, attr)
    elif name in VALID_CONFIG_SECTIONS:
        setattr(Configs, name, conf)


def _is_public_attribute(k, v):
    """Check if an attribute should be displayed in config summary."""
    if isinstance(v, (staticmethod, classmethod)):
        return False
    return not k.startswith('_')


def getConfigs(arguments=None):
    """Log all current configuration values."""
    lines = [
        '\n********** Configurations **********',
        f'\thome.path: {homepath}',
        f'\tmain.config: {main_config_path}',
        f'\targuments: {arguments}\n',
    ]
    for k, v in Configs.__dict__.items():
        if _is_public_attribute(k, v):
            lines.append(f'\tConfigs.{k}: {v}')
    _LOG.info('\n'.join(lines))

def _read_config_file(filename, cparser, opts,
                      child_process=False, expand=None):
    """Read a config file and return command-line style defaults."""
    config_defaults = []

    with open(filename, 'r') as cfile:
        cparser.read_file(cfile)
        if cparser.has_section('commandline'):
            for k, v in cparser.items('commandline'):
                config_defaults.append(f'--{k}')
                config_defaults.append(v)

        for section in cparser.sections():
            if section == 'commandline':
                continue
            section_ns = getattr(opts, section, None) or Namespace()
            for k, v in cparser.items(section):
                if expand and k == 'path':
                    v = os.path.join(expand, v)
                setattr(section_ns, k, v)
            setattr(opts, section, section_ns)
    return config_defaults

def validateConfigs():
    """Validate that required binary paths exist for the chosen methods."""
    invalid = []

    # BLASTN is always required (used in binning step)
    required = {'BLASTN': Configs.blastn_path}

    # Only check tools required by the selected methods
    if Configs.alignment_method == 'witch':
        required['WITCH'] = Configs.witch_path
    if Configs.placement_method == 'bscampp':
        required['BSCAMPP'] = Configs.bscampp_path
    if Configs.placement_method in ('pplacer-taxtastic', 'bscampp'):
        required['PPLACER'] = Configs.pplacer_path

    for method, path in required.items():
        if not path or not os.path.exists(path):
            invalid.append((method, path))
    return len(invalid) == 0, invalid

def buildConfigs(parser, cmdline_args, child_process=False, rerun=False):
    """Build the full configuration from config files and command-line args."""
    cparser = configparser.ConfigParser()
    cparser.optionxform = str

    # Load defaults from main.config
    default_args = Namespace()
    cmdline_default = _read_config_file(
        main_config_path, cparser, default_args, child_process=child_process)

    # Check for user config file
    args = parser.parse_args(cmdline_args)
    cmdline_user = []
    if getattr(args, 'config_file', None) is not None:
        Configs.config_file = args.config_file
        cmdline_user = _read_config_file(
            Configs.config_file, cparser, default_args,
            child_process=child_process)

    # Re-parse in priority order: [defaults, user config, command-line]
    args = parser.parse_args(
        cmdline_default + cmdline_user + cmdline_args,
        namespace=default_args)

    Configs.command = args.command
    for k in args.__dict__:
        k_attr = getattr(args, k)
        if k in Configs.__dict__:
            if k_attr is not None:
                setattr(Configs, k, k_attr)
        else:
            set_valid_configuration(k, k_attr)

    verbose = os.getenv('TIPP_LOGGING_LEVEL', 'info').upper()
    if verbose in LOGGING_LEVELS:
        Configs.verbose = verbose

    if Configs.command in ('abundance', 'detection'):
        if args.num_cpus > 0:
            Configs.num_cpus = min(os.cpu_count(), args.num_cpus)
        else:
            Configs.num_cpus = os.cpu_count()

        mode_presets = {
            'tipp3-fast': ('blast', 'bscampp'),
            'tipp3': ('witch', 'pplacer-taxtastic'),
        }
        amethod, pmethod = mode_presets.get(Configs.mode, ('blast', 'bscampp'))

        if Configs.alignment_method is None:
            Configs.alignment_method = amethod
        if Configs.placement_method is None:
            Configs.placement_method = pmethod

    elif Configs.command == 'download_refpkg':
        pass
    else:
        raise ValueError(
            f"Unknown subcommand: '{Configs.command}'. "
            "Valid subcommands: abundance, detection, download_refpkg")

    # Validate binary paths; retry once by regenerating config
    b_valid, invalid_paths = validateConfigs()
    if not b_valid and not rerun:
        _LOG.warning('Some required binaries are missing, '
                     're-initializing config...')
        init_config_file(homepath, rerun=True)
        buildConfigs(parser, cmdline_args, rerun=True)
    elif not b_valid and rerun:
        errmsg = 'Failed to find valid binaries for the following:\n'
        for method, path in invalid_paths:
            errmsg += f'\t{method}: {path}\n'
        _LOG.error(errmsg)
        raise ValueError(errmsg)

"""
Initialize the TIPP3 configuration file (~/.tipp3/main.config).
"""

import os, sys, shutil, logging
import configparser
from platform import platform

_log = logging.getLogger(__name__)


def find_main_config(homepath):
    """Locate existing main.config from the home.path pointer file."""
    with open(homepath, 'r') as f:
        _root_dir = f.read().strip()
    main_config_path = os.path.join(_root_dir, 'main.config')
    if os.path.exists(main_config_path):
        return _root_dir, main_config_path
    return None, None


def init_config_file(homepath, rerun=False, prioritize_user_software=True):
    """Create or regenerate ~/.tipp3/main.config from default.config.

    Detects installed binaries (blastn, pplacer, bscampp, witch.py) from
    the user's PATH and records their locations.
    """
    bypass_setup = True

    if not rerun:
        if os.path.exists(homepath):
            if os.stat(homepath).st_mtime >= os.stat(__file__).st_mtime:
                _root_dir, main_config_path = find_main_config(homepath)
                if _root_dir is None:
                    _log.info('home.path exists but main.config missing, '
                              'regenerating...')
                else:
                    return _root_dir, main_config_path
            else:
                _log.info('Found outdated home.path, regenerating...')
                os.remove(homepath)
        else:
            _log.info(f'home.path not found at {homepath}, creating...')
    else:
        _log.info('Re-initializing the config file...')

    _root_dir = ''
    if not bypass_setup:
        _root_dir = input('Create main.config file at [default: ~/.tipp3/]: ')

    if _root_dir == '':
        _root_dir = os.path.expanduser('~/.tipp3')
    else:
        _root_dir = os.path.abspath(_root_dir)
    main_config_path = os.path.join(_root_dir, 'main.config')
    _log.info(f'Initializing main configuration: {main_config_path}')

    os.makedirs(_root_dir, exist_ok=True)
    with open(homepath, 'w') as f:
        f.write(_root_dir)

    _config_path = os.path.join(os.path.dirname(__file__), 'default.config')
    if not os.path.exists(_config_path):
        raise FileNotFoundError(
            f"Default config file missing: {_config_path}. "
            "Please reinstall TIPP3 or redownload from GitHub.")

    if os.path.exists(main_config_path):
        _log.info(f'Overwriting existing config: {main_config_path}')

    cparser = configparser.ConfigParser()
    cparser.optionxform = str

    default_config = configparser.ConfigParser()
    with open(_config_path, 'r') as f:
        default_config.read_file(f)
    for section in default_config.sections():
        cparser.add_section(section)
        for k, v in default_config[section].items():
            cparser.set(section, k, v)

    platform_name = platform()
    tools_dir = os.path.join(os.path.dirname(__file__), 'tools')

    cparser.set('basic', 'pplacer_path',
                os.path.join(tools_dir, 'pplacer', 'pplacer'))

    if 'macos' in platform_name.lower():
        _log.warning("macOS detected. Some bundled binaries may not work. "
                     "Please install pplacer, blastn, etc. manually.")
    else:
        _log.info(f'Platform: {platform_name}')

    if prioritize_user_software:
        _log.info('Detecting installed software from PATH...')
        software_map = {
            'bscampp': 'bscampp_path',
            'witch.py': 'witch_path',
            'pplacer': 'pplacer_path',
            'blastn': 'blastn_path',
        }
        for binary, config_key in software_map.items():
            found = shutil.which(binary)
            _log.info(f'  {binary}: {found or "not found"}')
            if found:
                cparser.set('basic', config_key, found)

    with open(main_config_path, 'w') as f:
        cparser.write(f)
    _log.info(f'main.config written to {main_config_path}')
    _log.info(f'To customize, edit {main_config_path} directly.')
    return _root_dir, main_config_path

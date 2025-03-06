###########################################################################
#    Copyright 2012 Siavash Mirarab, Nam Nguyen, and Tandy Warnow.
#    This file is part of SEPP.
#
#    SEPP is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SEPP is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with SEPP.  If not, see <http://www.gnu.org/licenses/>.
###########################################################################

from operator import itemgetter
import logging
import os

'''
Updated @ 3.6.2025 by Chengze Shen

Major changes to suit TIPP3 pipeline.
'''
__version__ = "0.1"
_INSTALL_PATH = __path__[0]

__all__ = ['read_binning', 'read_alignment', 'read_placement',
        'tipp3_pipeline', 'jobs', 'refpkg_downloader']

def get_setup_path():
    return _INSTALL_PATH


def get_logging_level(logging_level='info'):
    #### logging level map ####
    logging_level_map = {
            'DEBUG': logging.DEBUG, 'INFO': logging.INFO,
            'WARNING': logging.WARNING, 'ERROR': logging.ERROR,
            'CRITICAL': logging.CRITICAL,
            }
    # obtain environment variable to determine the logging level,
    # if TIPP3_LOGGING_LEVEL is not empty
    env_level = os.getenv('TIPP_LOGGING_LEVEL')
    if env_level is not None:
        ll = env_level.upper()
    else:
        ll = logging_level.upper()
    return logging_level_map.get(ll, logging.INFO) 
    #return logging.DEBUG if _DEBUG else logging.INFO

__set_loggers = set()

def get_logger(name="tipp3", log_path=None, logging_level='info'):
    logger = logging.getLogger(name)
    if name not in __set_loggers:
        level = get_logging_level(logging_level)
        logging_formatter = logging.Formatter(
            ("[%(asctime)s] %(filename)s (line %(lineno)d):"
             " %(levelname) 8s: %(message)s"))
        logging_formatter.datefmt = '%H:%M:%S'
        logger.setLevel(level)
        if log_path == None:
            ch = logging.StreamHandler()
        else:
            # use FileHandler instead
            ch = logging.FileHandler(log_path, mode='a')

        ch.setLevel(level)
        ch.setFormatter(logging_formatter)
        logger.addHandler(ch)
        __set_loggers.add(name)
    return logger


#### unused for TIPP3 for now
#def reset_loggers():
#    global __set_loggers
#    __set_loggers = set()
#    import pkgutil
#    import sepp
#    package = sepp
#    for modl, name, _ in pkgutil.iter_modules(package.__path__):
#        logger = (getattr(getattr(sepp, name, None), "_LOG", None))
#        print("--- *", name, logger)
#        if logger:
#            setattr(getattr(sepp, name, None), "_LOG", get_logger(
#                "tipp3.%s" % name))


def log_exception(logger):
    """Logs the exception trace to the logObj as an error"""
    import traceback
    import io
    s = io.StringIO()
    traceback.print_exc(None, s)
    logger.error(s.getvalue())
    exit(1)

os.sys.setrecursionlimit(1000000)


def sort_by_value(d, reverse=False):
    return sorted(iter(d.items()), key=itemgetter(1), reverse=reverse)

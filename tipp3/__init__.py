"""
TIPP3: Taxonomic Identification and Phylogenetic Profiling version 3.

Originally adapted from SEPP (Mirarab, Nguyen, Warnow 2012).
Substantially rewritten for TIPP3 by Chengze Shen (2024-2025).
"""

import logging
import os
import sys

__version__ = "0.5"

__all__ = ['query_binning', 'query_alignment', 'query_placement',
           'tipp3_pipeline', 'jobs', 'refpkg_loader']


def get_logging_level(logging_level='info'):
    """Resolve the effective logging level from env var or argument."""
    level_map = {
        'DEBUG': logging.DEBUG, 'INFO': logging.INFO,
        'WARNING': logging.WARNING, 'ERROR': logging.ERROR,
        'CRITICAL': logging.CRITICAL,
    }
    env_level = os.getenv('TIPP_LOGGING_LEVEL')
    ll = env_level.upper() if env_level is not None else logging_level.upper()
    return level_map.get(ll, logging.INFO)


_configured_loggers = set()


def get_logger(name="tipp3", log_path=None, logging_level='info'):
    """Get or create a named logger with consistent formatting."""
    logger = logging.getLogger(name)
    if name not in _configured_loggers:
        level = get_logging_level(logging_level)
        formatter = logging.Formatter(
            "[%(asctime)s] %(filename)s (line %(lineno)d): "
            "%(levelname)8s: %(message)s",
            datefmt='%H:%M:%S')
        logger.setLevel(level)

        if log_path is None:
            handler = logging.StreamHandler()
        else:
            handler = logging.FileHandler(log_path, mode='a')

        handler.setLevel(level)
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        _configured_loggers.add(name)
    return logger


sys.setrecursionlimit(1000000)

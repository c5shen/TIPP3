import os
from tipp3.configs import Configs
from tipp3 import get_logger

_LOG = get_logger(__name__)

'''
Load TIPP3 reference package in
'''
def loadReferencePackage(refpkg_path, refpkg_version):
    refpkg = {}

    # sanity check for the existence of refpkg_path
    if not refpkg_path or not os.path.exists(refpkg_path):
        errmsg = 'Refpkg does not exist: {}'.format(refpkg_path)
        _LOG.error(errmsg)
        raise ValueError(errmsg)

    # refpkg dir path from commandline
    path = os.path.join(refpkg_path, refpkg_version)
    input = os.path.join(path, "file-map-for-tipp.txt")
    _LOG.info('Reading refpkg from {}'.format(path))

    # load exclusion list, if any
    exclusion = set() 
    try:
        raw = Configs.Refpkg.exclusion
        exclusion = set(raw.strip().split(','))
    except AttributeError:
        pass

    refpkg["genes"] = []
    with open(input) as f:
        for line in f.readlines():
            [key, val] = line.split('=')

            [key1, key2] = key.strip().split(':')
            val = os.path.join(path, val.strip())

            try:
                refpkg[key1][key2] = val
            except KeyError:
                refpkg[key1] = {}
                refpkg[key1][key2] = val

            if (key1 != "blast") and (key1 != "taxonomy"):
                refpkg["genes"].append(key1)

    # excluding marker genes if specified
    _LOG.info('Excluding markers (if exist): {}'.format(exclusion))
    refpkg["genes"] = set(refpkg["genes"]).difference(exclusion)
    refpkg["genes"] = list(refpkg["genes"])
    _LOG.info('Marker genes: {}'.format(refpkg["genes"]))
    _LOG.info('Number of marker genes: {}'.format(len(refpkg["genes"])))

    return refpkg

import os
from tipp3.configs import Configs
from tipp3 import get_logger

_LOG = get_logger(__name__)

'''
Load TIPP3 reference package in
'''
def load_reference_package():
    refpkg = {}

    # refpkg dir path from commandline
    path = os.path.join(Configs.refpkg_path, Configs.refpkg_version)
    input = os.path.join(path, "file-map-for-tipp.txt")
    _LOG.INFO('Reading refpkg from {}'.format(path))

    # load exclusion list, if any
    exclusion = set() 
    if getattr(Configs, 'Refpkg', None) != None:
        if 'exclusion' in Configs.Refpkg.__dict__:
            exclusion = set(Configs.Refpkg['exclusion'])

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
    refpkg["genes"] = set(refpkg["genes"]).difference(exclusion)
    refpkg["genes"] = list(refpkg["genes"])
    _LOG.INFO('Marker genes: {}'.format(refpkg["genes"]))
    _LOG.INFO('Number of marker genes: {}'.format(len(refpkg["genes"])))

    return refpkg

import os
import zipfile
from tipp3.configs import Configs
from tipp3 import get_logger, log_exception

_LOG = get_logger(__name__)

'''
Load TIPP3 reference package in
'''
def loadReferencePackage(refpkg_path, refpkg_version):
    refpkg = {}

    # sanity check for the existence of refpkg_path
    if not refpkg_path or not os.path.exists(refpkg_path):
        errmsg = ('Refpkg does not exist: {}'.format(refpkg_path) + 
            '\nPlease download reference package using subcommand ' +
            '\"download_refpkg\"')
        _LOG.error(errmsg)
        raise ValueError(errmsg)

    # refpkg dir path from commandline
    path = os.path.join(refpkg_path, refpkg_version)
    input = os.path.join(path, "file-map-for-tipp.txt")
    _LOG.info('Reading refpkg from {}'.format(path))

    # load exclusion list, if any
    exclusion = set() 
    try:
        raw = getattr(Configs, 'refpkg').exclusion
        exclusion = set(raw.strip().split(','))
    except AttributeError:
        pass

    refpkg["genes"] = []
    with open(input) as f:
        for line in f.readlines():
            [key, val] = line.split('=')

            [key1, key2] = key.strip().split(':')

            # hotfix before pushing a new version of TIPP3 refpkg
            # --> change all "taxonomy.table" to "all_taxon.taxonomy"
            if val == 'taxonomy.table':
                val = 'all_taxon.taxonomy'
            val = os.path.join(path, val.strip())

            try:
                refpkg[key1][key2] = val
            except KeyError:
                refpkg[key1] = {}
                refpkg[key1][key2] = val

            if (key1 != "blast") and (key1 != "taxonomy"):
                refpkg["genes"].append(key1)
    
    # add path variable to each marker gene refpkg
    # to use with pplacer-taxtastic
    for marker in refpkg["genes"]:
        marker_refpkg_path = os.path.join(path, f"{marker}.refpkg")
        refpkg[marker]['path'] = marker_refpkg_path

    # excluding marker genes if specified
    _LOG.info('Excluding markers (if exist): {}'.format(exclusion))
    refpkg["genes"] = set(refpkg["genes"]).difference(exclusion)
    refpkg["genes"] = list(refpkg["genes"])
    _LOG.info('Marker genes: {}'.format(refpkg["genes"]))
    _LOG.info('Number of marker genes: {}'.format(len(refpkg["genes"])))

    return refpkg

'''
Download the latest TIPP reference package
'''
def downloadReferencePackage(outdir, decompress=False):
    # hosted on Illinois Databank
    url = 'https://databank.illinois.edu/datafiles/yq7ry/download'
    _LOG.info(f"Downloading the latest TIPP reference package from {url}")

    # specify the name of the downloaded zipfile
    outpath = os.path.join(outdir, 'tipp3-refpkg.zip')
    if os.path.exists(outpath):
        _LOG.warning(
                f"{outpath} exists! Skipped download to avoid overwriting.")
    else:
        cmd = f"wget {url} -O {outpath}"
        os.system(cmd)

    # identify the extract directory name
    try:
        zf = zipfile.ZipFile(outpath)
    except FileNotFoundError as e:
        log_exception(_LOG)
    files = zf.namelist()
    dir = files[0].split('/')[0]

    # decompress if specified
    if decompress:
        _LOG.info(f"Decompressing {os.path.basename(outpath)} to {outdir}")
        if os.path.isdir(dir):
            _LOG.warning(f"\"{dir}/\" from the zipfile exists at {outdir}! "
                    "Skipped decompression to avoid overwriting.")
        else:
            cmd = f"unzip -d {outdir} {outpath}"
            os.system(cmd)
            _LOG.info(f"Decompressed refpkg directory name: {dir}")
    else:
        _LOG.info("Finished downloading. To use with TIPP3, please "
            f"decompress the downloaded file at {outpath}")
    return True

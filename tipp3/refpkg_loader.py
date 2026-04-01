"""
Load and download TIPP3 reference packages.
"""

import os, subprocess, zipfile
from tipp3.configs import Configs
from tipp3 import get_logger

_LOG = get_logger(__name__)


def loadReferencePackage(refpkg_path, refpkg_version):
    """Load a TIPP3 reference package from disk."""
    refpkg = {}

    if not refpkg_path or not os.path.exists(refpkg_path):
        raise FileNotFoundError(
            f"Reference package not found: {refpkg_path}\n"
            "Download it with: tipp3 download_refpkg -d <path> --decompress")

    path = os.path.join(refpkg_path, refpkg_version)
    filemap = os.path.join(path, "file-map-for-tipp.txt")
    if not os.path.exists(filemap):
        raise FileNotFoundError(
            f"file-map-for-tipp.txt not found in {path}. "
            f"Is '{refpkg_version}' the correct refpkg version?")
    _LOG.info(f'Loading refpkg from {path}')

    exclusion = set()
    try:
        raw = getattr(Configs, 'refpkg').exclusion
        exclusion = set(raw.strip().split(','))
    except AttributeError:
        pass

    refpkg["genes"] = []
    with open(filemap) as f:
        for line in f:
            line = line.strip()
            if not line or '=' not in line:
                continue
            key, val = line.split('=', 1)
            key1, key2 = key.strip().split(':', 1)

            if val.strip() == 'taxonomy.table':
                val = 'all_taxon.taxonomy'
            val = os.path.join(path, val.strip())

            if key1 not in refpkg:
                refpkg[key1] = {}
            refpkg[key1][key2] = val

            if key1 not in ("blast", "taxonomy"):
                refpkg["genes"].append(key1)

    for marker in refpkg["genes"]:
        refpkg[marker]['path'] = os.path.join(path, f"{marker}.refpkg")

    if exclusion:
        _LOG.info(f'Excluding markers: {exclusion}')
    refpkg["genes"] = list(set(refpkg["genes"]) - exclusion)
    _LOG.info(f'Loaded {len(refpkg["genes"])} marker genes: {refpkg["genes"]}')

    return refpkg


def downloadReferencePackage(outdir, decompress=False):
    """Download the latest TIPP3 reference package from Illinois Databank."""
    latest_version = 'tipp3-refpkg-1-2.zip'
    url = 'https://databank.illinois.edu/datafiles/sarfb/download'

    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, latest_version)

    if os.path.exists(outpath):
        _LOG.info(f"{outpath} already exists, skipping download.")
    else:
        _LOG.info(f"Downloading TIPP3 reference package from {url}")
        result = subprocess.run(
            ['wget', url, '-O', outpath],
            capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(
                f"Download failed (exit code {result.returncode}):\n"
                f"{result.stderr}")

    try:
        zf = zipfile.ZipFile(outpath)
    except (FileNotFoundError, zipfile.BadZipFile) as e:
        raise RuntimeError(
            f"Cannot open downloaded file {outpath}: {e}") from e

    files = zf.namelist()
    refpkg_dir = files[0].split('/')[0]

    if decompress:
        target = os.path.join(outdir, refpkg_dir)
        if os.path.isdir(target):
            _LOG.info(f"'{refpkg_dir}/' already exists, skipping extraction.")
        else:
            _LOG.info(f"Extracting {os.path.basename(outpath)} to {outdir}")
            zf.extractall(outdir)
            _LOG.info(f"Extracted refpkg: {refpkg_dir}")
    else:
        _LOG.info(f"Downloaded to {outpath}. "
                  f"To use, decompress with: unzip -d {outdir} {outpath}")

    zf.close()
    return True

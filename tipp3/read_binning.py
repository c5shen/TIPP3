'''
Read binning using BLASTN
'''

import os, time
from tipp3 import get_logger
from tipp3.configs import Configs
from tipp3.jobs import BlastnJob

_LOG = get_logger(__name__)

def read_binning(refpkg):
    database_path = refpkg['blast']['database']
    _LOG.INFO(f"Initializing BlastnJob with BLASTN database: {database_path}")

    blastn_outdir = os.path.join(Configs.outdir, 'blast_output')
    if not os.path.isdir(blastn_outdir):
        os.makedirs(blastn_outdir)

    # run BLASTN
    job = BlastnJob(path=Configs.blastn_path,
            query_path=Configs.query_path,
            database_path=database_path,
            outdir=blastn_outdir,
            num_threads=Configs.num_cpus)
    blastn_outpath = job.run()

    # process BLASTN output
    #TODO

'''
A collection of jobs for TIPP3.

Jobs are designed to be run standalone, as long as all parameters
are provided correctly.
'''

import os, shutil, subprocess, stat, re, traceback
from subprocess import Popen
from tipp3 import get_logger

_LOG = get_logger(__name__)

'''
A BLASTN job that will run BLASTN to bin reads against the reference package
marker genes
'''
class BlastnJob:
    def __init__(self, **kwargs):
        self.job_type = 'blastn'
        self.outfmt = 0

        # arguments for running BLASTN
        self.path = None 
        self.query_path = None
        self.database_path = None
        self.outdir = None
        self.num_threads = None 

        for k, v in kwargs.items():
            setattr(self, k, v)

    def run(self):
        #blast:database = blast/alignment.fasta.db
        #blast:seq-to-marker-map = blast/seq2marker.tab
        outpath = os.path.join(self.outdir, 'blast.alignment.out') 
        cmd = [self.path, '-db', self.database_path,
                '-outfmt', str(self.outfmt),
                '-query', self.query_path,
                '-out', outpath,
                '-num_threads', str(self.num_threads)]
        p = Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        retcode = p.returncode
        stdout, stderr = p.communicate()

        if retcode == 0:
            _LOG.INFO("BLASTN completed.")
            return outpath
        else:
            error_msg = ' '.join(["Error occurred when running BLASTN.",
                f"Return code: {retcode}"])
            print(error_msg)
            _LOG.ERROR(error_msg + '\n' + stderr.decode('utf-8'))
            exit(1)

'''
A BSCAMPP job that will replace the original PplacerJob
'''
class BscamppJob(ExternalSeppJob):
    def __init__(self, **kwargs):
        self.job_type = 'bscampp'
        ExternalSeppJob.__init__(self, self.job_type, **kwargs)
        
        # initialize parameters
        self.backbone_alignment_file = None
        self.tree_file = None
        self.model_file = None
        self.extended_alignment_file = None
        self.full_extended_alignment_file = None
        self.out_file = None

    def setup(self, backbone_alignment_file, tree_file, model_file,
            extended_alignment_file, full_extended_alignment_file,
            out_file, **kwargs):
        self.backbone_alignment_file = backbone_alignment_file
        self.tree_file = tree_file
        self.model_file = model_file
        self.extended_alignment_file = extended_alignment_file
        self.full_extended_alignment_file = full_extended_alignment_file
        self.out_file = out_file

    def partial_setup_for_subproblem(self, subproblem, model_file, i, **kwargs):
        assert isinstance(subproblem, sepp.problem.SeppProblem)
        self.backbone_alignment_file = sepp.filemgr.tempfile_for_subproblem(
                'bscampp.backbone.', subproblem, '.fasta')
        self.tree_file = sepp.filemgr.tempfile_for_subproblem(
                'bscampp.tree.', subproblem, '.tre')
        self.extended_alignment_file = \
                sepp.filemgr.tempfile_for_subproblem(
                        'bscampp.extended.{}.'.format(i),
                        subproblem, '.fasta')
        self.full_extended_alignment_file = \
                sepp.filemgr.tempfile_for_subproblem(
                        'bscampp.full.extended.{}.'.format(i),
                        subproblem, '.fasta')
        self.out_file = os.path.join(
                sepp.filemgr.tempdir_for_subproblem(subproblem),
                self.extended_alignment_file.replace('fasta', 'jplace'))
        assert isinstance(subproblem.subtree, PhylogeneticTree)
        subproblem.subtree.write_newick_to_path(self.tree_file,)

        self.model_file = model_file.name \
                if hasattr(model_file, 'name') else model_file
        self.tmpfilenbr = i
        self._kwargs = kwargs
    
    def get_invocation(self):
        invoc = ['python3', self.path, '-d', os.path.dirname(self.out_file)]
        if 'user_options' in self._kwargs:
            invoc.extend(self._kwargs['user_options'].split())

        invoc.extend(['-o', '.'.join(self.out_file.split('/')[-1].split('.')[:-1]),
                      '-t', self.tree_file,
                      '-i', self.model_file,
                      '-a', self.full_extended_alignment_file,
                      '--tmpfilenbr', str(self.tmpfilenbr),
                      '-b', '100', '--threads', '1'])
        return invoc
    
    def characterize_input(self):
        return ('backbone_alignment_file:%s, tree_file:%s, model_file:%s, '
                'full_extended_alignment_file:%s, output:%s') % (
                    self.backbone_alignment_file, self.tree_file,
                    self.model_file, self.full_extended_alignment_file,
                    self.out_file)

    def read_results(self):
        if self.fake_run:
            return None
        assert os.path.exists(self.out_file)
        assert os.stat(self.out_file)[stat.ST_SIZE] != 0
        return self.out_file

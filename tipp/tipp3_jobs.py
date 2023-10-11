from sepp.jobs import ExternalSeppJob
from sepp.scheduler import Job, JobError
from sepp import get_logger
from sepp.tree import PhylogeneticTree
from abc import abstractmethod
from subprocess import Popen
import sepp

import os, shutil, subprocess, stat, re traceback
import sepp.config

import io
try:
    filetypes = (io.IOBase, file)
except NameError:
    filetypes = io.IOBase

_LOG = get_logger(__name__)

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
        self.info_file = None
        self.extended_alignment_file = None
        self.full_extended_alignment_file = None
        self.out_file = None

    def setup(self, backbone_alignment_file, tree_file, info_file,
            extended_alignment_file, full_extended_alignment_file,
            out_file, **kwargs):
        self.backbone_alignment_file = backbone_alignment_file
        self.tree_file = tree_file
        self.info_file = info_file
        self.extended_alignment_file = extended_alignment_file
        self.full_extended_alignment_file = full_extended_alignment_file
        self.out_file = out_file

    def partial_setup_for_subproblem(self, subproblem, info_file, i, **kwargs):
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

        self.info_file = info_file.name \
                if hasattr(info_file, 'name') else info_file
        self._kwargs = kwargs
    
    def get_invocation(self):
        invoc = ['python3', self.path, '-d', os.path.dirname(self.out_file)]
        if 'user_options' in self._kwargs:
            invoc.extend(self._kwargs['user_options'].split())

        invoc.extend(['-o', '.'.join(self.out_file.split('/')[-1].split('.')[:-1]),
                      '-t', self.tree_file,
                      '-i', self.info_file,
                      '-a', self.full_extended_alignment_file,
                      '-b', '100', '--threads', '1'])
        return invoc
    
    def characterize_input(self):
        return ('backbone_alignment_file:%s, tree_file:%s, info_file:%s, '
                'full_extended_alignment_file:%s, output:%s') % (
                    self.backbone_alignment_file, self.tree_file,
                    self.info_file, self.full_extended_alignment_file,
                    self.out_file)

    def read_results(self):
        if self.fake_run:
            return None
        assert os.path.exists(self.out_file)
        assert os.stat(self.out_file)[stat.ST_SIZE] != 0
        return self.out_file

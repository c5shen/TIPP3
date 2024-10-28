'''
A collection of jobs for TIPP3.

Jobs are designed to be run standalone, as long as all parameters
are provided correctly.
'''

import os, shutil, subprocess, stat, re, traceback
import threading
from subprocess import Popen
from abc import abstractmethod

from tipp3 import get_logger
from tipp3.configs import Configs

_LOG = get_logger(__name__)

'''
function to streamline the logging of stdout and stderr output from a job run
to a target logging file
'''
def stream_to_file(stream, fptr, logging):
    if logging:
        for line in iter(stream.readline, ''):
            fptr.write(line)
            fptr.flush()
    stream.close()

'''
Template Class Job for running external jobs
'''
class Job(object):
    def __init__(self):
        self.job_type = ""
        self.errors = []
        self.b_ignore_error = False
        self.pid = -1
        self.returncode = 0

    def __call__(self):
        return self.run()

    def get_pid(self):
        return self.pid

    # run the job with given invocation defined in a child class
    # raise errors when encountered
    def run(self, stdin="", lock=None, logging=False):
        try:
            cmd, outpath = self.get_invocation()
            _LOG.info(f"Running job_type: {self.job_type}, output: {outpath}")

            binpath = cmd[0]
            # special case for "java -jar ..."
            if binpath == 'java':
                binpath = cmd[2]

            assert os.path.exists(binpath), \
                    ("executable for %s does not exist: %s" % 
                     (self.job_type, binpath))
            assert \
                (binpath.count("/") == 0 or os.path.exists(binpath)), \
                ("path for %s does not exist (%s)" %
                 (self.job_type, binpath))

            _LOG.debug("Arguments: %s", " ".join(
                (str(x) if x is not None else "?NoneType?"
                 for x in cmd)))
        
            # run the job using Popen with given command
            # read in stdout and stderr
            p = Popen(cmd, text=True, stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.pid = p.pid
            stdout, stderr = p.communicate(input=stdin)
            self.returncode = p.returncode

            #p = Popen(cmd, bufsize=1, text=True,
            #        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            #_outdir = os.path.dirname(os.path.realpath(outpath))
            #if logging:
            #    logfile = os.path.join(_outdir, 'runtime.txt')
            #    fptr = open(logfile, 'w', buffering=1)

            #    # initialize writing to logging file (if logging)
            #    stdout_thread = threading.Thread(target=stream_to_file,
            #            args=(p.stdout, fptr, logging))
            #    stderr_thread = threading.Thread(target=stream_to_file,
            #            args=(p.stderr, fptr, logging))
            #    stdout_thread.start()
            #    stderr_thread.start()

            #    # join threads
            #    stdout_thread.join()
            #    stderr_thread.join()

            # finish up the process run
            #stdout, stderr = p.communicate()
            #self.returncode = p.returncode

            if self.returncode == 0:
                if lock:
                    try:
                        lock.acquire()
                        _LOG.info(f"{self.job_type} completed, output: {outpath}")
                    finally:
                        lock.release()
                else:
                    _LOG.info(f"{self.job_type} completed, output: {outpath}")
                return outpath
            else:
                error_msg = ' '.join([f"Error occurred running {self.job_type}.",
                    f"Return code: {self.returncode}"])
                print(error_msg)
                if lock:
                    try:
                        lock.acquire()
                        _LOG.error(error_msg + '\nSTDOUT: ' + stdout.decode('utf-8') +
                                '\nSTDERR: ' + stderr.decode('utf-8'))
                    finally:
                        lock.release()
                else:
                    _LOG.error(error_msg + '\nSTDOUT: ' + stdout.decode('utf-8') +
                            '\nSTDERR: ' + stderr.decode('utf-8'))
                #_LOG.error(error_msg + '\n' + stdout)
                exit(1)
        except Exception:
            traceback.print_exc()
            raise

    # implement in child class
    # return: (cmd, outpath)
    @abstractmethod
    def get_invocation(self):
        raise NotImplementedError(
            "get_invocation() should be implemented by subclasses")

'''
A BLASTN job that will run BLASTN to bin reads against the reference package
marker genes
'''
class BlastnJob(Job):
    def __init__(self, **kwargs):
        Job.__init__(self)
        self.job_type = 'blastn'
        self.outfmt = 0

        # arguments for running BLASTN
        self.path = '' 
        self.query_path = ''
        self.database_path = ''
        self.outdir = ''
        self.num_threads = 1 

        for k, v in kwargs.items():
            setattr(self, k, v)

    def get_invocation(self):
        #blast:database = blast/alignment.fasta.db
        #blast:seq-to-marker-map = blast/seq2marker.tab
        self.outpath = os.path.join(self.outdir, 'blast.alignment.out') 
        cmd = [self.path, '-db', self.database_path,
                '-outfmt', str(self.outfmt),
                '-query', self.query_path,
                '-out', self.outpath,
                '-num_threads', str(self.num_threads)]
        return cmd, self.outpath

'''
A BSCAMPP job that will run BSCAMPP for placing aligned query reads 
'''
class BscamppJob(Job):
    def __init__(self, **kwargs):
        Job.__init__(self)
        self.job_type = 'bscampp'
        
        # initialize parameters
        self.path = ''
        self.query_alignment_path = ''
        self.backbone_tree_path = ''
        self.tree_model_path = ''
        self.outdir = ''
        self.num_cpus = 1
        self.subtreesize = 2000

        for k, v in kwargs.items():
            setattr(self, k, v)

    def get_invocation(self):
        self.outpath = os.path.join(self.outdir, 'placement.jplace')
        cmd = [self.path, '-b', str(self.subtreesize),
                '-i', self.tree_model_path,
                '-t', self.backbone_tree_path,
                '-d', self.outdir, '-o', 'placement',
                '-a', self.query_alignment_path,
                '--threads', str(self.num_cpus)]
        return cmd, self.outpath 

'''
A pplacer-taxtastic job that runs pplacer with the taxtastic refpkg
'''
class PplacerTaxtasticJob(Job):
    def __init__(self, **kwargs):
        Job.__init__(self)
        self.job_type = 'pplacer-taxtastic'

        self.path = ''
        self.query_alignment_path = ''
        self.refpkg_path = ''
        self.outdir = ''
        self.num_cpus = 1
        self.model_type = 'GTR'

        for k, v in kwargs.items():
            setattr(self, k, v)

    def get_invocation(self):
        self.outpath = os.path.join(self.outdir, 'placement.jplace')
        cmd = [self.path, '-m', self.model_type,
                '-c', self.refpkg_path,
                '-o', self.outpath,
                '-j', str(self.num_cpus),
                self.query_alignment_path]
        return cmd, self.outpath

'''
A JsonMerger job for TIPP to read taxonomic information from a .jplace file
with the corresponding taxonomic tree
'''
class TIPPJsonMergerJob(Job):
    def __init__(self, **kwargs):
        Job.__init__(self)
        self.job_type = 'tippjsonmerger'

        self.path = ''
        self.taxonomy_path = ''
        self.mapping_path = ''
        self.classification_path = ''
        self.outdir = ''

        for k, v in kwargs.items():
            setattr(self, k, v)

    def get_invocation(self):
        self.outpath = os.path.join(self.outdir, 'tippjsonmerger.jplace') 
        cmd = ['java', '-jar', self.path,
                '-', '-', self.outpath,
                '-t', self.taxonomy_path, '-m', self.mapping_path,
                '-p', '0.0', '-C', '0.0', '-c', self.classification_path]
        return cmd, self.classification_path

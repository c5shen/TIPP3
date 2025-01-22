'''
A collection of jobs for TIPP3.

Jobs are designed to be run standalone, as long as all parameters
are provided correctly.
'''

import os, shutil, subprocess, stat, re, traceback, shlex
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
    def run(self, stdin="", lock=None, logging=False, shell=False):
        try:
            cmd, outpath = self.get_invocation()
            _LOG.debug(f"Running job_type: {self.job_type}, output: {outpath}")

            # failsafe for NotImplemented jobs
            if len(cmd) == 0:
                raise ValueError(
                        f"{self.job_type} does not have a valid run command. "
                        "It might be due to (invalid input types, ).")

            binpath = cmd[0]
            # special case for "java -jar ..."
            if binpath == 'java':
                binpath = cmd[2]

            assert os.path.exists(binpath) or binpath == "gzip", \
                    ("executable for %s does not exist: %s" % 
                     (self.job_type, binpath))
            assert \
                (binpath.count("/") == 0 or os.path.exists(binpath)), \
                ("path for %s does not exist (%s)" %
                 (self.job_type, binpath))

            _LOG.debug("Arguments: %s", " ".join(
                (str(x) if x is not None else "?NoneType?"
                 for x in cmd)))
        
            # logging to a local file at the target outdir
            # logname: <outdir>/{self.job_type}.txt
            stdout, stderr = '', ''
            scmd = ' '.join(cmd)
            if logging:
                logpath = os.path.join(
                        os.path.dirname(outpath), f'{self.job_type}.txt')
                outlogging = open(logpath, 'w', 1)

                # deal with piping between multiple commands (if any)
                if '|' in scmd:
                    _stdout = subprocess.PIPE
                    subcmds = [shlex.split(x) for x in scmd.split('|')]
                    prev_p = Popen(subcmds[0], text=True, bufsize=1,
                            stdout=_stdout)
                    for i in range(1, len(subcmds)):
                        if i == len(subcmds) - 1:
                            _stdout = outlogging

                        p = Popen(subcmds[i], text=True, bufsize=1,
                            stdin=prev_p.stdout,
                            stdout=_stdout, stderr=subprocess.PIPE)
                        self.pid = p.pid
                        prev_p = p
                    stdout, stderr = p.communicate()
                    stdout = ''
                else:
                    p = Popen(cmd, text=True, bufsize=1,
                            stdin=subprocess.PIPE,
                            stdout=outlogging, stderr=subprocess.PIPE)
                    self.pid = p.pid
                    stdout, stderr = p.communicate(input=stdin)
                    stdout = ''
                outlogging.close()
            else:
                if '|' in scmd:
                    subcmds = [shlex.split(x) for x in scmd.split('|')]
                    prev_p = Popen(subcmds[0], text=True, bufsize=1,
                            stdout=subprocess.PIPE)
                    for i in range(1, len(subcmds)):
                        p = Popen(subcmds[i], text=True, bufsize=1,
                            stdin=prev_p.stdout,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        self.pid = p.pid
                        prev_p = p
                    stdout, stderr = p.communicate()
                else:
                    p = Popen(cmd, text=True, bufsize=1,
                            stdin=subprocess.PIPE,
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
                        _LOG.debug(f"{self.job_type} completed, output: {outpath}")
                    finally:
                        lock.release()
                else:
                    _LOG.debug(f"{self.job_type} completed, output: {outpath}")
                return outpath
            else:
                error_msg = ' '.join([f"Error occurred running {self.job_type}.",
                    f"Return code: {self.returncode}"])
                print(error_msg)
                if lock:
                    try:
                        lock.acquire()
                        _LOG.error(error_msg + '\nSTDOUT: ' + stdout +
                                '\nSTDERR: ' + stderr)
                    finally:
                        lock.release()
                else:
                    _LOG.error(error_msg + '\nSTDOUT: ' + stdout +
                            '\nSTDERR: ' + stderr)
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

        # check input type: if .fasta as suffix then it is fine
        # if .fasta.gz suffix, then use the file as stdin 
        name_parts = self.query_path.split('.')
        suffix = name_parts[-1]

        cmd = []

        # normal input type with fasta/fa files
        if suffix in ['fa', 'fasta']:
            cmd = [self.path, '-db', self.database_path,
                    '-outfmt', str(self.outfmt),
                    '-query', self.query_path,
                    '-out', self.outpath,
                    '-num_threads', str(self.num_threads)]
        # awk to process fq and fastq files as input
        elif suffix in ['fq', 'fastq']:
            cmd.extend(['awk',
                '\'NR%4==1 {print \">\"substr($0, 2)} NR%4==2 {print $0}\'',
                '|'])
        # gzip to process gzip and gz files as input
        elif suffix in ['gz', 'gzip']: 
            cmd.extend(['gzip', '-dc', self.query_path, '|'])
            # check if this is fasta/fa gzip or fastq/fq gzip
            if len(name_parts) > 2:
                suffix2 = name_parts[-2]

            # if we have another layer of fastq/fq to deal with, extend
            # the cmd to deal with that
            if suffix2 in ['fastq', 'fq']:
                cmd.extend(['awk',
                    '\'NR%4==1 {print \">\"substr($0, 2)} NR%4==2 {print $0}\'',
                    '|'])

            # piped query will go in to BLAST as stdin
            cmd.extend([self.path, '-db', self.database_path,
                    '-outfmt', str(self.outfmt),
                    '-query', '-',
                    '-out', self.outpath,
                    '-num_threads', str(self.num_threads)])
        # will raise ValueError when run, due to not having a recognizable
        # input type
        else:
            return [], self.outpath
        
        return cmd, self.outpath

'''
A WITCH alignment job that will run WITCH to align a set of query reads
to their target marker gene backbone alignment
'''
class WITCHAlignmentJob(Job):
    def __init__(self, **kwargs):
        Job.__init__(self)
        self.job_type = 'witch-alignment'

        # initialize parameters
        self.path = ''
        self.query_path = ''
        self.backbone_path = ''
        self.backbone_tree_path = ''
        self.outdir = ''
        self.num_cpus = 1

        # set parameters
        for k, v in kwargs.items():
            setattr(self, k, v)
        self.kwargs = kwargs
    
    def get_invocation(self):
        self.outpath = os.path.join(self.outdir, 'est.aln.masked.fasta')
        cmd = [self.path, '-o', 'est.aln.fasta',
                ]
        # extend from additional kwargs
        for k, v in self.kwargs.items():
            if k != 'path':
                param_name = k.replace('_', '-')
                cmd.extend([f"--{param_name}", str(v)])
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
        self.backbone_alignment_path = ''
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
                '-a', self.backbone_alignment_path,
                '-q', self.query_alignment_path,
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

'''
A generic job that is designed for univeral types of additional jobs that
a user/coder can add to TIPP3. just need to supplement the correct configuration
in main.config or user.config
NOTE: the binary executable path given needs to be executable
'''
class GenericJob(Job):
    def __init__(self, **kwargs):
        Job.__init__(self)
        self.job_type = 'generic'

        # note: [job_type, path, outpath] needs to be set from configuration file
        for k, v in kwargs.items():
            setattr(self, k, v)
    
    def get_invocation(self):
        outpath, cmd = None, None 
        if getattr(self, 'outpath') != None:
            outpath = self.outpath
        if getattr(self, 'cmd') != None:
            cmd = cmd.split()
        
        # raise Error for generic job that does not have a command 
        if not outpath or not cmd:
            _LOG.error(f"Generic Job {self.job_type} " 
                "do not have outpath and cmd in configuration.")
            exit(1)

        return cmd, outpath

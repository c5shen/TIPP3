"""
A collection of jobs for TIPP3.

Jobs are designed to be run standalone, as long as all parameters
are provided correctly.
"""

import os, shutil, subprocess, traceback
from subprocess import Popen
from abc import abstractmethod

from tipp3 import get_logger

_LOG = get_logger(__name__)


def _validate_executable(binpath, job_type):
    """Validate that a binary path exists and is executable."""
    if shutil.which(binpath) is not None:
        return
    if os.path.isfile(binpath) and os.access(binpath, os.X_OK):
        return
    if binpath == 'java' or binpath.endswith('.jar'):
        return
    raise FileNotFoundError(
        f"Executable for {job_type} not found: '{binpath}'. "
        "Please check your main.config or ensure the tool is installed.")


def _run_pipeline(cmd_groups, stdin_data="", logging_fobj=None):
    """Run a pipeline of commands connected by pipes.

    Args:
        cmd_groups: list of lists, e.g. [['gzip', '-dc', 'f.gz'], ['blastn', ...]]
        stdin_data: string data to feed to the first process's stdin
        logging_fobj: open file object to redirect final stdout, or None

    Returns:
        (stdout, stderr, returncode) of the last process in the pipeline
    """
    if not cmd_groups:
        raise ValueError("Empty command pipeline")

    processes = []
    prev_stdout = subprocess.PIPE

    for i, cmd in enumerate(cmd_groups):
        is_first = (i == 0)
        is_last = (i == len(cmd_groups) - 1)

        stdin_arg = subprocess.PIPE if is_first else processes[-1].stdout
        if is_last and logging_fobj is not None:
            stdout_arg = logging_fobj
        else:
            stdout_arg = subprocess.PIPE

        p = Popen(cmd, text=True, bufsize=1,
                  stdin=stdin_arg, stdout=stdout_arg,
                  stderr=subprocess.PIPE)
        # allow upstream process to receive SIGPIPE
        if not is_first and processes[-1].stdout:
            processes[-1].stdout.close()
        processes.append(p)

    # send stdin data to the first process
    if stdin_data and processes:
        try:
            processes[0].stdin.write(stdin_data)
        except BrokenPipeError:
            pass
        processes[0].stdin.close()

    # wait for all processes to complete
    for p in processes:
        p.wait()

    last = processes[-1]
    stderr = last.stderr.read() if last.stderr else ''
    stdout = ''
    if last.stdout and logging_fobj is None:
        stdout = last.stdout.read()

    return stdout, stderr, last.returncode


class Job(object):
    """Template class for running external jobs."""

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

    def run(self, stdin="", lock=None, logging=False, shell=False):
        """Run the job with the invocation defined in a child class."""
        try:
            cmd, outpath = self.get_invocation()
            _LOG.debug(f"Running job_type: {self.job_type}, output: {outpath}")

            if len(cmd) == 0:
                raise ValueError(
                    f"{self.job_type} does not have a valid run command. "
                    "Check that your input file format is supported "
                    "(.fa/.fasta/.fq/.fastq, optionally .gz compressed).")

            # cmd can be a flat list (single command) or list-of-lists (pipeline)
            if isinstance(cmd[0], list):
                cmd_groups = cmd
            else:
                cmd_groups = [cmd]

            # validate the primary executable of each command in the pipeline
            for group in cmd_groups:
                binpath = group[0]
                if binpath == 'java' and len(group) > 2:
                    _validate_executable(group[2], self.job_type)
                else:
                    _validate_executable(binpath, self.job_type)

            _LOG.debug("Command pipeline: %s",
                       " | ".join(" ".join(str(x) for x in g) for g in cmd_groups))

            log_fobj = None
            try:
                if logging:
                    logpath = os.path.join(
                        os.path.dirname(outpath), f'{self.job_type}.txt')
                    log_fobj = open(logpath, 'w', 1)

                stdout, stderr, self.returncode = _run_pipeline(
                    cmd_groups, stdin_data=stdin, logging_fobj=log_fobj)
            finally:
                if log_fobj is not None:
                    log_fobj.close()

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
                error_msg = (f"Error running {self.job_type}. "
                             f"Return code: {self.returncode}")
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
                raise RuntimeError(error_msg + '\nSTDERR: ' + stderr)
        except Exception:
            _LOG.error(traceback.format_exc())
            raise

    @abstractmethod
    def get_invocation(self):
        """Return (cmd, outpath). cmd is a list of args for a single command,
        or a list of lists for a pipeline of commands connected by pipes."""
        raise NotImplementedError(
            "get_invocation() should be implemented by subclasses")


class BlastnJob(Job):
    """Run BLASTN to bin reads against the reference package marker genes."""

    def __init__(self, **kwargs):
        Job.__init__(self)
        self.job_type = 'blastn'
        self.outfmt = 0
        self.path = ''
        self.query_path = ''
        self.database_path = ''
        self.outdir = ''
        self.num_threads = 1

        for k, v in kwargs.items():
            setattr(self, k, v)

    @staticmethod
    def _awk_fastq_to_fasta_cmd():
        """Return an awk command (as arg list) that converts FASTQ to FASTA."""
        return ['awk', 'NR%4==1 {print ">"substr($0, 2)} NR%4==2 {print $0}']

    def _blastn_cmd(self, query_source='-'):
        """Return the core blastn command as a list."""
        return [self.path, '-db', self.database_path,
                '-outfmt', str(self.outfmt),
                '-query', query_source,
                '-out', self.outpath,
                '-num_threads', str(self.num_threads)]

    def get_invocation(self):
        self.outpath = os.path.join(self.outdir, 'blast.alignment.out')

        name_parts = self.query_path.split('.')
        suffix = name_parts[-1].lower()

        # FASTA files - direct input
        if suffix in ('fa', 'fasta'):
            cmd = self._blastn_cmd(query_source=self.query_path)
            return cmd, self.outpath

        # FASTQ files - convert to FASTA via awk, then pipe to blastn
        if suffix in ('fq', 'fastq'):
            pipeline = [
                self._awk_fastq_to_fasta_cmd() + [self.query_path],
                self._blastn_cmd(query_source='-'),
            ]
            return pipeline, self.outpath

        # Gzipped files - decompress, optionally convert FASTQ, then blastn
        if suffix in ('gz', 'gzip'):
            pipeline = [['gzip', '-dc', self.query_path]]

            suffix2 = name_parts[-2].lower() if len(name_parts) > 2 else ''
            if suffix2 in ('fastq', 'fq'):
                pipeline.append(self._awk_fastq_to_fasta_cmd())

            pipeline.append(self._blastn_cmd(query_source='-'))
            return pipeline, self.outpath

        _LOG.warning(f"Unrecognized input format: '{suffix}' from {self.query_path}")
        return [], self.outpath

class WITCHAlignmentJob(Job):
    """Run WITCH to align query reads to a marker gene backbone alignment."""

    def __init__(self, **kwargs):
        Job.__init__(self)
        self.job_type = 'witch-alignment'
        self.path = ''
        self.query_path = ''
        self.backbone_path = ''
        self.backbone_tree_path = ''
        self.outdir = ''
        self.num_cpus = 1

        for k, v in kwargs.items():
            setattr(self, k, v)
        self.kwargs = kwargs

    def get_invocation(self):
        self.outpath = os.path.join(self.outdir, 'est.aln.masked.fasta')
        cmd = [self.path, '-o', 'est.aln.fasta']
        for k, v in self.kwargs.items():
            if k != 'path':
                param_name = k.replace('_', '-')
                cmd.extend([f"--{param_name}", str(v)])
        return cmd, self.outpath

class BscamppJob(Job):
    """Run BSCAMPP for placing aligned query reads."""
    def __init__(self,
            path, query_alignment_path, backbone_alignment_path,
            backbone_tree_path, tree_model_path, outdir,
            base_method, num_cpus, **kwargs):
        Job.__init__(self)
        self.job_type = 'bscampp'
        
        # initialize parameters
        self.path = path
        self.query_alignment_path = query_alignment_path
        self.backbone_alignment_path = backbone_alignment_path
        self.backbone_tree_path = backbone_tree_path
        self.tree_model_path = tree_model_path
        self.outdir = outdir
        self.num_cpus = num_cpus
        self.kwargs = kwargs

        if base_method is not None:
            self.kwargs['placement_method'] = base_method

    def get_invocation(self):
        self.outpath = os.path.join(self.outdir, 'placement.jplace')
        cmd = [self.path,
               '-q', self.query_alignment_path,
               '-a', self.backbone_alignment_path,
               '-t', self.backbone_tree_path,
               '-i', self.tree_model_path,
               '-d', self.outdir,
               '-o', 'placement',
               '--num-cpus', str(self.num_cpus)]
        for k, v in self.kwargs.items():
            if k == 'support_value':
                continue
            param = k.replace('_', '-')
            cmd.extend([f'--{param}', str(v)])
        return cmd, self.outpath


class PplacerTaxtasticJob(Job):
    """Run pplacer with the taxtastic refpkg."""
    def __init__(self, path, query_alignment_path, refpkg_path, outdir,
                 num_cpus, **kwargs):
        Job.__init__(self)
        self.job_type = 'pplacer-taxtastic'
        self.path = path
        self.query_alignment_path = query_alignment_path
        self.refpkg_path = refpkg_path
        self.outdir = outdir
        self.num_cpus = num_cpus
        self.model_type = 'GTR'
        self.kwargs = kwargs

    def get_invocation(self):
        self.outpath = os.path.join(self.outdir, 'placement.jplace')
        cmd = [self.path,
               '-m', self.model_type,
               '-c', self.refpkg_path,
               '-o', self.outpath,
               '-j', str(self.num_cpus),
               self.query_alignment_path]
        return cmd, self.outpath



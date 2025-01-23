'''
Aligning query reads to their assigned marker genes, using either BLASTN
output or some other alignment method
'''

import os, time, shutil
from tipp3 import get_logger
from tipp3.configs import Configs
from tipp3.helpers.alignment_tools import Alignment
from tipp3.jobs import WITCHAlignmentJob

_LOG = get_logger(__name__)

'''
Align query reads to marker genes with the user-defined alignment methods
'''
def queryAlignment(refpkg, query_paths):
    # query_alignments only retain query alignments
    # the backbone part is unnecessary
    _LOG.info(f"Started alignment method: {Configs.alignment_method}")

    query_alignment_paths = {}
    additional_kwargs = {}
    try:
        additional_kwargs = getattr(Configs, Configs.alignment_method).__dict__
    except (AttributeError, ValueError) as e:
        pass
    _LOG.info("(From main.config) parameter options for " \
            f"{Configs.alignment_method}: {additional_kwargs}")

    # removing temporary files that we do not care
    temp_files = ['est.aln.fasta', 'tree_decomp', 'debug.txt']

    idx = 0; total = len(query_paths)
    for marker, query_path in query_paths.items():
        _LOG.info(f"Aligning query reads for {marker}")
        idx += 1
        alignment_dir = os.path.join(Configs.outdir, 'query_alignments',
                marker)
        if not os.path.isdir(alignment_dir):
            os.makedirs(alignment_dir)

        # multi-level detect to avoid rerunning:
        # first detect if 'est.aln.fasta' exists
        d1_path = os.path.join(alignment_dir, 'est.aln.masked.fasta')
        d2_path = os.path.join(alignment_dir, 'est.aln.masked.fasta.lz4')
        skip_alignment = False
        if os.path.exists(d1_path) and os.stat(d1_path).st_size > 0:
            # lz4 the file and remove d1_path
            _LOG.info(f"Found existing alignment: {d1_path}, compressing...")
            os.system(f"lz4 --rm -f -q {d1_path} {d2_path}")
            skip_alignment = True
        elif os.path.exists(d2_path) and os.stat(d2_path).st_size > 0:
            # keep the lz4 file and avoid re-running the alignment step
            _LOG.info(f"Found existing compressed alignment: {d2_path}")
            skip_alignment = True

        # not skipping alignment and need to generate the alignment file
        # to extract query alignments
        if not skip_alignment:
            # obtain the backbone alignment and tree paths for alignment
            backbone_path = refpkg[marker]['alignment']
            backbone_tree_path = refpkg[marker]['alignment-decomposition-tree']

            # identify job type - different alignment method(s)
            if Configs.alignment_method == 'witch':
                alignment_job = WITCHAlignmentJob(path=Configs.witch_path,
                        query_path=query_path,
                        backbone_path=backbone_path,
                        backbone_tree_path=backbone_tree_path,
                        outdir=alignment_dir,
                        num_cpus=Configs.num_cpus,
                        **additional_kwargs)
            else:
                raise NotImplementedError(
                        f"Alignment method {Configs.alignment_method} " \
                        " is not implemented yet.")
            raw_alignment_path = alignment_job.run(logging=True)
            # zip the output masked file
            os.system(f"lz4 --rm -f -q {raw_alignment_path} {d2_path}")

        # extract query alignment from {d2_path}
        queryunaln = Alignment(); queryunaln.read_file_object(query_path)
        query_keys = list(queryunaln.keys())
        aln = Alignment(); aln.read_file_object(d2_path)
        subaln = aln.sub_alignment(query_keys)
        # write to local
        query_alignment_path = os.path.join(alignment_dir,
                'est.aln.masked.queries.fasta')
        subaln.write(query_alignment_path, 'FASTA')
        query_alignment_paths[marker] = query_alignment_path

        # clean up some temp files
        for tf in temp_files:
            toremove = os.path.join(alignment_dir, tf)
            if not os.path.exists(toremove):
                continue
            if os.path.isfile(toremove):
                os.remove(toremove)
            else:
                shutil.rmtree(toremove)
        _LOG.info(f"({idx}/{total}) Finished aligning query reads for {marker}")
    return query_alignment_paths

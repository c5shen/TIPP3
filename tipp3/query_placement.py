'''
Place reads to marker gene taxonomic tree using a placement method
'''

import os, time, shutil
from tipp3 import get_logger
from tipp3.configs import Configs
from tipp3.jobs import BscamppJob, PplacerTaxtasticJob

_LOG = get_logger(__name__)


'''
Place query reads to their assigned marker gene taxonomic trees using the
user-defined phylogenetic placement method
'''
def queryPlacement(refpkg, query_alignment_paths):
    _LOG.info(f"Started placement method: {Configs.placement_method}")

    query_placement_paths = {}
    additional_kwargs = {}
    try:
        #additional_kwargs = Configs.BSCAMPP.__dict__
        additional_kwargs = getattr(Configs, Configs.placement_method).__dict__
    except (AttributeError, ValueError) as e:
        pass
    _LOG.info("(From main.config) parameter options for " \
            f"{Configs.placement_method}: {additional_kwargs}")

    # failsafe when no alignments are given
    if len(query_alignment_paths) == 0:
        _LOG.warning("No alignments are given for placements, returning...")
        return query_placement_paths
    
    # removing temporary files that we do not care
    temp_folders = ['tmp0']

    idx = 0; total = len(query_alignment_paths) 
    for marker, query_alignment_path in query_alignment_paths.items():
        _LOG.info("Placing aligned query reads for " \
                f"{marker}") #from {query_alignment_path}")
        idx += 1
        placement_dir = os.path.join(Configs.outdir, 'query_placements',
                marker)
        if not os.path.isdir(placement_dir):
            os.makedirs(placement_dir)

        # detect results from previous run, if exists return
        detect_path = os.path.join(placement_dir, 'placement.jplace')
        if os.path.exists(detect_path) and os.stat(detect_path).st_size > 0:
            _LOG.info(f"Found existing placement: {detect_path}")
            query_placement_paths[marker] = detect_path
            continue

        # identify job type - different placement methods
        # pplacer-taxtastic
        if Configs.placement_method == 'pplacer-taxtastic':
            placement_job = PplacerTaxtasticJob(path=Configs.pplacer_path,
                    query_alignment_path=query_alignment_path,
                    refpkg_path=refpkg[marker]['path'],
                    outdir=placement_dir,
                    num_cpus=Configs.num_cpus,
                    **additional_kwargs)

        # default: bscampp Job
        else:
            backbone_alignment_path = refpkg[marker]['alignment']

            placement_job = BscamppJob(path=Configs.bscampp_path,
                    query_alignment_path=query_alignment_path,
                    backbone_alignment_path=backbone_alignment_path,
                    backbone_tree_path=refpkg[marker]['additional-raxml-br-tree'],
                    tree_model_path=refpkg[marker]['additional-raxml-model-file'],
                    outdir=placement_dir,
                    num_cpus=Configs.num_cpus,
                    **additional_kwargs) 
        placement_path = placement_job.run(logging=True)
        query_placement_paths[marker] = placement_path

        # clean up some temp files
        for tf in temp_folders:
            toremove = os.path.join(placement_dir, tf)
            if os.path.exists(toremove):
                shutil.rmtree(toremove)
        _LOG.info(f"({idx}/{total}) Finished placing aligned query reads for {marker}")
    return query_placement_paths

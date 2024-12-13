'''
Aligning query reads to their assigned marker genes, using either BLASTN
output or some other alignment method
'''

import os, time
from tipp3 import get_logger
from tipp3.configs import Configs

_LOG = get_logger(__name__)

'''
Align query reads to marker genes with the user-defined alignment methods
'''
def queryAlignment(refpkg, query_paths):
    # query_alignments only retain query alignments
    # the backbone part is unnecessary
    query_alignment_paths = defaultdict(str)

    return query_alignment_paths

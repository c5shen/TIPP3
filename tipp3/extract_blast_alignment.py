'''
Extract and return the blast alignment paths
'''

import os, sys, re
from collections import defaultdict

from tipp3 import get_logger
from tipp3.helpers.alignment_tools import ExtendedAlignment

_LOG = get_logger(__name__)

global complement_map
complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

'''
Reverse-complement given sequences. Retain characters that are not {A, T, C, G}
'''
def reverseComplement(seq):
    try:
        char_list = [complement_map[c] for c in seq]
    except KeyError:
        char_list = [complement_map[c] if c in complement_map else c for c in seq]
    new_seq = ''.join(char_list)
    return new_seq[::-1]

'''
Given an aligned string (represented by upper/lower case letters and gaps,
return a condensed version that have the lowercase letters from both sides
compressed to front/back.
'''
def compressInsertions(seq):
    p = re.compile(r'[A-Z]+')
    alns = [(m.start(), m.end()) for m in p.finditer(seq)]
    # do not perform such task if there is no aligned column at all
    if len(alns) == 0:
        return seq

    # first occurrence of aligned position defines the back of front search space
    # i.e., [start, end)
    f_start, f_end = 0, alns[0][0]
    f_len = f_end - f_start

    # last occurrence of aligned position defines the front of the back search space
    b_start, b_end = alns[-1][1], len(seq)
    b_len = b_end - b_start

    # simplest way of compression: remove all gaps and add them back
    f_str_ins = seq[f_start:f_end].replace('-', '')
    f_len_ins = len(f_str_ins)
    f_str = f_str_ins + '-' * (f_len - f_len_ins)

    b_str_ins = seq[b_start:b_end].replace('-', '')
    b_len_ins = len(b_str_ins)
    b_str = '-' * (b_len - b_len_ins) + b_str_ins

    # combine the compressed front/back with the remaining sequence
    combined = f_str + seq[f_end:b_start] + b_str

    return combined

'''
Obtain retained columns of a target sequence
'''
def getRetainedColumns(aln, target_keys):
    target_to_retained_columns = {}
    aln_len = aln.sequence_length()
    
    # we only care about where each character appears in the backbone columns
    for target in target_keys:
        retained_columns = defaultdict(int)
        idx = 0
        for i in range(aln_len):
            if aln[target][i] != '-':
                retained_columns[idx] = i
                idx += 1
        target_to_retained_columns[target] = retained_columns
    return target_to_retained_columns

'''
Obtain aligned columns of query seq to target
'''
def getAlignedColumns(taxon, seq_len, aln_seq, target_seq, target_start):
    # obtain aligned columns of the query sequence (to the target sequence)
    qidx = 0; tidx = target_start
    aligned_columns = [-1] * seq_len
    for i in range(len(aln_seq)):
        # if current is a deletion, increment tidx
        if aln_seq[i] == '-':
            tidx += 1
            continue
        else:
            # current is an insertion, do not increment tidx and map to -1
            if target_seq[i] == '-':
                aligned_columns[qidx] = -1
                qidx += 1
            # current is a match, increment both indexes 
            else:
                # current 
                aligned_columns[qidx] = tidx
                qidx += 1
                tidx += 1
    return aligned_columns

#'''
#Obtain alignment_graph and backtrace to find the query merged alignment
#'''
#def getBacktrace(aligned_columns, retained_columns):
#    # target start marks from where blast considers the alignment of the target sequence
#    alignment_graph = [[0 for _j in range(backbone_length + 1)] for _i in range(seq_len + 1)]
#    weights = defaultdict(int)
#    min_col_ind, max_col_ind = backbone_length + 1, -1
#
#    for i in range(len(aligned_columns)):
#        subset_col = aligned_columns[i]
#        # ignore insertions
#        if subset_col == -1:
#            continue
#        j = retained_columns[subset_col]
#        weights[(i, j)] = 1   # meaning query[i] is aligned to backbone position j
#        min_col_ind = min(min_col_ind, j)
#        max_col_ind = max(max_col_ind, j)
#    # print(weights)
#
#    backtrace = [[0 for _j in range(backbone_length + 1)] for _i in range(seq_len + 1)]
#    for i in range(0, seq_len + 1):
#        for j in range(min_col_ind, max_col_ind + 2):
#            if i == 0 or j == min_col_ind:
#                alignment_graph[i][j] = 0   # on the edge of the alignment graph
#                continue
#            cur_max = 0; cur_bt = 0
#            w = weights.get((i-1, j-1), 0)
#            values = [alignment_graph[i-1][j-1] + w,
#                      alignment_graph[i-1][j],
#                      alignment_graph[i][j-1]]
#            for _ind, val in enumerate(values):
#                if _ind == 0 and w <= 0:
#                    cur_bt = 1    # in the case position i-1 in query is not mapped to j-1 in backbone, prefer going "upward"
#                    continue
#                if val > cur_max:
#                    cur_max = val
#                    cur_bt = _ind
#            alignment_graph[i][j] = cur_max
#            backtrace[i][j] = cur_bt
#    return backtrace, min_col_ind, max_col_ind

#'''
#Get combined aligned query sequence from backtrace
#'''
#def getCombinedSequence(seq, seq_len, backtrace, min_col_ind, max_col_ind):
#    # retrieve results from backtrace
#    result = []; i, j = seq_len, max_col_ind + 1
#    while i > 0 and j > min_col_ind:
#        bt = backtrace[i][j]
#        if bt == 0:
#            result.append(seq[i-1].upper())
#            i -= 1; j -= 1
#        elif bt == 1:
#            # query "insertion", lowercase
#            result.append(seq[i-1].lower())
#            i -= 1
#        elif bt == 2:
#            # query "deletion", gap
#            result.append('-')
#            j -= 1
#        else:
#            raise ValueError
#    while i > 0:
#        result.append(seq[i-1].lower()); i -= 1
#    while j > min_col_ind:
#        result.append('-'); j -= 1
#
#    # reverse result to be in correct order
#    result = result[::-1]
#
#    # append '-' to start and end of results to fill up the alignment
#    result = ['-'] * min_col_ind + result + ['-'] * (backbone_length - max_col_ind - 1)
#
#    # result
#    combined = compressInsertions(''.join(result))
#    return combined

'''
Convert combined sequence to an Extended Alignment object
'''
def getQueryExtendedAlignment(taxon, combined):
    # get an extended alignment object with just the query sequence and its updatd indexes
    query = ExtendedAlignment([])
    query[taxon] = combined; query._reset_col_names()
    insertion = -1; regular = 0
    for i in range(len(combined)):
        if combined[i].islower():
            query._col_labels[i] = insertion; insertion -= 1
        else:
            query._col_labels[i] = regular; regular += 1
    return query

'''
Directly construct the query alignment to backbone by setting characters based on retained_columns
'''
def constructCombined(seq, backbone_length, aligned_columns, retained_columns):
    ret = ['-'] * backbone_length   # initialize as gaps, later modify to be characters from seq
    for i in range(len(aligned_columns)):
        target_col = aligned_columns[i]
        if target_col == -1:
            continue
        # assign seq[i] to ret[retained_columns[target_col]]
        ret[retained_columns[target_col]] = seq[i]
    return ''.join(ret)

def extractionRunner(marker, inpath, bbpath, outpath):
    # get mapped queries to marker gene
    aln = ExtendedAlignment([])
    aln.read_file_object(bbpath)
    backbone_length = aln.sequence_length()

    mapped = defaultdict(list)
    with open(inpath, 'r') as f:
        lines = f.read().strip().split('\n')
        for i in range(0, len(lines), 2):
            taxon = lines[i].split('>')[1]
            items = eval(lines[i+1])
            mapped[taxon] = items

    # compute retained columns just for once for the target marker gene sequences that appear
    target_keys = set()
    for taxon in mapped.keys():
        target_keys.add(mapped[taxon][0][0])
    target_to_retained_columns = getRetainedColumns(aln, target_keys)

    #for taxon in tqdm(mapped.keys(), total=len(mapped), **tqdm_styles):
    for taxon in mapped.keys():
        aln_seq = mapped[taxon][3]
        target = mapped[taxon][0][0]
        target_seq = mapped[taxon][6]

        seq = aln_seq.replace('-', '')
        seq_len = len(seq)
        # deal with reverse-complement
        # target sequence is reverse-complemented
        #   - we need to restore that
        #   - and make the query reverse-complemented
        if mapped[taxon][4] > mapped[taxon][5]:
            target_start = mapped[taxon][5] - 1
            aln_seq = reverseComplement(aln_seq)
            target_seq = reverseComplement(target_seq)
            seq = reverseComplement(seq)
        else:
            target_start = mapped[taxon][4] - 1

        # obtain aligned columns of taxon
        aligned_columns = getAlignedColumns(taxon, seq_len, aln_seq, 
                target_seq, target_start)
        # obtain retained columns of target
        retained_columns = target_to_retained_columns[target]
        # directly construct the aligned query sequence by adding gaps where necessary and ignoring insertion characters
        # this will be the masked version
        combined = constructCombined(seq, backbone_length, aligned_columns, retained_columns)
        aln[taxon] = combined
    # write a masked version to local (masking lower case letters)
    # NOTE: only retaining query keys to save space
    subaln = aln.sub_alignment(list(mapped.keys()))
    subaln.write(outpath, 'FASTA')

    ## write a query-only version
    #subaln = aln.sub_alignment(list(mapped.keys()))
    #subaln.write(outqueriespath, 'FASTA')

    ## write a backbone-only version
    #bbset = set(aln.keys()).difference(set(mapped.keys()))
    #bbaln = aln.sub_alignment(list(bbset))
    #bbaln.write(outbackbonepath, 'FASTA')

'''
entry point for extracting BLAST alignment
'''
def extractBlastAlignment(refpkg, workdir, query_blast_paths):
    _LOG.info("Started extracting query pairwise alignments from BLAST output") 

    # each marker's BLAST result path is determined in query_blast_paths
    alignment_outdir = workdir + '/query_alignments'
    if len(query_blast_paths) == 0:
        _LOG.warning('No query reads are assigned to any marker '
        'gene according BLAST.')
        return None

    ret = {}
    for marker, inpath in query_blast_paths.items():
        bbpath = refpkg[marker]['alignment']
        outdir = os.path.join(alignment_outdir, marker)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        # NOTE: only keep the query alignments (backbone part is unnecessary)
        outpath = os.path.join(outdir, 'est.aln.masked.queries.fasta')
        ret[marker] = outpath

        if os.path.exists(outpath) and os.stat(outpath).st_size > 0:
            _LOG.info(f"Skipping {marker}, found existing alignment: {outpath}")
            continue

        #outqueriespath = os.path.join(outdir, 'est.aln.masked.queries.fasta')
        #outbackbonepath = os.path.join(outdir, 'est.aln.masked.backbone.fasta')
        extractionRunner(marker, inpath, bbpath, outpath)
        _LOG.debug('Writing BLAST alignment for marker gene {} to {}'.format(
            marker, outpath))
    return ret

    #inpath = sys.argv[1]    # categorized by marker gene
    #bbpath = sys.argv[2]    # marker gene backbone alignment
    #outpath = sys.argv[3]   # output with query alignments merged
    #outqueriespath = sys.argv[4]    # output with only queries
    #outbackbonepath = sys.argv[5]    # output with only backbone sequences
    #runner(inpath, bbpath, outpath, outqueriespath, outbackbonepath)


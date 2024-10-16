'''
Read binning using BLASTN
'''

import os, time
from collections import defaultdict

from tipp3 import get_logger
from tipp3.configs import Configs
from tipp3.jobs import BlastnJob
from tipp3.extract_blast_alignment import extractBlastAlignment

_LOG = get_logger(__name__)

# simple reverse complement map 
global complement_map
complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

'''
Map reads to refpkg using BLASTN, then process the reads to extract
their assignments and/or alignments to marker genes
'''
def readBinning(refpkg):
    database_path = refpkg['blast']['database']
    _LOG.info(f"Initializing BlastnJob with BLASTN database: {database_path}")

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
    query_blast_paths, query_aln = processBlastnOutput(refpkg, blastn_outpath,
            blastn_outdir)

    # write query reads to local for alignments, truncate if neccessary
    query_outdir = os.path.join(Configs.outdir, 'queries')
    query_paths = splitQueries(refpkg, query_aln, query_outdir)

    # if Configs.alignment_method is BLAST, extract alignment at this step
    query_alignment_paths = {}
    if Configs.alignment_method == 'blast':
        query_alignment_paths = extractBlastAlignment(refpkg,
                Configs.outdir, query_blast_paths)

    return query_paths, query_alignment_paths

'''
process the BLASTN output, extract assignment information and/or alignment
'''
def processBlastnOutput(refpkg, blastn_outpath, blastn_outdir):
    # seq-to-marker map
    map_path = refpkg['blast']['seq-to-marker-map']
    gene_mapping = readGeneMapping(map_path)

    # read blastn output to memory
    threshold = 50
    try:
        threshold = int(Configs.BLAST.threshold)
    except (AttributeError, ValueError) as e:
        pass
    _LOG.info(f"Filtering BLASTN results: >= {threshold}bp coverage.")

    query_aln = readBlastnOutput(blastn_outpath, gene_mapping, threshold)
    query_blast_paths = {}
    
    # read in marker genes to use
    markers = refpkg['genes']

    # write to blastn_outdir as checkpoints
    write_order = ['source_taxon', 'qstart', 'qend', 'qaln', 'sstart', 'send', 'saln']
    marker_fptr = {}
    for marker in markers:
        _path = '{}/queries-{}.out'.format(blastn_outdir, marker) 
        marker_fptr[marker] = {}
        marker_fptr[marker]['path'] = _path
        marker_fptr[marker]['fptr'] = open(_path, 'w')
        marker_fptr[marker]['count'] = 0
    
    for k, v in query_aln.items():
        marker = v['source_taxon'][1]
        # only deal with 
        if marker in marker_fptr:
            towrite = [v[_x] for _x in write_order]
            marker_fptr[marker]['fptr'].write('>{}\n{}\n'.format(k, towrite))
            marker_fptr[marker]['count'] += 1

    for marker in markers:
        marker_fptr[marker]['fptr'].close()
        _path = marker_fptr[marker]['path']
        # remove empty marker gene outputs
        if marker_fptr[marker]['count'] == 0:
            _LOG.debug(f'{marker} has no assigned queries, removing {_path}')
            os.system('rm {}'.format(_path))
        else:
            query_blast_paths[marker] = marker_fptr[marker]['path']
    return query_blast_paths, query_aln


######################## HELPER FUNCTIONS ###########################

'''
helper function to obtain the reverse complement of a sequence
'''
def reverseComplement(seq):
    try:
        char_list = [complement_map[c] for c in seq]
    except KeyError:
        char_list = [complement_map[c] if c in complement_map else c for c in seq]
    new_seq = ''.join(char_list)
    return new_seq[::-1]

'''
helper function to read marker gene mapping
'''
def readGeneMapping(map_path, delimiter='\t'):
    mapping = {}
    with open(map_path, 'r') as f:
        for line in f:
            parts = line.strip().split(delimiter)
            mapping[parts[0]] = parts
    return mapping

'''
helper function to read in blastn output
'''
def readBlastnOutput(blastn_outpath, gene_mapping, threshold):
    query_aln = defaultdict(dict)
    with open(blastn_outpath, 'r') as f:
        line = f.readline()
        b_consider, scan_mode, cur_taxon = False, False, ''

        # initialization of the following
        # source_taxon, qcov, qstart, qend, qaln, sstart, send, saln
        cur_aln = {'source_taxon': ('', ''), 'qcov': -1,
                'qstart': -1, 'qend': -1, 'qaln': '',
                'sstart': -1, 'send': -1, 'saln': ''}

        while line:
            line = line.split('\n')[0].strip()
            # start of a new query entry, starting reading
            if line.startswith('Query='):
                taxon = line.split('= ')[1]
                b_consider = True
                cur_taxon = taxon
            # reach the end of a query entry
            elif line.startswith('Effective search space'):
                # check if (1) qcov >= threshold
                #          (2) cur_aln has better coverage than previous one
                if cur_taxon != '':
                    cur_aln = updateQueryAlignment(
                            query_aln, cur_taxon, cur_aln, threshold)
                    cur_taxon = ''
                    scan_mode = False
            # is actively reading in entries for a query
            elif b_consider:
                if line.startswith('>'):
                    # cur_aln is not empty, flush current results before moving
                    # to the next one
                    if cur_aln['source_taxon'][1] != '':
                        cur_aln = updateQueryAlignment(
                                query_aln, cur_taxon, cur_aln, threshold)
                    source = line.split('>')[1]
                    # obtain corresponding marker gene from gene_mapping
                    scan_mode = True
                    cur_aln['source_taxon'] = (source, gene_mapping[source][1])

                    # skip scanning since we just start, move to next line
                    line = f.readline()
                    continue
                
                # is actively scanning for alignment part
                if scan_mode:
                    # query alignment part
                    if line.startswith('Query'):
                        # parts[1]: current qstart
                        # parts[2]: current qaln
                        # parts[3]: current qend
                        parts = line.split()
                        if cur_aln['qstart'] == -1:
                            cur_aln['qstart'] = int(parts[1])
                        cur_aln['qaln'] += parts[2]
                        cur_aln['qend'] = int(parts[3])
                    # target alignment part
                    elif line.startswith('Sbjct'):
                        parts = line.split()
                        if cur_aln['sstart'] == -1:
                            cur_aln['sstart'] = int(parts[1])
                        cur_aln['saln'] += parts[2]
                        cur_aln['send'] = int(parts[3])
            # move to next line
            line = f.readline()
    return query_aln

'''
helper function to update query_aln with a new entry recorded
replace the old one if the new one has longer coverage
returns a completely new entry for the following scan
'''
def updateQueryAlignment(query_aln, cur_taxon, cur_aln, threshold):
    # compute current query coverage
    cur_aln['qcov'] = abs(cur_aln['qstart'] - cur_aln['qend']) + 1
    if cur_aln['qcov'] >= threshold:
        # new entry
        if len(query_aln[cur_taxon]) == 0:
            query_aln[cur_taxon] = cur_aln
            # compare current qcov to the one already stored
        elif cur_aln['qcov'] > query_aln[cur_taxon]['qcov']:
            query_aln[cur_taxon] = cur_aln
    return {'source_taxon': ('', ''), 'qcov': -1,
        'qstart': -1, 'qend': -1, 'qaln': '',
        'sstart': -1, 'send': -1, 'saln': ''}

'''
helper function to slit queries to marker genes and write to local for
alignment later
'''
def splitQueries(refpkg, query_aln, query_outdir):
    markers = refpkg['genes']

    if not os.path.isdir(query_outdir):
        os.makedirs(query_outdir)

    # initialize assigned queries for each marker gene
    ret = {}
    marker_fptr = {}
    for marker in markers:
        _outpath = os.path.join(query_outdir, '{}.queries.fasta'.format(marker))
        ret[marker] = _outpath
        marker_fptr[marker] = open(_outpath, 'w')

    # write blastn alignment output info
    info_path = os.path.join(query_outdir, 'blast-binned.out')
    info_f = open(info_path, 'w', buffering=1)
    info_f.write('qseqid,sseqid,marker,trim_qstart,trim_qend,qlen\n')

    # go over each entry and write query to corresponding marker gene output
    for taxon, v in query_aln.items():
        sseqid, smarker = v['source_taxon']
        qstart, qend, sstart, send = v['qstart'], v['qend'], v['sstart'], v['send']

        # retain only the part that's mapped by BLASTN
        seq = v['qaln'].replace('-', '').upper()

        # reverse complement
        if sstart > send:
            rev_seq = reverseComplement(seq)
        else:
            rev_seq = seq
        
        # write to corresponding marker gene
        marker_fptr[smarker].write('>{}\n{}\n'.format(taxon, rev_seq))
        info_f.write('{},{},{},{},{},{}\n'.format(taxon, sseqid, smarker,
            qstart, qend, qend - qstart + 1))

    for marker in markers:
        marker_fptr[marker].close()
    info_f.close()

    # remove empty queries.fasta
    for marker, _outpath in ret.items():
        if os.stat(_outpath).st_size == 0:
            _LOG.debug(f'Removing redundant {marker}: {_outpath}')
            os.remove(_outpath)
    return ret


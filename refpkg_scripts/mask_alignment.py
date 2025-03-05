#!/usr/bin/env python

import os, sys, gzip
from operator import add

def read_fasta(path, suffix=None):
    ret = {}
    name = None
    seq_list = list()
    if not suffix:
        f = open(path, 'r')
    elif suffix == 'gz':
        f = gzip.open(path, 'rt')
    else:
        raise NotImplementedError
        
    for line_number, line in enumerate(f):
        if line.startswith('>'):
            if name:
                ret[name] = ''.join(seq_list)
                seq_list = list()
            name = line[1:].strip()
        else:
            seq = ''.join(line.strip().split())
            seq_list.append(seq)
    if name:
        ret[name] = ''.join(seq_list)
    f.close()
    return ret

#1: original alignment file
#2: output masked alignment file
#3: masked proportion, range (0, 1]

inpath = sys.argv[1]
outpath = sys.argv[2]
thres = float(sys.argv[3])

# also allow for gz file reading
suffix = None
if inpath.split('.')[-1] in ['gzip', 'gz']:
    suffix = 'gz'

# record input alignment sequences
aln = read_fasta(inpath, suffix=suffix)
seqs = list(aln.values())
num_seq, seqlen = len(seqs), len(seqs[0])
gaps = [0] * seqlen

for seq in seqs:
    row = [c == '-' for c in seq]
    gaps = list(map(add, gaps, row))

# normalize by num_seq
for i in range(len(gaps)):
    gaps[i] /= num_seq

# retain columns that has gaps[i] < thres
retained = []
for i, x in enumerate(gaps):
    if x < thres:
        retained.append(i)

with open(outpath, 'w') as f:
    for name, seq in aln.items():
        f.write(f'>{name}\n')
        f.write(''.join([seq[i] for i in retained]) + '\n')
import os, sys
import numpy as np
from collections import defaultdict
from alignment_tools import Alignment

# indir
#indir = '../genome_known_51'
indir = sys.argv[1]
assert os.path.isdir(indir)

# genomes to use
#genome_file = '../genome_known_50.txt'
genome_file = sys.argv[2]
assert os.path.exists(genome_file)
with open(genome_file, 'r') as f:
    genomes = f.read().strip().split('\n')

# outdir
outdir = './genomes'
if not os.path.isdir(outdir):
    os.makedirs(outdir)

# obtain accession to taxid map
_map_file = sys.argv[3]
assert os.path.exists(_map_file)
accession_to_taxid = {}
with open(_map_file, 'r') as f:
    for line in f:
        taxid, acc = line.strip().split(',')
        accession_to_taxid[acc] = taxid

genomes_tsv = []
dnatypes_tsv = []
abundances_tsv = []
abundances = defaultdict(int)

for genome in genomes:
    inpath = os.path.join(indir, f'{genome}.zip')
    if not os.path.isfile(inpath):
        raise ValueError(f'{genome} file missing')

    outpath = os.path.join(outdir, f'{genome}.fna')
    if not os.path.exists(outpath):
        # deal with the file
        # this creates a directory "ncbi_dataset/data/[accession]"
        os.makedirs('tmp')
        os.system(f'cp {inpath} tmp/{genome}.zip; cd tmp; unzip {genome}.zip')
        aln_path = f'tmp/ncbi_dataset/data/{genome}/' + \
                os.popen(f'ls tmp/ncbi_dataset/data/{genome}').read().split('\n')[0]
        os.system(f'cp {aln_path} {outdir}/{genome}.fna')
        os.system('rm -r tmp')
    # read in the genome length to use for abundance calculation
    aln = Alignment(); aln.read_file_object(outpath)
    genome_size = sum([len(seq) for seq in aln.values()])
    taxid = accession_to_taxid[genome]
    genomes_tsv.append([taxid, outpath])
    dnatypes_tsv.append([taxid, 'complete sequence', 'circular'])
    abundances[taxid] += genome_size

# normalize abundances to be sum to 100(%)
#total_size = sum(abundances.values())
#for taxid in abundances.keys():
#    # round to 2 decimals
#    abundances[taxid] = round(abundances[taxid] / total_size * 100, 2)

# write genomes_tsv
with open('genomes.tsv', 'w') as f:
    for item in genomes_tsv:
        f.write('\t'.join(item) + '\n')

# write dnatypes_tsv
with open('dnatypes.tsv', 'w') as f:
    for item in dnatypes_tsv:
        f.write('\t'.join(item) + '\n')

# write abundances_tsv (as a suffix so we can control the number of reads
# to generate later)
with open('suffix_abundances.tsv', 'w') as f:
    total = sum(abundances.values())
    for k in abundances.keys():
        abundances[k] = abundances[k] * 100 / total
    #print(sum(abundances.values()), abundances)
    assert np.isclose(sum(abundances.values()), 100, rtol=1e-5)

    for k, v in abundances.items():
        f.write('{0}\t{1:.3f}\n'.format(k, v))

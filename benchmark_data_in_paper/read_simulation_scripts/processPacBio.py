import os, sys

# directory containing all generated read files (separated by entries in
# reference fasta file)
indir = sys.argv[1]

# output path containing all reads combined together
outpath = sys.argv[2]

if os.path.isdir(indir):
    raw_files = os.popen('ls {}'.format(indir)).read().strip().split('\n')
    
    # only use fastq files
    files = [f for f in raw_files if '.fastq' in f]
    refs = [f for f in raw_files if '.ref' in f]

    out = open(outpath, 'w', buffering=1)
    
    for i in range(len(files)):
        print(files[i])
        inpath = os.path.join(indir, files[i])

        # obtain serial number from refpath: should be the first line
        # right after ">"
        refpath = os.path.join(indir, refs[i])
        serial = None
        with open(refpath, 'r') as f:
            line = f.readline()
            serial = line[1:].split(' ')[0]

        with open(inpath, 'r') as f:
            # read in 4 lines a time
            lines = f.read().strip().split('\n')
            for j in range(0, len(lines), 4):
                idx = j // 4
                seq = lines[j+1] 
                if '.fq' in outpath or '.fastq' in outpath:
                    out.write('@{}-{}/1\n'.format(serial, idx))
                else:
                    out.write('>{}-{}/1\n'.format(serial, idx))

                out.write('{}\n'.format(seq))

                # write fastq information as well
                if '.fq' in outpath or '.fastq' in outpath:
                    out.write('+\n')
                    out.write('{}\n'.format(lines[j+3]))
    out.close()

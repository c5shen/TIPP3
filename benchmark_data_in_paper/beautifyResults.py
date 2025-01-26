import os, sys

queries = ['cami_long', 'cami_short',
        'known_exp2', 'mixed_exp2', 'novel_exp2',
        'known_pacbio_exp2', 'mixed_pacbio_exp2', 'novel_pacbio_exp2',
        'known_nanopore_exp2_2', 'mixed_nanopore_exp2_2', 'novel_nanopore_exp2_2',
        'known_metaphlan_exp2', 'mixed_metaphlan_exp2_2', 'novel_metaphlan_exp2_2',
        ]

method_map = {
        'bracken_filtered': 'bracken',
        'bracken_all_reads': 'bracken',
        'kraken_filtered': 'kraken_post',
        'kraken_all_reads': 'kraken_post',
        'metaphlan_all_reads': 'metaphlan_latest',
        'motus': 'motus',
        'metabuli_filtered': 'metabuli',
        'metabuli_all_reads': 'metabuli',
        'tipp3': 'pp-taxtastic-[90]',
        'tipp3-fast': 'bscampp-1000-[95]',
        }
tag_map = {
    'filtered': 'DATASETS_SMALL',
    'all_reads': 'all_reads'
    }
ranks = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']

rawdir = '/projects/illinois/eng/cs/warnow/chengze5/batch-SCAMPP/exp_tax_identification/placement_results/PROFILES'
done = 0
for querydir in queries:
    for method, oriname in method_map.items():
        # tag for filtered/all_reads
        intag = 'DATASETS_SMALL'
        outtag = None
        if len(method.split('_')) > 1:
            suffix = '_'.join(method.split('_')[1:])
            outtag = suffix
            intag = tag_map[suffix]
        elif method == 'motus':
            intag = 'all_reads'

        # alignment type --> only care for tipp3-fast, since it uses BLAST
        alndir = 'witch-alignment'
        if method == 'tipp3-fast':
            alndir = 'blast-alignment'

        # output directory name
        outdir = os.path.join('PROFILES', querydir, method)

        # input directory name --> only change is for TIPP3/TIPP3-fast
        indir = oriname
        if method == 'tipp3':
            indir = 'pp-taxtastic'
        elif method == 'tipp3-fast':
            indir = 'bscampp-1000'

        for rank in ranks:
            inpath = os.path.join(rawdir, querydir, alndir, indir,
                    '{}.{}.abundance.{}.tsv'.format(
                        oriname, intag, rank))
            #print(inpath, os.path.exists(inpath)); exit()
            if not os.path.exists(inpath):
                print('Skipping {} {}'.format(querydir, method))
                break
            if not os.path.isdir(outdir):
                os.makedirs(outdir)
                done += 1
            outpath = os.path.join(outdir,
                    '{}.abundance.{}.tsv'.format(
                    method, rank))
            os.system('cp {} {}'.format(inpath, outpath))
print(done)

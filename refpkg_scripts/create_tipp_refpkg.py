#!/usr/bin/env python
import os, sys, time
from subprocess import Popen, PIPE, STDOUT

import concurrent.futures
from concurrent.futures import ProcessPoolExecutor

def runner(args):
    cmd, marker, outdir, alignment_path, is_gzipped = args
    os.system(cmd)

    # est.aln.nuc.fasta derives from est.aln.nuc.fasta.gz, so no need to keep an extra copy
    if is_gzipped:
        os.system(f"rm {alignment_path}")

    # create the return statement
    to_write_map = {
            'alignment': 'est.aln.nuc.fasta',
            'size': 'est.aln.nuc.fasta.size',
            #'hmm': 'est.aln.nuc.fasta.hmm',
            'alignment-decomposition-tree': 'est.fasttree.tre',
            'placement-tree': 'est.fasttree.tre',
            'placement-tree-stats': 'est.fasttree.log',
            'seq-to-taxid-map': 'species.updated.mapping',
            'taxonomy': 'all_taxon.taxonomy',
            'additional-raxml-br-tree': 'est.raxml.bestTree.rooted',
            'additional-raxml-model-file': 'est.raxml.bestModel',
            }
    return marker, to_write_map

def runner2(args):
    to_write_map, gene, outdir, alignment_path = args

    # make sure which filetype we are reading (compressed or uncompressed)
    suffix = alignment_path.strip().split('.')[-1]

    for k, v in to_write_map.items():
        if k == 'size':
            # write size
            if suffix in ['gz', 'gzip']:
                cmd = ['gzip', '-d', alignment_path, '-c', '|', 'wc', '-l']
            else:
                cmd = ['cat', alignment_path, '|', 'wc', '-l']
            
            size = os.popen(' '.join(cmd)).read().strip().split()[0]
            size = int(size) // 2    # two lines per sequence
            size_path = os.path.join(outdir, 'est.aln.nuc.fasta.size')
            with open(size_path, 'w') as f:
                f.write('{}\n'.format(size))
        #elif k == 'hmm':
        #    # create an HMM (UPP way)
        #    hmm_path = os.path.join(dataset_pkgdir, 'est.aln.nuc.fasta.hmm')
        #    hmmcmd = 'hmmbuild --cpu 1 --dna --ere 0.59 --symfrac 0.0 --informat afa -o /dev/null {} {}'.format(
        #            hmm_path, alignment_path)
        #    os.system(hmmcmd)
    return gene

def clean_taxonomy_table(inpath, outpath):
    # kept columns
    columns = ['tax_id', 'parent_id', 'rank', 'tax_name',
        'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    lines = []
    with open(inpath, 'r') as f:
        col_to_idx = {}
        for line in f:
            parts = list(eval(line.strip()))
            # first line, map idx to column name
            if parts[0] == 'tax_id':
                for i, col in enumerate(parts):
                    if col in columns:
                        col_to_idx[col] = i
            else:
                row = [parts[col_to_idx[col]] for i, col in enumerate(columns)]
                lines.append(row)
    with open(outpath, 'w') as f:
        # write header
        f.write(','.join(columns) + '\n')
        for row in lines:
            # if empty, write NA
            f.write(','.join([x if x != '' else 'NA' for x in row]) + '\n')

def main():
    #1: working directory that contains all gene folders
    #2: ncbi_taxonomy.db path
    #3: output directory to write the refpkg to
    #4: number of threads to use for multiprocessing
    workdir = sys.argv[1]
    taxdb_path = sys.argv[2]
    basedir = sys.argv[3]
    num_threads = int(sys.argv[4])

    # default versioning and metadata. Change these as you like
    name = 'Chengze Shen'
    text = ("TIPP3 refpkg. Tree topology estimated using RAxML-ng with GTRGAMMA, "
            "branch lengths estimated with both RAxML-ng and FastTree-2.")
    version = '1.1'
    version = "markers-v4"
    
    # identify all genes from workdir
    # NOTE: the folder names should match gene names
    files = os.popen(f"ls {workdir}").read().strip().split('\n')
    genes = [f for f in files if os.path.isdir(f"{workdir}/{f}")]
    
    pkgdir = os.path.join(basedir, version)
    if not os.path.isdir(pkgdir):
        os.makedirs(pkgdir)

    # initialize multiprocessing
    pool = ProcessPoolExecutor(max_workers=num_threads)

    # create a "file-map-for-tipp.txt" file at the dir of the pkg
    filemap = open(os.path.join(pkgdir, 'file-map-for-tipp.txt'), 'w')

    # taxonomy dir
    taxdir = os.path.join(pkgdir, 'taxonomy')
    if not os.path.isdir(taxdir):
        os.makedirs(taxdir)

    # aggregate all unique taxids
    unique_taxids = set()
    for gene in genes:
        indir = os.path.join(workdir, gene)
        species_path = os.path.join(indir, 'species.updated.txt')
        with open(species_path, 'r') as f:
            for line in f:
                if not line.startswith('tax_id'):
                    unique_taxids.add(line.strip())
    all_taxid_path = os.path.join(taxdir, 'species.updated.txt')
    with open(all_taxid_path, 'w') as f:
        f.write('tax_id\n')
        f.write('\n'.join(list(unique_taxids)) + '\n')
    print(f"total unique taxids: {len(unique_taxids)}")
    
    # create taxtable
    taxtable_path = os.path.join(taxdir, 'taxonomy.table')
    cmd = f"taxit taxtable {taxdb_path} -i {all_taxid_path} -o {taxtable_path}"
    os.system(cmd)
    # clean up to create all_taxon.taxonomy
    taxonomy_path = os.path.join(taxdir, 'all_taxon.taxonomy')
    clean_taxonomy_table(taxtable_path, taxonomy_path)
    try:
        os.remove(taxtable_path)
    except IOError:
        pass
    filemap.write('taxonomy:taxonomy = taxonomy/all_taxon.taxonomy\n')

    # create a blast/ directory with all blast related items
    blastdir = os.path.join(pkgdir, 'blast')
    if not os.path.isdir(blastdir):
        os.makedirs(blastdir)
    # run makeblastdb command
    # FIRST, aggregate all sequences from each gene to form a single alignment.fasta file
    # THEN, create a seq2marker.tab that maps each sequence to their corresponding gene
    print("Aggregating all sequences to create a BLAST database...")
    blast_aln_path = os.path.join(blastdir, 'alignment.fasta')
    # remove old blast_aln_path to avoid appending extra lines
    if os.path.exists(blast_aln_path):
        os.system(f"rm {blast_aln_path}")
    seq2marker_path = os.path.join(blastdir, 'seq2marker.tab')
    for gene in genes:
        indir = os.path.join(workdir, gene)
        aln_path = os.path.join(indir, 'est.aln.nuc.fasta')
        gzipped_path = os.path.join(indir, 'est.aln.nuc.fasta.gz')
        # use cat to append alignment to the alignment.fasta file
        cmd = "awk '/^>/ {{print; next}} {{gsub(/-/, \"\"); print}}' >> {}".format(blast_aln_path)
        # use awk to append a tab and the gene name to each line from 
        cmd2 = "awk '/^>/ {{print substr($0, 2) \"\t\" \"{}\"}}' >> {}".format(gene, seq2marker_path)
        if not os.path.exists(aln_path):
            if not os.path.exists(gzipped_path):
                raise FileNotFoundError(f"Neither {aln_path} nor {gzipped_path} can be found.")
            else:
                # gzipped version
                cmd = f"gzip -d {gzipped_path} -c | " + cmd
                cmd2 = f"gzip -d {gzipped_path} -c | " + cmd2
        else:
            # regular file version
            cmd = f"cat {aln_path} | " + cmd
            cmd2 = f"cat {aln_path} | " + cmd2
        os.system(cmd)
        os.system(cmd2)
    print("Running makeblastdb...")
    cmd = f"makeblastdb -in {blast_aln_path} -out {blast_aln_path}.db -dbtype nucl"
    os.system(cmd)
    filemap.write('blast:database = blast/alignment.fasta.db\n')
    filemap.write('blast:seq-to-marker-map = blast/seq2marker.tab\n')

    futures = []
    for gene in genes:
        # using the taxonomy tree generated before alignment
        indir = os.path.join(workdir, gene)
        seq_info_path = os.path.join(indir, 'species.updated.mapping')
        aln_path = os.path.join(indir, 'est.aln.nuc.fasta')
        gzipped_path = os.path.join(indir, 'est.aln.nuc.fasta.gz')
        
        # if regular aln_path does not exist, gunzip to obtain it and remove the copy later
        is_gzipped = False
        if not os.path.exists(aln_path):
            if not os.path.exists(gzipped_path):
                raise FileNotFoundError(f"Neither {aln_path} nor {gzipped_path} can be found.")
            else:
                is_gzipped = True
                cmd = f"gzip -d {gzipped_path} -c > {aln_path}"
                os.system(cmd)

        ############ FastTree-2 branch length re-estimation
        tree_path = os.path.join(indir, 'est.fasttree.tre')
        tree_stat_path = os.path.join(indir, 'est.fasttree.log')

        outdir = os.path.join(pkgdir, f'{gene}.refpkg')
        cmd = ' '.join([
            f"taxit create -c -a '{name}' -d '{text}' -r {version}",
            f"--aln-fast {aln_path}",
            f"--tree-file {tree_path} --tree-stats {tree_stat_path}",
            f"-P {outdir} -l {gene} --no-reroot"])

        arg = (cmd, gene, outdir, aln_path, is_gzipped)
        futures.append(pool.submit(runner, arg))

    gene_to_map = {}
    for future in concurrent.futures.as_completed(futures):
        print(f"done creating refpkg: {gene}")
        gene, to_write_map = future.result()
        gene_to_map[gene] = to_write_map
        outdir = os.path.join(pkgdir, f'{gene}.refpkg')

        # copy additional files to outdir
        indir = os.path.join(workdir, gene)
        seq_info_path = os.path.join(indir, 'species.updated.mapping')
        os.system(f"cp {seq_info_path} {outdir}/")

        taxonomy_path = os.path.join(indir, 'all_taxon.taxonomy')
        os.system(f"cp {taxonomy_path} {outdir}/")

        raxml_model_path = os.path.join(indir, 'est.raxml.bestModel')
        raxml_tree_path = os.path.join(indir, 'est.raxml.bestTree.rooted')
        to_copy = [raxml_model_path, raxml_tree_path]
        for item in to_copy:
            os.system(f"cp {item} {outdir}/")

    sorted_genes = sorted(list(gene_to_map.keys()))
    for gene in sorted_genes:
        to_write_map = gene_to_map[gene]
        for k, v in to_write_map.items():
            filemap.write(f"{gene}:{k} = {gene}.refpkg/{v}\n")
    
    # MP-version for finding alignment size (and HMM creation, currently not used)
    futures = []
    for gene in sorted_genes:
        to_write_map = gene_to_map[gene]
        outdir = os.path.join(pkgdir, f'{gene}.refpkg')
        alignment_path = os.path.join(outdir, 'est.aln.nuc.fasta')
        args = (to_write_map, gene, outdir, alignment_path)
        futures.append(pool.submit(runner2, args))
    for future in concurrent.futures.as_completed(futures):
        gene = future.result()
        print(f"done computing num sequences: {gene}")

    filemap.close()
    pool.shutdown()

if __name__ == "__main__":
    main()

#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=cs
#SBATCH --mem=128GB

module load python/3.9.16
module load java/11
time=/usr/bin/time

# example to run TIPP3 with a given refpkg and a set of query reads
bin=../tipp3.py
inpath=data/nanopore.small.queries.fasta
refpkg=/home/chengze5/tallis/tipp3/tipp3-refpkg
t=16

scenario=1
if [[ $1 != "" ]]; then
    scenario=$1
fi

export TIPP_LOGGING_LEVEL=debug

if [[ $scenario == 1 ]]; then
    # BLAST alignment, BSCAMPP placement (TIPP3-fast)
    outdir=tipp3_scenario1
    $bin -i ${inpath} --reference-package ${refpkg} --outdir ${outdir} \
        --alignment-method blast --placement-method bscampp \
        -t $t --keeptemp
elif [[ $scenario == 2 ]]; then
    # BLAST alignment, pplacer-taxtastic placement
    outdir=tipp3_scenario2
    $bin -i ${inpath} --reference-package ${refpkg} --outdir ${outdir} \
        --alignment-method blast --placement-method pplacer-taxtastic \
        -t $t --keeptemp
elif [[ $scenario == 3 ]]; then
    # WITCH alignment, pplacer-taxtastic placement (TIPP3)
    outdir=tipp3_scenario3
    $bin -i ${inpath} --reference-package ${refpkg} --outdir ${outdir} \
        --alignment-method witch --placement-method pplacer-taxtastic \
        -t $t --keeptemp
elif [[ $scenario == 4 ]]; then
    # TIPP3-fast (BLAST+BSCAMPP) with .gz input type (fasta file)
    inpath=data/nanopore.queries.fasta.gz
    outdir=tipp3_scenario4
    $bin -i ${inpath} --reference-package ${refpkg} --outdir ${outdir} \
        --alignment-method blast --placement-method bscampp \
        -t $t --keeptemp
elif [[ $scenario == 5 ]]; then
    # TIPP3-fast (BLAST+BSCAMPP) with .gz input type (fasta file)
    inpath=data/illumina.small2.queries.fq.gz
    outdir=tipp3_scenario5
    $bin -i ${inpath} --reference-package ${refpkg} --outdir ${outdir} \
        --alignment-method blast --placement-method bscampp \
        -t $t --keeptemp
fi

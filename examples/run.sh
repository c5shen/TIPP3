#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=eng-instruction
##SBATCH --partition=cs
#SBATCH --mem=128GB

# example to run TIPP3 with a given refpkg and a set of query reads
bin=../run_tipp3.py
inpath=data/illumina.small.queries.fasta.gz
t=16

# NOTICE: supplement your own path of refpkg here
#refpkg=$HOME/Desktop/Research/phd_project/tipp3/tipp3-refpkg/
refpkg=$HOME/tallis/tipp3/tipp3-refpkg/
#refpkg=../refpkg_scripts/custom_tipp_refpkg

scenario=1
if [[ $1 != "" ]]; then
    scenario=$1
fi

export TIPP_LOGGING_LEVEL=debug

if [[ $scenario == 1 ]]; then
    # BLAST alignment, BSCAMPP placement (TIPP3-fast)
    outdir=tipp3_scenario1
    $bin abundance -i ${inpath} --reference-package ${refpkg} --outdir ${outdir} \
        --alignment-method blast --placement-method bscampp \
        -t $t
elif [[ $scenario == 2 ]]; then
    # BLAST alignment, pplacer-taxtastic placement
    outdir=tipp3_scenario2
    $bin abundance -i ${inpath} --reference-package ${refpkg} --outdir ${outdir} \
        --alignment-method blast --placement-method pplacer-taxtastic \
        -t $t
elif [[ $scenario == 3 ]]; then
    # WITCH alignment, pplacer-taxtastic placement (TIPP3)
    # also keep temporary files
    outdir=tipp3_scenario3
    $bin abundance -i ${inpath} --reference-package ${refpkg} --outdir ${outdir} \
        --alignment-method witch --placement-method pplacer-taxtastic \
        -t $t --keeptemp
elif [[ $scenario == 4 ]]; then
    # TIPP3-fast for species detection, with a custom detection threshold B=0.3
    outdir=tipp3_scenario4
    detection_threshold=0.3
    $bin detection -i ${inpath} --reference-package ${refpkg} --outdir ${outdir} \
        -t $t -B ${detection_threshold} --keeptemp
fi
#elif [[ $scenario == 4 ]]; then
#    # TIPP3-fast (BLAST+BSCAMPP) with .gz input type (fasta file)
#    inpath=data/nanopore.queries.fasta.gz
#    outdir=tipp3_scenario4
#    $bin abundance -i ${inpath} --reference-package ${refpkg} --outdir ${outdir} \
#        --alignment-method blast --placement-method bscampp \
#        -t $t
#elif [[ $scenario == 5 ]]; then
#    # TIPP3-fast (BLAST+BSCAMPP) with .gz input type (fastq file)
#    inpath=data/illumina.small.queries.fq.gz
#    outdir=tipp3_scenario5
#    $bin abundance -i ${inpath} --reference-package ${refpkg} --outdir ${outdir} \
#        --alignment-method blast --placement-method bscampp \
#        -t $t
#elif [[ $scenario == 6 ]]; then
#    # TIPP3-fast with species detection
#    outdir=tipp3_scenario6
#    $bin detection -i ${inpath} --reference-package ${refpkg} --outdir ${outdir} \
#        --alignment-method blast --placement-method bscampp \
#        -t $t
#fi

#!/bin/bash

# example to run TIPP3 with a given refpkg and a set of query reads
bin=../tipp3.py
#inpath=data/RplO_COG0200.queries.fasta
inpath=data/small_set.queries.fasta
#inpath=data/testing_known_2_dataset21.queries.fasta
refpkg=/home/chengze5/tallis/tipp3/tipp3-refpkg
outdir=./tipp3_output
t=16

scenario=1
if [[ $1 != "" ]]; then
    scenario=$1
fi

if [[ $scenario == 1 ]]; then
    $bin -i ${inpath} --reference-package ${refpkg} --outdir ${outdir} \
        --alignment-method blast --placement-method bscampp \
        -t $t
elif [[ $scenario == 2 ]]; then
    $bin -i ${inpath} --reference-package ${refpkg} --outdir ${outdir} \
        --alignment-method blast --placement-method pplacer-taxtastic \
        -t $t
fi

#!/bin/bash

# example to run TIPP3 with a given refpkg and a set of query reads
bin=../tipp3.py
inpath=data/RplO_COG0200.queries.fasta
refpkg=/home/chengze5/tallis/tipp3/tipp3-refpkg
outdir=./tipp3_output

$bin -i ${inpath} --reference-package ${refpkg} --outdir ${outdir} \
    --alignment-method blast

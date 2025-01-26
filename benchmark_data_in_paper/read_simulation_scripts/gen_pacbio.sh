#!/bin/bash

time=/usr/bin/time
seed=4221705

inpath=../genome_known_50.fasta
output=testing_known_pacbio_2
outdir=$(pwd)/$output
model=$HOME/tallis/softwares/pbsim-1.0.3-Linux-amd64/data/model_qc_clr

if [ ! -d $output ]; then
    mkdir -p $output
fi

{ $time -v pbsim ${inpath} \
        --prefix ${outdir}/${output} \
        --depth 20 \
        --length-min 400 \
        --length-mean 3000 \
        --accuracy-mean 0.78 \
        --accuracy-sd 0.07 \
        --seed $seed \
        --model_qc $model ; } 2> runtime_gen_pacbio.txt

# aggregate all reads
python3 processPacBio.py ${output} ${output}.fasta
python3 processPacBio.py ${output} ${output}.fq

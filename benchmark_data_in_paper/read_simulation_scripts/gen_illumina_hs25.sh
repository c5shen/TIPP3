#!/bin/bash

time=/usr/bin/time

# SAMPLE INPUT PATH TO A FASTA FILE CONTAINING ALL GENOMES TO SIMULATE
# READS FROM
inpath=genome_known_50.fasta
output=testing_known_2_dataset2

# NOTE - this is using the command for the novel_100 and novel_33,
# for some reason the parameters are different for known_51, but I think
# it shouldn't be (and the original parameters are missing mlen and std)
{ $time -v art_illumina -ss HS25 -sam \
    -i ${inpath} -p -l 150 -f 20 -m 200 -s 10 \
    -o $output ; } 2> runtime_gen_illumina_hs25.txt

fq2fa ${output}1.fq ${output}1.fasta

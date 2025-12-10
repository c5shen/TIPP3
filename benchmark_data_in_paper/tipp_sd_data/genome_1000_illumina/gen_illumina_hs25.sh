#!/bin/bash
#SBATCH --time=03:59:59
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=secondary
#SBATCH --mem=64GB
#SBATCH -o %j.out
#SBATCH -e %j.out
#SBATCH --mail-user=chengze5@illinois.edu
#SBATCH --mail-type=END

time=/usr/bin/time

queryset=testing_detection_1000_cov_1
coverage=1
inpath=../genome_known_1000.fasta
output=${queryset}_dataset2

# NOTE - this is using the command for the novel_100 and novel_33,
# for some reason the parameters are different for known_51, but I think
# it shouldn't be (and the original parameters are missing mlen and std)
{ $time -v art_illumina -ss HS25 -sam \
    -i ${inpath} -p -l 150 -f $coverage -m 200 -s 10 \
    -o $output ; } 2> runtime_gen_illumina_hs25.txt

fq2fa ${output}1.fq ${output}1.fasta

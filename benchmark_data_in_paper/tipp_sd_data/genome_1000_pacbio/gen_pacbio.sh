#!/bin/bash
#SBATCH --time=03:59:59
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=secondary
#SBATCH --mem=64GB
#SBATCH -o ../logs/%j.out
#SBATCH -e ../logs/%j.out
#SBATCH --mail-user=chengze5@illinois.edu
#SBATCH --mail-type=END

time=/usr/bin/time
seed=86810585

#inpath=../genome_known_1000.fasta
name=testing_detection_1000_pacbio_cov_1
model=$HOME/tallis/softwares/pbsim-1.0.3-Linux-amd64/data/model_qc_clr
coverage=1

idx=0
inpath=../genome_known_1000_${idx}.fasta
while [ -s $inpath ]; do
    output=${name}_${idx}
    outdir=$(pwd)/$output
    if [ ! -d $output ]; then
        mkdir -p $output
    fi

    { $time -v pbsim ${inpath} \
            --prefix ${outdir}/${output} \
            --depth $coverage \
            --length-min 400 \
            --length-mean 3000 \
            --accuracy-mean 0.78 \
            --accuracy-sd 0.07 \
            --seed $seed \
            --model_qc $model ; } 2> runtime_gen_pacbio_${idx}.txt

    # aggregate all reads
    python3 ../processPacBio.py ${output} ${output}.fasta
    python3 ../processPacBio.py ${output} ${output}.fq
    # remove the outdir
    rm -r $outdir
    # increment idx and update inpath
    idx=$(($idx + 1))
    inpath=../genome_known_1000_${idx}.fasta
done

# combine outputs to a single fasta file
idx=0
outpath=${name}_${idx}.fasta
finalpath=${name}.fasta
touch $finalpath
while [ -s $outpath ]; do
    cat $outpath >> $finalpath
    idx=$(($idx + 1))
    outpath=${name}_${idx}.fasta
done

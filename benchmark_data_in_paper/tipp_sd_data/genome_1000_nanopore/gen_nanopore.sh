#!/bin/bash
#SBATCH --account=warnow-cs-eng
#SBATCH --time=03:59:59
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --partition=secondary
#SBATCH --mem=128GB
#SBATCH -o ../logs/%j.out
#SBATCH -e ../logs/%j.out
#SBATCH --mail-user=chengze5@illinois.edu
#SBATCH --mail-type=END

time=/usr/bin/time
source ~/.bash_profile

#module load anaconda/2023-Jul/3
module load anaconda3
source activate nanosim-env
#export PYTHONPATH=~/.conda/envs/nanosim-env/lib/python3.8/site-packages:$PYTHONPATH
pythonbin=~/.conda/envs/nanosim-env/bin/python

seed=67109123
t=16

output=testing_detection_1000_nanopore_cov_1
outdir=$(pwd)/$output
model=../nanosim-metagenome-pretrained-model/training
bin=$HOME/tallis/softwares/NanoSim-v3.2.2/src/simulator.py

# number of reads to simulate 
num_reads=1000000

if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi

# run createNanoporeGenomes.py to create genomes.tsv, dnatypes.tsv
genomedir=../genome_known_1000
genomelist=../genome_known_1000.txt
accmap=../taxid_to_accession.csv
python3 ../createNanoporeGenomes.py ${genomedir} ${genomelist} ${accmap}

genome_path=genomes.tsv
abundance_path=abundances.tsv
dnatype_path=dnatypes.tsv

# create abundance_path
echo -e "Size\t${num_reads}" > ${abundance_path}
cat suffix_abundances.tsv >> ${abundance_path}

# output FASTQ file and convert later
{ $time -v $pythonbin $bin metagenome -gl ${genome_path} -a ${abundance_path} \
    -dl ${dnatype_path} -c $model -o $output/nanopore --seed ${seed} \
    --fastq -t $t ; } 2> ./runtime_gen_nanopore.txt

# output files --> (used) nanopore_sample0_aligned_reads.fastq
#               -> (not used) nanopore_sample0_unaligned_reads.fastq
cp ${output}/nanopore_sample0_aligned_reads.fastq ./${output}.fq
fq2fa ./${output}.fq ./${output}.fasta

# get motus long reads
#./motus_prep_long.sh ${output}.fq

# zip output dir and remove
zip -r ${output}.zip ${output}/; rm -r ${output}

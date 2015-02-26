#! /bin/sh
#SBATCH -p general
#SBATCH --mem=10G
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -J rapsearch_%a
#SBATCH -o out/rapsearch_%a.out
#SBATCH -e err/rapsearch_%a.err
#SBATCH --time=48:00:00
#SBATCH --array=1-384

module load centos6/biopython-1.62b_python-2.7.3
#module load centos6/numpy-1.7.1_python-2.7.3
module load bio/bowtie2-2.1.0
module load bio/samtools-0.1.18

if [ ! -f tmp/annot/allannot_parsed_$SLURM_ARRAY_TASK_ID.fasta.m8 ];
then
python create_annotations.py tmp/fasta_files/$SLURM_ARRAY_TASK_ID.fasta ../uniref90 ../uniref50 /n/huttenhower_lab/tools/rapsearch2 10
fi

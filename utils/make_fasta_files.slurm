#! /bin/sh
#SBATCH -p serial_requeue
#SBATCH --mem=2G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J make_fastas_buccal_mucosa
#SBATCH -o out/make_fastas_buccal_mucosa.out
#SBATCH -e err/make_fastas_buccal_mucosa.err
#SBATCH --time=2:00:00

module load centos6/biopython-1.62b_python-2.7.3
#module load centos6/numpy-1.7.1_python-2.7.3
module load bio/bowtie2-2.1.0
module load bio/samtools-0.1.18

python map_runner.py -m input/mapper_buccal_mucosa.txt -w 2

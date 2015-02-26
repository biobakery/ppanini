#! /bin/sh
#SBATCH -p general
#SBATCH --mem=1G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J make_sams_stool
#SBATCH -o out/make_sams_stool.out
#SBATCH -e err/make_sams_stool.err
#SBATCH --time=144:00:00

module load centos6/biopython-1.62b_python-2.7.3
#module load centos6/numpy-1.7.1_python-2.7.3
module load bio/bowtie2-2.1.0
module load bio/samtools-0.1.18

python map_runner.py -m input/mapper_stool.txt -w 1

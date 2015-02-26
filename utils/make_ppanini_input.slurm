#! /bin/sh
#SBATCH -p general
#SBATCH --mem=20G
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -J anares_preppanini
#SBATCH -o out/anares_preppanini.out
#SBATCH -e err/anares_preppanini.err
#SBATCH --time=24:00:00

module load centos6/biopython-1.62b_python-2.7.3
#module load centos6/numpy-1.7.1_python-2.7.3
module load bio/bowtie2-2.1.0
module load bio/samtools-0.1.18

python preprocess_genedata_table.py -m input/mapper_anterior_nares.txt -p 10 -umap ../map_UniRef90_50_jan0115.txt -u90 ../uniref90 -u50 ../uniref50 -w 3 -o anterior_nares_ppanini_input.txt -r /n/huttenhower_lab/tools/rapsearch2


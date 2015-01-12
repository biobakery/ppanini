import os
import sys
import pdb
import re
import argparse
import numpy
import subprocess

def parse_mapper(mapper_file):

	mapper = {}
	mapper_foo = open(mapper_file)
	mapper_foo = mapper_foo.readlines()

	header = [re.sub('[\t\n\r]', '', i) for i in mapper_foo[0].split('\t')]

	for i, val_i in enumerate(mapper_foo):
		split_val_i = [re.sub('[\t\n\r]', '', i) for i in val_i.split('\t')]
		mapper[i] = {}
		for j, val_j in enumerate(split_val_i):
			mapper[i][header[j]] = val_j

	return mapper

def get_abundance(mapper, all_paths):

	for sample in mapper:
		assembly_x = mapper[sample]['#ASSEMBLIES']
		reads_x = mapper[sample]['#READS']
		os.system('bowtie2-build --quiet '+all_paths['assemblies']+'/'+mapper[sample]['#ASSEMBLIES']+' tmp/'+mapper[sample]['#ASSEMBLIES']+'_index')
		os.system('tar -xOvf '+all_paths['reads']+'/'+mapper[sample]['#READS']+' | \
				   bowtie2 -x tmp/'+mapper[sample]['#ASSEMBLIES']+'_index -U - --no-unal --very_sensitive | \
				   samtools view -bS - > tmp/'+mapper[sample]['#ASSEMBLIES']+'.bam')
		os.system('samtools sort tmp/'+mapper[sample]['#ASSEMBLIES']+'.bam '+'tmp/'+mapper[sample]['#ASSEMBLIES']+'.sorted')
		os.system('samtools index tmp/'+mapper[sample]['#ASSEMBLIES']+'.sorted.bam')
		os.system('samtools idxstats tmp/'+mapper[sample]['#ASSEMBLIES']+'.sorted.bam > '+'tmp/'+mapper[sample]['#ASSEMBLIES']+'.txt')
		mapper[sample]['abund_file'] = 'tmp/'+mapper[sample]['#ASSEMBLIES']+'.txt'

	return read_abundance_tables(mapper, all_paths)
		
def read_abundance_tables(mapper, all_paths):

	pass
		#tar -xOvf reads__file | bowtie2 -x filename_index -U - --no-unal --very_sensitive | \
		#samtools view -bS - > newfile.bam
		#samtools sort filename newfilename.sorted
		#samtools index newfilename.sorted.bam
		#samtools idxstats newfilename > statsfilename


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-a', '--assemblies_folder', help='Folder containing ORF-called contig assemblies of samples')
	parser.add_argument('-r', '--reads_folder', help='Folder containing reads from samples')
	parser.add_argument('-n', '--niche_present', choices=[True, False], help='Is Niche-classification present? [True or False]]')
	parser.add_argument('-m', '--mapper_file', help='Mapper file associating read files with their corresponding assembly files and, if available, GFF3 files and NICHE informations')
	parser.add_argument('-g', '--gff3_folder', help='Folder containing GFF3 files')
	parser.add_argument('-u90x', '--uniref90_xml', help='UniRef90 XML file')
	parser.add_argument('-u90', '--uniref90_fasta', help='UniRef90 fasta file')
	parser.add_argument('-u50', '--uniref50_fasta', help='UniRef50 fasta file')
	parser.add_argument('-rs2', '--rapsearch2_path', help='Path to run RAPSEARCH2')
	parser.add_argument('-st', '--samtools_path', help='Path to run SAMTOOLS')
	parser.add_argument('-bwt2', '--bowtie2_path', help='Path to run bowtie2_path')

	args = parser.parse_args()

	all_paths = {'assemblies': args.assemblies_folder, \
				 'reads': args.reads_folder, \
				 'uniref90x': args.uniref90_xml, \
				 'uniref90': args.uniref90_fasta, \
				 'uniref50': args.uniref50_fasta, \
				 'rapsearch2': args.rapsearch2_path, \
				 'samtools': args.samtools_path, \
				 'bowtie2': args.bowtie2_path}

	if args.gff3_folder is not None:
		all_paths['gff3s'] = args.gff3_folder

	mapper = parse_mapper(args.mapper_file)
	get_abundance(mapper, all_paths)
	get_annotations(mapper, all_paths)


	[metadata, uniref_gis, gis_unannotated, gene_ids, data_matrix] = read_gene_table(args.input_table)
	all_centroids = get_centroids(gis_unannotated, args.fasta_folder, metadata, args.usearch_folder, uniref_gis)
	centroids_data_matrix = get_centroids_table(gene_ids, all_centroids, data_matrix, metadata)
	[centroid_prev_abund, all_prevalence, all_mean_abund, flag] = get_prevalence_abundance(centroids_data_matrix, metadata)

	if flag:
		#Niche-classification present
		imp_centroids = get_important_niche_centroids(centroid_prev_abund, all_prevalence, all_mean_abund, args.output_folder)
	else:
		imp_centroids = get_important_centroids(centroid_prev_abund, all_prevalence, all_mean_abund, args.output_folder)

	
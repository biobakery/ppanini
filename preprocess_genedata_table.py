import os
import sys
import pdb
import re
import argparse
import numpy
import subprocess
import multiprocessing

from src import annotations_mod
from src import abundance_mod


def parse_mapper(mapper_file):
	'''Reads the mapper_file that contains all the paths and metadata for input

	Input: mapper_file = path_to_mapper_file
	Format_1: #SAMPLE\t?NICHE\tFASTA\tBAMS\n
	Format_2: #SAMPLE\t?NICHE\tFASTA\tSAMS\n
	Format_3: #SAMPLE\tREADS\tASSEMBLIES\tGFF3S\t?NICHE\n
	
	?NICHE: optional metadata

	Output: mapper_final = {sample:{'READS': path_to_reads_file, #Format_3
		   					 'CONTIG_ASSEMBLIES': path_to_assemblies, #Format_3
		   					 'FASTAS': path_to_fasta_file, #Format_1 or Format_2
		   					 'SAMS': path_to_sam_file, #Format_2
		   					 'BAMS': path_to_bam_file, #Format_1
		   					 'GFF3S': path_to_gff3_file, #Format_3
		   					 'NICHE': niche_metadata}, ...}'''
# #SAMPLE<\t>?READS<\t>?ASSEMBLIES<\t>?NICHE<\t>?FASTA<\t>?GFF3S<\t>?SAMS<\t>?BAMS
	mapper = {}
	mapper_foo = open(mapper_file)
	mapper_foo = mapper_foo.readlines()

	header = [re.sub('[#\t\n\r]', '', i) for i in mapper_foo[0].split('\t')]

	niche_flag = [True for i in header if i == 'NICHE']
	
	for i, val_i in enumerate(mapper_foo[1:]):
		split_val_i = [re.sub('[\t\n\r]', '', i) for i in val_i.split('\t')]
		mapper[i] = {}
		for j, val_j in enumerate(split_val_i):
			mapper[i][header[j]] = val_j

	mapper_final = {}

	for i in mapper:
		mapper_final[mapper[i]['SAMPLES']] = mapper[i]

	return [mapper_final, niche_flag]

def read_id_mapping(mapping_file):
	umap90_50 = {}
	
	with open(mapping_file) as foo:
		for line in foo.readlines():
			split_i = [re.sub('[\t\n\r]', '', i) for i in line.split('\t')]
			umap90_50[split_i[0]] = split_i[1]
	return umap90_50

def generate_gene_table(abundance_dict, annotations_dict, all_paths, niche_flag, mapper, output_table):
	'''Writes the final gene abundance table with annotations for PPANINI input
	
	Input: abundance_dict = {sample: {gene: abundance, ...},...}	   
		   annotations_dict = {sample: {gene: annotation, ...},...}
		   all_paths = {'uniref90': path_to_uniref90_index, 
		   				'uniref50': path_to_uniref50_index, 
		   				'umap90_50': path_to_uniref90_uniref50_mapping}
		   niche_flag = Boolean [True, False] If niches present in mapper_file
		   mapper = {sample:{'READS': path_to_reads_file, 
		   					 'CONTIG_ASSEMBLIES': path_to_assemblies, 
		   					 'FASTAS': path_to_fasta_file,
		   					 'SAMS': path_to_sam_file,
		   					 'BAMS': path_to_bam_file,
		   					 'GFF3S': path_to_gff3_file,
		   					 'NICHE': niche_metadata}, ...}
		   output_table = path_to_output_table'''

	samples = abundance_dict.keys()
	fasta_row = [mapper[i]['FASTAS'] for i in samples]

	if niche_flag:
		niche_row = [mapper[i]['NICHE'] for i in samples]

	with open(output_table, 'w') as foo:
		print 'Writing files'
		if niche_flag:
			foo.writelines([str.join('\t', ['#NICHE']+niche_row)+'\n']) #header

		foo.writelines([str.join('\t', ['#FASTAS']+fasta_row)+'\n']) #header
		foo.writelines([str.join('\t', ['#GENES']+samples)+'\n']) #header
		
		for i, sample in enumerate(samples):
			for gene in abundance_dict[sample]:
		
				abund_x_i = abundance_dict[sample][gene]
				data_row = numpy.zeros(len(samples))
				data_row[i] = abund_x_i
				str_data_row = [str(ele) for ele in data_row]
				
				if gene in annotations_dict[sample]: 
					annot_x_i = annotations_dict[sample][gene]
					if annot_x.startswith('UniRef90'):
						umap_i = re.sub('[\t\n\r]','',subprocess.check_output('grep -w '+annot_x_i+' '+all_paths['uniref_map']).split('\t')[-1])
						annot_x = annot_x_i + '|' + umap90_50[annot_x_i]
					else:
						annot_x = 'UniRef90_unknown|' + annot_x_i 
				else:
					annot_x = 'UniRef90_unknown|UniRef50_unknown'
				
				foo.writelines([str.join('\t', [gene+'|'+annot_x]+str_data_row)+'\n'])


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-m', '--mapper_file', help='Mapper file associating read files with their corresponding assembly files and, if available, GFF3 files and NICHE informations')
	parser.add_argument('-p', '--processes', help='Folder containing GFF3 files', default=4)
	parser.add_argument('-umap', '--uniref90_50', help='UniRef90 XML file')
	parser.add_argument('-u90', '--uniref90_fasta', help='UniRef90 INDEX file')
	parser.add_argument('-u50', '--uniref50_fasta', help='UniRef50 INDEX file')
	parser.add_argument('-w', '--workflow', choices=['1', '2', '3'], help='Workflow type Choices:[1, 2, 3]; \
																   1: BAM and FASTA files; \
																   2: SAM and FASTA FILES; \
																   3: CONTIG ASSEMBLIES, READS and GFF3 files')
	parser.add_argument('-o', '--output_table', help='Gene Table to write', default=sys.stdout)

	args = parser.parse_args()

	nprocesses = int(args.processes)

	try:
		os.mkdir('tmp')
	except:
		pass

	all_paths = {'uniref_map': args.uniref90_50, \
				 'uniref90': args.uniref90_fasta, \
				 'uniref50': args.uniref50_fasta, }
	
	print 'Step1: Reading MAPPER FILE, located at: ' + args.mapper_file

	[mapper, niche_flag] = parse_mapper(args.mapper_file)

	print 'Step2: Calculating abundance for genes via bowtie2/samtools'
	abundance_dict = abundance_mod.get_abundance(mapper, nprocesses, int(args.workflow)) #sample:contig: mapped reads/seq_length
	
	print 'Step3: Annotating genes via RAPSEARCH2'
	if int(args.workflow) == 3:
		[annotations_dict, gene_contig_mapper] = annotations_mod.get_annotations_fromcontigs(mapper, all_paths, nprocesses)
		##Update abundance_dict from {sample:{contig:abundance}} to {sample:{gene:abundance}} via gene_contig_mapper
		abundance_dict = abundance_mod.update_abundance_dict(abundance_dict, gene_contig_mapper)
	else:
		annotations_dict = annotations_mod.get_annotations_fromgenes(mapper, all_paths, nprocesses)

	print 'Calling generate_gene_table'
	
	generate_gene_table(abundance_dict, annotations_dict, all_paths, niche_flag, mapper, args.output_table)
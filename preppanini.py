
import os
import sys
import pdb
import re
import argparse
import numpy
import subprocess
import multiprocessing

from src import create_idxstats
from src import create_fastas
from src import create_annotations
from src import write_ppanini_table
from utils import utilities

basename = ''

def parse_mapper(mapper_file):
	'''Reads the mapper_file that contains all the paths and metadata for input

	Input: mapper_file = path_to_mapper_file
	Format_1: #SAMPLE\t?NICHE\tFASTAS\tBAMS\n
	Format_2: #SAMPLE\t?NICHE\tFASTAS\tSAMS\n
	Format_3: #SAMPLE\tREADS\tASSEMBLIES\tGFF3S\t?NICHE\n

	?NICHE: optional metadata

	Output: mapper_final = {sample:{'READS': path_to_reads_file, #Format_3
	                                             'CONTIG_ASSEMBLIES': path_to_assemblies, #Format_3
	                                             'FASTAS': path_to_fasta_file, #Format_1 or Format_2
	                                             'SAMS': path_to_sam_file, #Format_2
	                                             'BAMS': path_to_bam_file, #Format_1
	                                             'GFF3S': path_to_gff3_file, #Format_3
	                                             'NICHE': niche_metadata}, ...}'''
	# #SAMPLE<\t>?READS<\t>?ASSEMBLIES<\t>?NICHE<\t>?FASTAS<\t>?GFF3S<\t>?SAMS<\t>?BAMS
	mapper = {}
	mapper_foo = open(mapper_file)
	mapper_foo = mapper_foo.readlines()

	header = [re.sub('[#\t\n\r]', '', i).upper() for i in mapper_foo[0].split('\t')]

	niche_flag = [True for i in header if i == 'NICHE']
	gff3_flag = [True for i in header if i == 'GFF3S']

	for i, val_i in enumerate(mapper_foo[1:]):
	    split_val_i = [re.sub('[\t\n\r]', '', i) for i in val_i.split('\t')]
	    mapper[split_val_i[0]] = {}
	    for j, val_j in enumerate(split_val_i):
	            mapper[split_val_i[0]][header[j]] = val_j
	
	return [mapper, niche_flag, gff3_flag]

def abundance_module(mapper, workflow):
	
	keywords = {2: 'SAMS', 1: 'BAMS'}
	abundance_method = {3: create_idxstats.generate_abundance_viabwt2, \
						2: create_idxstats.generate_abundance_viasam, \
						1: create_idxstats.generate_abundance_viabam}

	if workflow == 3:
		for sample in mapper:
			mapper[sample]['abundance_file'] = abundance_method[workflow](mapper[sample]['CONTIG_ASSEMBLIES'], mapper[sample]['READS'], sample)[0]
	else:
		for sample in mapper:
			mapper[sample]['abundance_file'] = abundance_method[workflow](mapper[sample][keywords[workflow]], sample)[0]
	return mapper
	
def annotation_module(mapper, all_paths, nprocesses):
	## create pangenome insert
	annotations_dict = create_annotations.generate_annotation(mapper[sample]['FASTAS'], all_paths, nprocesses, basename)

	return annotations_dict

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-m', '--mapper_file', help='Mapper file associating read files with their corresponding assembly files and, if available, GFF3 files and NICHE informations', required=True)
	parser.add_argument('-p', '--processes', help='Number of threads', default=4)
	parser.add_argument('--uniref90_50', help='IDMAPPING FILE for UniRef90 ID to UniRef50 ID mapping')
	parser.add_argument('--uniref90', help='UniRef90 INDEX file')
	parser.add_argument('--uniref50', help='UniRef50 INDEX file')
	parser.add_argument('--diamond', help='Path to folder containing DIAMOND') #add to be in path??
	parser.add_argument('--usearch', help='Path to folder containing USEARCH') #add to be in path?
	parser.add_argument('-w', '--workflow', choices=['1', '2', '3'], help='Workflow type Choices:[1, 2, 3]; \
	                                                                        1: BAM and FASTA files; \
	                                                                        2: SAM and FASTA FILES; \
	                                                                        3: CONTIG ASSEMBLIES, READS and GFF3 files', required=True)
	parser.add_argument('--annotation_only', default=False, help='Perform annotation only')
	parser.add_argument('--abundance_only', default=False, help='Perform abundance only')
	parser.add_argument('--write_tables_only', default=False, help='Write table only (Annotations and Abundance files exist under /n/annot/ and /tmp/idxstats/')
	parser.add_argument('--basename', help='BASENAME for all the output files')
	parser.add_argument('-o', '--output_table', help='Gene Table to write', default=sys.stdout)

	args = parser.parse_args()
	annotation_only = bool(args.annotation_only)
	abundance_only = bool(args.abundance_only)
	write_tables_only = bool(args.write_tables_only)

	nprocesses = int(args.processes)
	workflow = int(args.workflow)
	mapper_file = args.mapper_file
	basename = args.basename

	if not basename:
		basename = mapper_file.split('.')[0].split('/')[-1]

	gene_contig_mapper = {}
	n = 1

	#################################
	###Initiating folder hierarchy
	#################################
	utilities.create_folders(['tmp', 'tmp/fasta_files'])
	#################################

	all_paths = {'uniref_map': args.uniref90_50, \
	             'uniref90': args.uniref90, \
	             'uniref50': args.uniref50, \
	             'diamond': args.diamond, \
	             'usearch': args.usearch}

	[mapper, niche_flag, gff3_flag] = parse_mapper(mapper_file)

	###################################
	##MODULE1: abundance_only or ALL
	###################################
	if not annotation_only and not write_tables_only: 
		print 'Step'+str(n)+': Creating abundance tables via SAMTOOLS'
		mapper = abundance_module(mapper, workflow)
		n +=1

	###################################
	##MODULE2: annotation_only or ALL
	###################################
	if not write_tables_only and not abundance_only: 
		print 'Step'+str(n)+': Mapping genes against UniRef90 and UniRef50 via DIAMOND'
		if workflow == 3:
			print 'Step'+str(n)+': Pulling genes from contigs'
			n +=1
			for sample in mapper:
				mapper[sample]['FASTAS'] = 'tmp/fasta_files/'+sample+'.fasta'
				gene_contig_mapper[sample] = create_fastas.pullgenes_fromcontigs(mapper[sample]['CONTIG_ASSEMBLIES'], \
														  mapper[sample]['GFF3S'], \
														  mapper[sample]['FASTAS'])
		annotations_dict = annotation_module(mapper, all_paths, nprocesses)
		n +=1
	
	###################################
	##MODULE3: write_tables_only or ALL
	###################################
	if not annotation_only and not abundance_only: 

		if write_tables_only:
			#create annotations_dict
			annotations_dict = {}
			for sample in mapper:
				annotations_dict[sample] = create_annotations.read_dict('tmp/annot/'+sample+'.m8')
			#create abundance_files keys in mapper
			for sample in mapper:
				mapper[sample]['abundance_file'] = 'tmp/idxstats/'+sample+'.txt'
			#mapper for gene_contig for workflow3?
			if workflow == 3:
				gene_contig_mapper = {}
				for sample in mapper:
					mapper[sample]['FASTAS'] = 'tmp/fasta_files/'+sample+'.fasta'
					[gene_contig_mapper_i, gene_start_stop_i, contig_gene_mapper_i] = create_fastas.read_gff3(mapper[sample]['GFF3S'])
					gene_contig_mapper[sample] = gene_contig_mapper_i

		print 'Step'+str(n)+': Compiling Abundance dictionary'
		abundance_dict = create_idxstats.read_abundance_tables(mapper)
		n +=1
		
		if workflow == 3:
			abundance_dict = create_idxstats.update_abundance_dict(abundance_dict, gene_contig_mapper)
		
		print 'Step'+str(n)+': Writing to PPANINI-input format table...'
		n +=1
		write_ppanini_table.generate_gene_table(abundance_dict, annotations_dict, all_paths, niche_flag, mapper, args.output_table)

import os
import sys
import pdb
import re
import argparse
import numpy
import subprocess
import multiprocessing
import logging

from src import quantify_genes
from src import utilities
from src import annotate_genes
from src import write_ppanini_table
from utils import utilities

basename = ''
logger = logging.getLogger(__name__)

def parse_mapper(mapper_file):
	'''Reads the mapper_file that contains all the paths and metadata for input

	Input: mapper_file = path_to_mapper_file
	Format_1: #SAMPLE\t?NICHE\tFAAS\tBAMS\n
	Format_2: #SAMPLE\t?NICHE\tFAAS\tSAMS\n
	Format_3: #SAMPLE\tREADS\tASSEMBLIES\tGFF3S\t?NICHE\n

	?NICHE: optional metadata

	Output: mapper_final = {sample:{'READS': path_to_reads_file, #Format_3
	                                             'CONTIG_ASSEMBLIES': path_to_assemblies, #Format_3
	                                             'FAAS/FNAS': path_to_fasta_file, #Format_1 or Format_2
	                                             'SAMS': path_to_sam_file, #Format_2
	                                             'BAMS': path_to_bam_file, #Format_1
	                                             'GFF3S': path_to_gff3_file, #Format_3
	                                             'NICHE': niche_metadata}, ...}'''
	# #SAMPLE<\t>?READS<\t>?ASSEMBLIES<\t>?NICHE<\t>?FASTAS<\t>?GFF3S<\t>?SAMS<\t>?BAMS
	mapper = {}
	mapper_foo = open(mapper_file)

	header = []
	for line in mapper_foo:
		if line.startswith('#'):
			header = [re.sub('[#\t\n\r]', '', i).upper() for i in line.split('\t')]
		else:
			split_val_i = [re.sub('[\t\n\r]', '', i) for i in line.split('\t')]
			mapper[split_val_i[0]] = {}
			for i, val in enumerate(split_val_i):
				mapper[split_val_i[0]][header[i]] = val
	flags = {}
	for i in header:
		flags[i] = True
	labels = ['SAMS', 'BAMS', 'NICHE', 'GFF3S', 'ABUNDANCE_TABLES', 'ANNOTATION', 'READS', 'CONTIG_ASSEMBLIES', 'FAAS', 'FNAS']
	for i in labels:
		if not i in flags:
			flags[i] = False
	
	return [mapper, flags]

def abundance_module(mapper, workflow):
	
	keywords = {2: 'SAMS', 1: 'BAMS'}
	abundance_method = {3: quantify_genes.generate_abundance_viabwt2, \
						2: quantify_genes.generate_abundance_viasam, \
						1: quantify_genes.generate_abundance_viabam}

	if workflow == 3:
		for sample in mapper:
			mapper[sample]['abundance_file'] = abundance_method[workflow](mapper[sample]['CONTIG_ASSEMBLIES'], mapper[sample]['READS'], sample)[0]
	else:
		for sample in mapper:
			mapper[sample]['abundance_file'] = abundance_method[workflow](mapper[sample][keywords[workflow]], sample)[0]
	return mapper
	
def annotation_module(mapper, all_paths, nprocesses):
	## create pangenome insert
	annotations_dict = annotate_genes.generate_annotation(mapper[sample]['FAAS'], all_paths, nprocesses, basename)

	return annotations_dict

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-m', '--mapper_file', help='Mapper file containing paths to data', required=True)
	parser.add_argument('--basename', help='BASENAME for all the output files')
	parser.add_argument('--bypass_abundance', default=False, action='store_true', help='Bypass quantifying abundance')
	parser.add_argument('--bypass_annotation', default=False, action='store_true', help='Bypass annotating genes')
	parser.add_argument('--bypass_clust', default=False, action='store_true', help='Bypass annotating genes')
	parser.add_argument('--bypass_write_table', default=False, action='store_true', help='Bypass writing table')
	parser.add_argument('--usearch', default=False, help='Path to USEARCH') #add to be in path?
	parser.add_argument('--vsearch', default=False, help='Path to VSEARCH') #add to be in path?
	parser.add_argument('--diamond', default=False, help='Path to DIAMOND') #add to be in path??
	parser.add_argument('--rapsearch', default=False, help='Path to RAPSEARCH') #add to be in path??
	parser.add_argument('--threads', help='Number of threads', default=1)
	parser.add_argument('--uniref90', help='UniRef90 INDEX file')
	parser.add_argument('--uniref50', help='UniRef50 INDEX file')
	parser.add_argument('--to_normalize', default=False, action='store_true', help='Default HUMAnN2 table; if sam-idxstats table; enable')


	parser.add_argument('--uniref90_50', help='IDMAPPING FILE for UniRef90 ID to UniRef50 ID mapping')
	parser.add_argument('-w', '--workflow', choices=['1', '2', '3'], help='Workflow type Choices:[1, 2, 3]; \
	                                                                        1: BAM and FASTA files; \
	                                                                        2: SAM and FASTA FILES; \
	                                                                        3: CONTIG ASSEMBLIES, READS and GFF3 files', required=True)
	parser.add_argument('--annotation_only', default=False, help='Perform annotation only')
	parser.add_argument('--abundance_only', default=False, help='Perform abundance only')
	parser.add_argument('--write_tables_only', default=False, help='Write table only (Annotations and Abundance files exist under /n/annot/ and /tmp/idxstats/')
	
	# parser.add_argument('-o', '--output_table', help='Gene Table to write', default=sys.stdout)
	
	parser.add_argument('--log_level',default='DEBUG', help='Choices: [DEBUG, INFO, WARNING, ERROR, CRITICAL]')

	args = parser.parse_args()
	nprocesses = int(args.threads)
	workflow = int(args.workflow)
	mapper_file = args.mapper_file
	basename = args.basename

	#flags
	norm_flag = args.norm_flag
	annotation_only = bool(args.annotation_only)
	abundance_only = bool(args.abundance_only)
	write_tables_only = bool(args.write_tables_only)
	
	if not basename:
		basename = mapper_file.split('.')[0].split('/')[-1]
	
	output_folder = basename
	output_table = output_folder+'/'+basename+'_ppanini.txt'

	utilities.create_folders([output_folder, output_folder+'/tmp'])

	log_file = output_folder+'/'+basename+'.log'
	logging.basicConfig(filename=log_file, \
						format='%(asctime)s - %(name)s - %(levelname)s: %(message)s', \
						level=getattr(logging, args.log_level), \
						filemode='w', \
						datefmt='%m/%d/%Y %I:%M:%S %p')
	
	[mapper, flags] = parse_mapper(mapper_file)

	paths_dict = {}

	if flags['GFF3S']:
		paths_dict['FASTAS'] = output_folder+'/tmp/fasta_files'
		utilities.create_folders([paths_dict['FASTA_FILES']])
		for sample in mapper:
			mapper[sample]['FNAS'] = paths_dict['FASTA_FILES']+'/'+sample+'.fna'
			mapper[sample]['FAAS'] = paths_dict['FASTA_FILES']+'/'+sample+'.faa'

		for sample in mapper:
				utilities.pullgenes_fromcontigs(mapper[sample]['CONTIG_ASSEMBLIES'], \
											    mapper[sample]['GFF3S'], \
											    mapper[sample]['FNAS'], \
												mapper[sample]['FAAS'])
	if not args.bypass_abundance or not flags['ABUNDANCE_TABLES']:

		paths_dict['ABUNDANCE_TMP'] = output_folder+'/tmp/abundance_tmp'
		paths_dict['ABUNDANCE_TABLES'] = output_folder+'/tmp/abundance_tables'
		utilities.create_folders([paths_dict['ABUNDANCE_TMP'], paths_dict['ABUNDANCE_TABLES']])
	
		if flags['BAMS']:
			for sample in mapper:
				mapper[sample]['ABUNDANCE_TABLES'] = quantify_genes.generate_abundance_viabam(mapper[sample]['BAMS'], \
																							  sample, \
																							  paths_dict)
		elif flags['SAMS']:
			for sample in mapper:
				mapper[sample]['ABUNDANCE_TABLES'] = quantify_genes.generate_abundance_viasam(mapper[sample]['SAMS'], \
																							  sample, \
																							  paths_dict)
		elif flags['CONTIG_ASSEMBLIES']:
			for sample in mapper:
				try:
					mapper[sample]['ABUNDANCE_TABLES'] = quantify_genes.generate_abundance_viabam(mapper[sample]['FNAS'], \
																							  mapper[sample]['READS'], \
																							  sample, \
																							  paths_dict)
				except KeyError:
					raise Exception('No FNAS found. Either provide in mapper or add GFF3s to pull from contigs')
		else:
			raise Exception('No Abundance information; See --bypass_abundance if no calculation required')		
	
	if not args.bypass_annotation or not flags['ANNOTATION']:
		paths_dict['ANNOTATION_TMP'] = output_folder+'/tmp/annotation_tmp'

		whole_genome_catalog = paths_dict['ANNOTATION_TMP']+'/'+basename+'.fasta'
		utilities.create_folders([paths_dict['ANNOTATION_TMP']])
		annotation_table = paths_dict['ANNOTATION_TMP']+'/'+basename+'.annotations.txt'
		# mapper['ANNOTATION'] = annotation_table
		#To ensure all samples have FAAs
		for sample in mapper:
			if not 'FAAS' in mapper[sample]:
				if not 'FNAS' in mapper[sample]:
					raise Exception('No FNAS or FAAS found. Either provide in mapper or add GFF3s to pull from contigs')
				else:
					paths_dict['FASTAS'] = output_folder+'/tmp/fasta_files'
					utilities.create_folders([paths_dict['FASTA_FILES']])
					mapper[sample]['FAAS'] = paths_dict['FASTA_FILES']+'/'+sample+'.faa'
					faa = utilities.read_fasta(mapper[sample]['FNAS'])
					utilities.write_fasta(faa, mapper[sample]['FAAS'])
		os.system('cat '+' '.join([mapper[sample]['FAAS'] for sample in mapper])+' > '+whole_genome_catalog)

		genome_catalog = '' #The catalog run against UniRef90
		
		if not args.bypass_clust:
			clust_method = 'vsearch'
			gene_centroids_file_path = paths_dict['ANNOTATION_TMP']+'/'+basename+'.centroids.fasta'
			gene_centroid_clusters_file_path = paths_dict['ANNOTATION_TMP']+'/'+basename+'.uc'
			if args.usearch:
				clust_method = args.usearch
				utilities.run_uclust(clust_method, \
									 whole_genome_catalog, \
									 gene_centroids_file_path, \
									 gene_centroid_clusters_file_path, \
									 0.9, \
									 nprocesses)
			else:
				if args.vsearch:
					clust_method = args.vsearch #assumes vsearch in path if not provided
				utilities.run_vclust(clust_method, \
								 	 whole_genome_catalog, \
								 	 gene_centroids_file_path, \
								 	 gene_centroid_clusters_file_path, \
								 	 0.9, \
								 	 nprocesses)
				genome_catalog = gene_centroids_file_path
		else:
			genome_catalog = whole_genome_catalog

		##Run UniRef now
		out_u90_fname = paths_dict['ANNOTATION_TMP']+'/'+basename+'.centroids.u90'
		out_u50_fname = paths_dict['ANNOTATION_TMP']+'/'+basename+'.centroids.u50'
		u50_gene_input = paths_dict['ANNOTATION_TMP']+'/'+basename+'.centroids.u50input.fasta'

		search_method = 'diamond'
		if args.rapsearch:
			search_method = args.rapsearch
		else:
			if args.diamond:
				search_method = args.diamond
		
		run_diamond(genome_catalog, args.uniref90, out_u90_fname, nprocesses)
		centroid_sequences = utilities.read_fasta(genome_catalog)
		[centroid_annotations90, diamond50_seqs] = parse_annotation_table(u90_out_fname+'.m8', centroid_sequences, 90.0)
		centroid_annotations = centroid_annotations90

		if not diamond50_seqs == {}:
			utilities.write_fasta(diamond50_seqs, u50_gene_input, True)
			run_diamond(u50_gene_input, args.uniref50, u50_out_fname, nprocesses)
			[centroid_annotations50, diamondukn_seqs] = parse_annotation_table(u50_out_fname+'.m8', diamond50_seqs, 50.0)
		
		for gid in centroid_annotations50:
			if gid not in centroid_annotations:
				centroid_annotations[gid] = sample_annotations50[gid]
			else:
				raise Exception('GeneID for UniRef50 exists in UniRef90;\nSearching twice for the same gene;\nError in parse_annotation_table\n')
		
		if not args.bypass_clust:
			centroid_gis = get_clusters_dict(gene_centroid_clusters_file_path)
			annotation_dict = annotate_genes.get_annotations_dict(centroid_annotations, centroid_gis)
		else:
			annotation_dict = centroid_annotations

		utilities.write_dict(annotation_dict, annotation_table)

	if not args.bypass_write_table:
		if args.to_normalize:
			abundance_dict = quantify_genes.read_abundance_tables(mapper, True)
		else:
			abundance_dict = quantify_genes.read_abundance_tables(mapper, False)

		try:
			annotation_dict = annotation_dict
		except: #if Annotations have not been performed
			paths_dict['ANNOTATION_TMP'] = output_folder+'/tmp/annotation_tmp'
			utilities.create_folders([paths_dict['ANNOTATION_TMP']])
			annotation_table = paths_dict['ANNOTATION_TMP']+'/'+basename+'.annotations.txt'
			paths_annots = []
			try:
				for sample in mapper:
					if not mapper[sample]['ANNOTATION'] in paths_annots:
						paths_annots += [mapper[sample]['ANNOTATION']]
			except KeyError:
				raise Exception('Annotations have not been performed\nAdd ANNOTATION to mapper file\n')
			os.system('cat '+' '.join(paths_annots)+' > '+annotation_table)
			annotation_dict = utilities.read_dict(annotation_table)

		generate_gene_table(abundance_dict, annotations_dict, niche_flag, mapper, output_table):

	# gene_contig_mapper = {}

	# all_paths = {'uniref_map': args.uniref90_50, \
	#              'uniref90': args.uniref90, \
	#              'uniref50': args.uniref50, \
	#              'diamond': args.diamond, \
	#              'usearch': args.usearch}

	# ###################################
	# ##MODULE1: abundance_only or ALL
	# ###################################
	# if not annotation_only and not write_tables_only: 
	# 	print 'Step'+str(n)+': Creating abundance tables via SAMTOOLS'
	# 	mapper = abundance_module(mapper, workflow)
	# 	n +=1

	# ###################################
	# ##MODULE2: annotation_only or ALL
	# ###################################
	# if not write_tables_only and not abundance_only: 
	# 	print 'Step'+str(n)+': Mapping genes against UniRef90 and UniRef50 via DIAMOND'
	# 	n +=1
		
	# 	if workflow == 3:
	# 		print 'Step'+str(n)+': Pulling genes from contigs'
	# 		n +=1
	# 		for sample in mapper:
	# 			gene_contig_mapper[sample] = utilities.pullgenes_fromcontigs(mapper[sample]['CONTIG_ASSEMBLIES'], \
	# 													  mapper[sample]['GFF3S'], \
	# 													  mapper[sample]['FAAS'])
	# 	annotations_dict = annotation_module(mapper, all_paths, nprocesses)
	
	# ###################################
	# ##MODULE3: write_tables_only or ALL
	# ###################################
	# if not annotation_only and not abundance_only: 
	# 	if write_tables_only:
	# 		#create annotations_dict
	# 		annotations_dict = {}
	# 		for sample in mapper:
	# 			annotations_dict[sample] = utilties.read_dict(mapper[sample]['ANNOTS'])
	# 		#create abundance_files keys in mapper
	# 		for sample in mapper:
	# 			mapper[sample]['abundance_file'] = mapper[sample]['ABUNDS']
	# 		#mapper for gene_contig for workflow3?
	# 		if workflow == 3:
	# 			gene_contig_mapper = {}
	# 			for sample in mapper:
	# 				# if not 'FASTAS' in mapper[sample]:
	# 				# 	mapper[sample]['FASTAS'] = 'tmp/fasta_files/'+sample+'.fasta'
	# 				[gene_contig_mapper_i, gene_start_stop_i, contig_gene_mapper_i] = utilities.read_gff3(mapper[sample]['GFF3S'])
	# 				gene_contig_mapper[sample] = gene_contig_mapper_i
	# 	print 'Step'+str(n)+': Compiling Abundance dictionary'
	# 	abundance_dict = quantify_genes.read_abundance_tables(mapper, norm_flag)
	# 	n +=1
	# 	if workflow == 3:
	# 		abundance_dict = quantify_genes.update_abundance_dict(abundance_dict, gene_contig_mapper)
		
	# 	print 'Step'+str(n)+': Writing to PPANINI-input format table...'
	# 	n +=1
	# 	write_ppanini_table.generate_gene_table(abundance_dict, annotations_dict, all_paths, niche_flag, mapper, output_table)
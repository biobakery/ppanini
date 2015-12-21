import pdb
import os
import sys
import re
import argparse
import numpy
import subprocess
import multiprocessing
import logging

from .. import quantify_genes
from .. import utilities
from .. import annotate_genes
from .. import write_ppanini_table

basename = ''
logger = logging.getLogger(__name__)

def check_user_options(args):
	if not args.bypass_annotation: 
		if not args.uniref90:
			sys.exit("Please bypass annotation or provide uniref90 database path!")
	elif args.uniref90:
		sys.exit("bypass annotation or provide uniref90 database path!")
	return
		
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
	logger.debug('Parsing Mapper file: '+mapper_file)

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

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-m', '--mapper-file', dest='mapper_file', help='Mapper file containing paths to data', required=True)
	parser.add_argument('--basename', help='BASENAME for all the output files')
	parser.add_argument('--bypass-abundance', dest= 'bypass_abundance', default=False, action='store_true', help='Bypass quantifying abundance')
	parser.add_argument('--bypass-annotation', dest= 'bypass_annotation', default=False, action='store_true', help='Bypass annotating genes')
	parser.add_argument('--bypass-clust', dest='bypass_clust', default=False, action='store_true', help='Bypass annotating genes')
	parser.add_argument('--bypass-write-table', dest= 'bypass_write_table', default=False, action='store_true', help='Bypass writing table')
	parser.add_argument('--usearch', default=False, help='Path to USEARCH') #add to be in path?
	parser.add_argument('--vsearch', default=False, help='Path to VSEARCH') #add to be in path?
	parser.add_argument('--diamond', default=False, help='Path to DIAMOND') #add to be in path??
	parser.add_argument('--rapsearch', default=False, help='Path to RAPSEARCH') #add to be in path??
	parser.add_argument('--threads', help='Number of threads', default=1)
	parser.add_argument('--uniref90', help='UniRef90 INDEX file')
	parser.add_argument('--to-normalize', dest='to_normalize', default=False, action='store_true', help='Default HUMAnN2 table; if sam-idxstats table; enable')
	parser.add_argument('--log-level', dest='log_level', default='DEBUG', help='Choices: [DEBUG, INFO, WARNING, ERROR, CRITICAL]')

	args = parser.parse_args()
	nprocesses = int(args.threads)
	mapper_file = args.mapper_file
	basename = args.basename
	check_user_options(args)
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
		logger.debug('Pulling genes from contigs')
		paths_dict['FASTA_FILES'] = output_folder+'/tmp/fasta_files'
		utilities.create_folders([paths_dict['FASTA_FILES']])
		for sample in mapper:
			mapper[sample]['FNAS'] = paths_dict['FASTA_FILES']+'/'+sample+'.fna'
			mapper[sample]['FAAS'] = paths_dict['FASTA_FILES']+'/'+sample+'.faa'
		for sample in mapper:
				utilities.pullgenes_fromcontigs(mapper[sample]['CONTIG_ASSEMBLIES'], \
												mapper[sample]['GFF3S'], \
												mapper[sample]['FNAS'], \
												mapper[sample]['FAAS'])

	if not args.bypass_abundance and not flags['ABUNDANCE_TABLES']:
		logger.debug('Quantifying gene abundance in samples')
		paths_dict['ABUNDANCE_TMP'] = output_folder+'/tmp/abundance_tmp'
		paths_dict['ABUNDANCE_TABLES'] = output_folder+'/tmp/abundance_tables'
		utilities.create_folders([paths_dict['ABUNDANCE_TMP'], paths_dict['ABUNDANCE_TABLES']])
	
		if flags['BAMS']:
			logger.debug('Running SAMTOOLS on BAM FILES')
			for sample in mapper:
				mapper[sample]['ABUNDANCE_TABLES'] = quantify_genes.generate_abundance_viabam(mapper[sample]['BAMS'], \
																							  sample, \
																							  paths_dict)
		elif flags['SAMS']:
			logger.debug('Running SAMTOOLS on SAM FILES')
			for sample in mapper:
				mapper[sample]['ABUNDANCE_TABLES'] = quantify_genes.generate_abundance_viasam(mapper[sample]['SAMS'], \
																							  sample, \
																							  paths_dict)
		elif flags['CONTIG_ASSEMBLIES']:
			logger.debug('Running BOWTIE2 on READS against CONTIG_ASSEMBLIES')
			for sample in mapper:
				try:
					mapper[sample]['ABUNDANCE_TABLES'] = quantify_genes.generate_abundance_viabwt2(mapper[sample]['FNAS'], \
																							  mapper[sample]['READS'], \
																							  sample, \
																							  paths_dict)
				except KeyError:
					logger.debug('Gene fasta files not found; Please provide them as FNAs or GFF3 files to pull genes from contigs')
					raise Exception('No FNAS found. Either provide in mapper or add GFF3s to pull from contigs')
		else:
			raise Exception('No Abundance information; See --bypass_abundance if no calculation required')		
	
	if not args.bypass_annotation and not flags['ANNOTATION']:
		logger.debug('Running ANNOTATION MODULE')
		paths_dict['ANNOTATION_TMP'] = output_folder+'/tmp/annotation_tmp'

		whole_genome_catalog = paths_dict['ANNOTATION_TMP']+'/'+basename+'.fasta'
		utilities.create_folders([paths_dict['ANNOTATION_TMP']])
		annotation_table = paths_dict['ANNOTATION_TMP']+'/'+basename+'.annotations.txt'

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
			logger.debug('Running CLUSTERING before ANNOTATION')
			#clust_method = 'vsearch'
			gene_centroids_file_path = paths_dict['ANNOTATION_TMP']+'/'+basename+'.centroids.fasta'
			gene_centroid_clusters_file_path = paths_dict['ANNOTATION_TMP']+'/'+basename+'.uc'
			if args.usearch:
				clust_method = args.usearch
				annotate_genes.run_uclust(clust_method, \
									 whole_genome_catalog, \
									 gene_centroids_file_path, \
									 gene_centroid_clusters_file_path, \
									 0.9, \
									 nprocesses)
			elif args.vsearch:
					clust_method = args.vsearch #assumes vsearch in path if not provided
					annotate_genes.run_vclust(clust_method, \
											whole_genome_catalog, \
											gene_centroids_file_path, \
											gene_centroid_clusters_file_path, \
											0.9)
			else:
				raise Exception('No clustering software found: Please use --vsearch or --usearch to specify the path to software')
			genome_catalog = gene_centroids_file_path
		else:
			logger.debug('BYPASSING CLUSTERING')
			genome_catalog = whole_genome_catalog

		##Run UniRef now
		out_u90_fname = paths_dict['ANNOTATION_TMP']+'/'+basename+'.centroids.u90'

		#paths_dict['diamond'] = 'diamond'
		if args.rapsearch:
			paths_dict['rapsearch'] = args.rapsearch
			annotate_genes.run_rapsearch(genome_catalog, args.uniref90, out_u90_fname, nprocesses, paths_dict)
		elif args.diamond:
			paths_dict['diamond'] = args.diamond
			annotate_genes.run_diamond(genome_catalog, args.uniref90, out_u90_fname, nprocesses, paths_dict)
		else:
			raise Exception('No similarity search software found: Please use --diamond or --rapsearch to specify the path to software')
		logger.debug('Running SEARCH against UniRef')
		
		
		centroid_sequences = utilities.read_fasta(genome_catalog)
		[centroid_annotations90, diamond50_seqs] = annotate_genes.parse_annotation_table(out_u90_fname+'.m8', centroid_sequences, 90.0)
		centroid_annotations = centroid_annotations90

		if not args.bypass_clust:
			centroid_gis = annotate_genes.get_clusters_dict(gene_centroid_clusters_file_path)
			annotation_dict = annotate_genes.get_annotations_dict(centroid_annotations, centroid_gis)
		else:
			annotation_dict = centroid_annotations
		logger.debug('Annotation DICT CREATED')
		utilities.write_dict(annotation_dict, annotation_table)

	if not args.bypass_write_table:
		
		logger.debug('Reading abundance tables')
		if args.to_normalize:
			[abundance_dict, samples] = quantify_genes.read_abundance_tables(mapper, True)
		else:
			[abundance_dict, samples] = quantify_genes.read_abundance_tables(mapper, False)
		logger.debug('Abundance tables parsed')

		try:
			annotation_dict = annotation_dict
		except: #if Annotations have not been performed
			paths_dict['ANNOTATION_TMP'] = output_folder + '/tmp/annotation_tmp'
			utilities.create_folders([paths_dict['ANNOTATION_TMP']])
			annotation_table = paths_dict['ANNOTATION_TMP'] + '/' + basename + '.annotations.txt'
			paths_annots = []
			try:
				for sample in mapper:
					if not mapper[sample]['ANNOTATION'] in paths_annots:
						paths_annots += [mapper[sample]['ANNOTATION']]
			except KeyError:
				raise Exception('Annotations have not been performed\nAdd ANNOTATION to mapper file\n')
			os.system('cat '+' '.join(paths_annots)+' > ' + annotation_table)
			annotation_dict = utilities.read_dict(annotation_table)
			logger.debug('Annotations parsed')
		logger.debug('Writing PPANINI table')
		write_ppanini_table.generate_gene_table(abundance_dict, annotation_dict, flags['NICHE'], mapper, output_table, samples)

if __name__ == '__main__':
	main()

import os
import sys
import pdb
import re
import argparse
import numpy
import logging

from src import utilities
from src import annotate_genes

logger = logging.getLogger(__name__)

basename = ''
temp_folder = ''

numpy.seterr(divide='ignore', invalid='ignore')

def read_gene_table(gene_table_fname):
	'''Returns the different elements from the gene table

	Input: gene_table_fname = Filename of the gene_table

	Output: metadata = [metadata strings]; Rows with # as first character in table
			uniref_gis = {UniRef_XYZ: [list of gene ids]}
			gis_unannotated = {sample_id: [gene ids]}
			gene_ids = [List of all gene ids]
			data_matrix = The abundance table [[0, 0,.],[0, 0,.],...]'''
	
	logger.debug('read_gene_table')

	gene_table = open(gene_table_fname) 
	metadata = []
	gis_unannotated, uniref_gis, gis_dm = {}, {}, {} # Cluster IDs: [Gene IDs]
	samples = []
	uniref_dm = {}
	gi_dm = {}
	for line in gene_table:
		if line.startswith('#'):
			metadata += [line]
		else:
			split_i = line.split('\t')
			annot = split_i[0].split('|') #geneID column split 
			u90_annot = [i for i in annot if 'UniRef90' in i][0]
			
			data_row = numpy.array([float(i) for i in split_i[1:]])
			
			if 'UniRef90_unknown' == u90_annot:
				try: #same name
					gis_dm[annot[0]] += data_row
				except:
					gis_dm[annot[0]] = data_row
			else: #same uniref90 id
				try:
					uniref_dm[u90_annot] += data_row
				except KeyError:
					uniref_dm[u90_annot] = data_row			
	return [uniref_dm, gis_dm, metadata]

def get_centroids_fromUCLUST(gene_centroid_clusters_file_path, genes):

	logger.debug('get_centroids_fromUCLUST')
	
	#genes list of unannotated genes

	cluster_dict = {}
	genes_clustered = []
	foo = open(gene_centroid_clusters_file_path)
	for line in foo:
		if line.startswith('H'):
			split_line = [re.sub('[\r\t\n]','', i) for i in line.split('\t')[-2:]]
			if sum([1 for i in split_line if i in genes]):
				try:
					cluster_dict[split_line[1]] += [split_line[0]]
					genes_clustered += [split_line[0]]
				except:
					cluster_dict[split_line[1]] = [split_line[0], split_line[1]]
					genes_clustered += [split_line[0], split_line[1]]
	for i in genes:
		if not i in genes_clustered:
			cluster_dict[i] = [i]
	return cluster_dict

def get_centroids(uniref_dm, gi_dm, usearch_folder, uclust_file, gene_catalog, nprocesses):
	'''Returns the dict of all centroids containing clusters of gene IDs

	Input:	gis_unannotated = {sample_id: [gene ids]}
			metadata = [metadata strings]; Rows with # as first character in table
			uniref_gis = {UniRef_XYZ: [list of gene ids]}
			gis_unannotated = {sample_id: [gene ids]}
			usearch_folder = Location of the USEARCH program or in path if not provided
			uclust_file = path to USEARCH_UCLUST FILE from PREPPANINI

	Output: all_centroids = {gene_centroid : [List of gene ids]}'''
	
	logger.debug('get_centroids')

	centroids_fasta = {}

	if not uclust_file:
		centroid_gis = get_clusters(gene_catalog, usearch_folder, nprocesses) #all UniRef90_unknowns are clustered across samples
	else:
		centroid_gis = get_centroids_fromUCLUST(uclust_file, gi_dm.keys())

	gc_dm = {}
	for centroid in centroid_gis:
		for gene in centroid_gis[centroid]:
			try:
				gc_dm[centroid] += gi_dm[gene]
			except:
				gc_dm[centroid] = gi_dm[gene]
	
	for centroid in uniref_gis:
		gc_dm[centroid] = uniref_gis[centroid]
	
	return gc_dm


def get_clusters(centroids_fasta, args, nprocesses): #ONLY FOR THE UNIREF UNANNOTATED
	'''Returns the dict of unannotated gene centroids containing clusters of genes at 90% similarity

	Input:	centroids_fasta = {UniRef90_unknown: [List of gene ids]}
			usearch_folder = Location of the USEARCH program

	Output: centroid_gis = {gene_centroid: [List of genes in the cluster]}'''

	logger.debug('get_clusters')

	allgenes_file_path = centroids_fasta #temp_folder+'/'+basename+'_centroids_for_clustering.fasta'
	gene_centroids_file_path = temp_folder+'/'+basename+'_centroids.fasta'
	gene_centroid_clusters_file_path = temp_folder+'/'+basename+'_clusters.uc'

	# utilities.write_fasta(centroids_fasta, allgenes_file_path, True) #ensures all clustering happens to FAAs
	
	clust_method = 'vsearch'
	if args.usearch:
		clust_method = args.usearch
		annotate_genes.run_uclust(clust_method, allgenes_file_path, gene_centroids_file_path, gene_centroid_clusters_file_path, 0.9, nprocesses)
	else:
		if args.vsearch:
			clust_method = args.usearch
		annotate_genes.run_vclust(clust_method, allgenes_file_path, gene_centroids_file_path, gene_centroid_clusters_file_path, 0.9, nprocesses)
		
	
	centroid_gis = annotate_genes.get_clusters_dict(gene_centroid_clusters_file_path)

	return centroid_gis


def get_centroids_table(all_centroids, metadata):
	'''Returns data matrix containing gene centroids and abundance per sample

	Input:	metadata = [metadata strings]; Rows with # as first character in table
			gene_ids = [List of all gene ids]
			data_matrix = The abundance table [[0, 0,.],[0, 0,.],...]
			all_centroids = {gene_centroid : [List of gene ids]}

	Output: centroids_data_matrix = {gene_centroid: [Gene centroid abundance across samples]}'''
	
	logger.debug('get_centroids_table')
	# pdb.set_trace()
	centroids_data_matrix = []
	centroids_list = []
	
	for centroid in all_centroids:
		centroids_list +=[centroid]
		centroids_data_matrix +=[all_centroids[centroid]]

	#NORMALIZATION PER SAMPLE
	centroids_data_matrix = numpy.array(centroids_data_matrix)
	norm_data_matrix = centroids_data_matrix/sum(centroids_data_matrix)
	norm_data_matrix = norm_data_matrix*1e6

	gene_centroids_table_file_path = temp_folder+'/'+basename+'_gene_centroids_table.txt'
	
	with open(gene_centroids_table_file_path,'w') as foo:
		foo.writelines(metadata)
		# foo.writelines([str.join('\t', ['#SAMPLES'] + metadata[-1].split('\t')[1:])])
		for i, val in enumerate(centroids_list):
			foo.writelines([str.join('\t', [val] + [str(j) for j in norm_data_matrix[i]]) + '\n'])
			

	return [norm_data_matrix, centroids_list]


def get_prevalence_abundance(centroids_data_matrix, centroids_list, metadata):
	'''Returns the dict of centroids with their prevalence and abundance

	Input:	centroids_data_matrix = {gene_centroid: [Gene centroid abundance across samples]}
			metadata = [metadata strings]; Rows with # as first character in table

	Output: centroid_prev_abund = {centroid: {'abund': mean abundance, 'prev': prevalence}}
			all_prevalence = [List of all observed gene centroid prevalence values (>0) across samples]
			all_abund = [List of all calculated gene centroid abundance across samples]
			flag = True (if NICHE PRESENT) or False(if NICHE ABSENT)'''
	
	logger.debug('get_prevalence_abundance')

	centroid_prev_abund_file_path = temp_folder+'/'+basename+'_centroid_prev_abund.txt'
	
	[niche_line, ind] = utilities.is_present(metadata, '#NICHE')
	
	if niche_line:
		niche_flag = True
		[centroid_prev_abund, all_alpha_prev, all_mean_abund] = get_niche_prevalence_abundance (centroids_data_matrix, centroids_list, niche_line)
		return [centroid_prev_abund, all_alpha_prev, all_mean_abund, niche_flag]
	else:
		niche_flag = False
		centroid_prev_abund = {}
		all_prevalence = [] 
		all_abund = []
		
		for centroid in centroids_data_matrix:
			#abund only where the gene is present in sample
			abund_i = [i for i in  centroids_data_matrix[centroid] if i > 0]
			centroid_prev_abund[centroid] = {'abund': numpy.mean(numpy.array(abund_i)), \
											 'prev': float(sum(numpy.array(centroids_data_matrix[centroid]) > 0)/float(len(centroids_data_matrix[centroid])))}
			
			all_prevalence += [centroid_prev_abund[centroid]['prev']]
			all_abund += abund_i

		write_prev_abund_matrix(centroid_prev_abund, centroid_prev_abund_file_path)
		
		return [centroid_prev_abund, all_prevalence, all_abund, niche_flag]

def get_niche_prevalence_abundance(centroids_data_matrix, centroids_list, niche_line):
	'''Returns the dict of centroids with their prevalence and abundance

	Input:	centroids_data_matrix = {gene_centroid: [Gene centroid abundance across samples]}
			metadata = [metadata strings]; Rows with # as first character in table
			line = line from mapper file containing NICHE information

	Output: centroid_prev_abund = {centroid: {'abund': float_mean abundance, 
											  'a_prev': {'NICHE**': prevalence of centroid across samples within niche, ...}},
											  'b_prev': float_median of alpha prevalence within niches for centroid}}
			all_alpha_prev = {niche: List of all observed gene centroid alpha prevalence values (>0) across samples within each niche}
			all_mean_abund = [List of all calculated mean gene centroid abundance across samples]'''

	logger.debug('get_niche_prevalence_abundance')
	
	centroid_prev_abund_file_path = temp_folder+'/'+basename+'_centroid_prev_abund.txt'
	
	niches = {}
	
	split_i = [re.sub('[\r\t\n]', '', i) for i in niche_line.split('\t')[1:]]
	for i, val in enumerate(split_i):
		try:
			niches[val] += [i]
		except KeyError:
			niches[val] = [i]
	
	niches_label = niches.keys()
	centroid_prev_abund = {}

	all_alpha_prev = {}
	all_alpha_abund = []

	for niche in niches_label:
		all_alpha_prev[niche] = []

	for i, centroid in enumerate(centroids_list):
		centroid_prev_abund[centroid] = {}
		a_prev = {}
		a_abund = {}
		for niche in niches_label:
			a_prev[niche] = 0
			a_abund[niche] = []
 			for ind in niches[niche]:
 				check_i = int(centroids_data_matrix[i][ind] > 0) #present in a sample
				a_prev[niche] += check_i
				if check_i:
					a_abund[niche] += [centroids_data_matrix[i][ind]]

			a_prev[niche] = float(a_prev[niche])/float(len(niches[niche]))
			all_alpha_prev[niche] += [a_prev[niche]]

			# all_alpha_abund[niche] += a_abund[niche]
			if not a_abund[niche]: #if gene completely empty throughout samples
				a_abund[niche] = [0]
			a_abund[niche] = numpy.mean(numpy.array(a_abund[niche]))
		
		max_mean_abund = max(a_abund.values())
		all_alpha_abund += [max_mean_abund]
		b_prev = numpy.median(a_prev.values())

		centroid_prev_abund[centroid]['a_prev'] = a_prev
		centroid_prev_abund[centroid]['abund'] = max_mean_abund #a_abund
		centroid_prev_abund[centroid]['beta_prevalence'] = b_prev
	
	dict_to_print = {}
	for centroid in centroid_prev_abund:
		dict_to_print[centroid] = {'beta_prevalence': centroid_prev_abund[centroid]['beta_prevalence'], \
								   'mean_abundance': centroid_prev_abund[centroid]['abund']}
		for niche in centroid_prev_abund[centroid]['a_prev']:
			dict_to_print[centroid]['alpha_prevalence_'+niche] = centroid_prev_abund[centroid]['a_prev'][niche]
	write_prev_abund_matrix(dict_to_print, centroid_prev_abund_file_path)
	
	return [centroid_prev_abund, all_alpha_prev, all_alpha_abund]


def get_important_niche_centroids(centroid_prev_abund, all_alpha_prev, all_alpha_abund, tshld, output_folder):
	'''Returns the dict of important gene centroids [>= 10th percentile of alpha_prevalence and mean abundance]

	Input:	centroid_prev_abund = {centroid: {'abund': mean abundance, 'prev': prevalence}}
			all_alpha_preva = {niche: List of all observed gene centroid alpha prevalence values (>0) across samples within each niche}
			all_mean_abund = [List of all calculated mean gene centroid abundance across samples]
			output_folder = Location of the results folder

	Output: imp_centroids = {centroid: {'mean_abundance': mean abundance, 
										'beta_prevalence': median of alpha prevalences observed in each niche, 
										'alpha_prevalence_NICHEX': alpha_prevalence for niche X, ...}}'''
	logger.debug('get_important_niche_centroids')

	imp_centroid_prev_abund_file_path = basename+'_imp_centroid_prev_abund.txt'
	tshld_prev = {}
	
	tshld_abund = numpy.percentile(all_alpha_abund, tshld[0])
	
	logger.debug('get_important_niche_centroids: tshld_abund:'+str(tshld_abund))
	
	for niche in all_alpha_prev:
		tshld_prev[niche] = 2*numpy.std(all_alpha_prev[niche])/numpy.sqrt(len(all_alpha_prev[niche]))#, tshlds['prev']) 
	
		logger.debug('get_important_niche_centroids: tshld_prev:'+str(niche)+':'+str(tshld_prev[niche]))
	
	imp_centroids = {}

	for centroid in centroid_prev_abund:
		abund_check = centroid_prev_abund[centroid]['abund'] >= tshld_abund
		prev_check = sum([centroid_prev_abund[centroid]['a_prev'][niche]- tshld_prev[niche]> tshld[1] for niche in centroid_prev_abund[centroid]['a_prev']])
		
		if abund_check and prev_check:
			imp_centroids[centroid]={'mean_abundance': centroid_prev_abund[centroid]['abund'], \
									 'beta_prevalence': centroid_prev_abund[centroid]['beta_prevalence']}
			for niche in centroid_prev_abund[centroid]['a_prev']:
				imp_centroids[centroid]['alpha_prevalence_' + niche] = centroid_prev_abund[centroid]['a_prev'][niche]
	
	write_prev_abund_matrix(imp_centroids, output_folder + '/' + imp_centroid_prev_abund_file_path)
	
	return imp_centroids


def get_important_centroids(centroid_prev_abund, all_prevalence, all_abund, tshld, output_folder):
	'''Returns the dict of important gene centroids [value-2SE(prevalence and abundance) >0.1]

	Input:	centroid_prev_abund = {centroid: {'abund': mean abundance, 'prev': prevalence}}
			all_prevalence = [List of all observed gene centroid prevalence values (>0) across samples]
			all_mean_abund = [List of all calculated mean gene centroid abundance across samples]
			output_folder = Location of the results folder

	Output: imp_centroids = {centroid: {'abund': mean abundance, 'prev': prevalence}}'''
	logger.debug('get_important_centroids')
	imp_centroid_prev_abund_file_path = basename+'_imp_centroid_prev_abund.txt'
	tshld_prev = 2*numpy.std(numpy.array(all_prevalence))/numpy.sqrt(len(all_prevalence))

	tshld_abund = numpy.percentile(all_abund, tshld[0])
	
	imp_centroids = {}

	for centroid in centroid_prev_abund:
		abund_check = centroid_prev_abund[centroid]['abund'] >= tshld_abund
		prev_check = centroid_prev_abund[centroid]['prev'] - tshld_prev > tshld[1]
		if abund_check and prev_check:
			imp_centroids[centroid]={'mean_abundance': centroid_prev_abund[centroid]['abund'], \
									 'prevalence': centroid_prev_abund[centroid]['prev']}
	
	write_prev_abund_matrix(imp_centroids, output_folder + '/' + imp_centroid_prev_abund_file_path)
	
	return imp_centroids

def write_prev_abund_matrix(centroid_prev_abund, out_file):
	'''Writes the centroids prevalence and abundance information in text file

	Input: centroid_prev_abund = {centroids: {'abund': mean abundance, ...}}
		   out_file = output_filename

	Output: Writes the centroids dictionary to the output_filename'''
	logger.debug('write_prev_abund_matrix')
	

	keys = []
	for i in centroid_prev_abund:
		keys = centroid_prev_abund[i].keys()
		break
	
	with open(out_file,'w') as foo:
		foo.writelines(['Centroids\t' + str.join('\t', keys) + '\n'])
		for centroid in centroid_prev_abund:
			foo.writelines([str.join('\t', [centroid] + [str(centroid_prev_abund[centroid][key]) for key in keys]) + '\n'])

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--input_table', help='REQUIRED: Gene abundance table with metadata', required=True)
	parser.add_argument('-o','--output_folder', help='Folder containing results', default=False)
	parser.add_argument('--gene_catalog', default=False, help='GENE CATALOG')
	parser.add_argument('--uc', default=False, help='UCLUST file containg centroids and clustered genes')
	parser.add_argument('--usearch', default=False, help='Path to USEARCH') #add to be in path?
	parser.add_argument('--vsearch', default=False, help='Path to VSEARCH') #add to be in path?
	parser.add_argument('--basename', help='BASENAME for all the output files')
	parser.add_argument('--log_level',default='DEBUG', help='Choices: [DEBUG, INFO, WARNING, ERROR, CRITICAL]')
	parser.add_argument('--threads', default=1, help='Number of threads')
	parser.add_argument('--tshld_abund', default=75, help='[X] Percentile Cutoff for Abundance; Default=75th')
	parser.add_argument('--tshld_prev', default=0.1, help='Threshold: val-2*SE > tshld_prev')

	args = parser.parse_args()
	nprocesses = args.threads

	tshld = [float(args.tshld_abund), float(args.tshld_prev)]

	basename = args.basename
	input_table = args.input_table
	uclust_file = args.uc
	gene_catalog = args.gene_catalog

	if not basename:
		basename = input_table.split('.')[0].split('/')[-1]

	if not args.output_folder:
		output_folder = basename
	else:
		output_folder = args.output_folder

	temp_folder = output_folder+'/'+basename+'_temp'

	utilities.create_folders([temp_folder, output_folder])

	log_file = output_folder+'/'+basename+'.log'
	logging.basicConfig(filename=log_file, \
						format='%(asctime)s - %(name)s - %(levelname)s: %(message)s', \
						level=getattr(logging, args.log_level), \
						filemode='w', \
						datefmt='%m/%d/%Y %I:%M:%S %p')

	[uniref_dm, gi_dm, metadata]= read_gene_table(input_table)
	all_centroids = get_centroids(uniref_dm, gi_dm, usearch_folder, uclust_file, gene_catalog, nprocesses)
	[centroids_data_matrix, centroids_list] = get_centroids_table(all_centroids, metadata)
	[centroid_prev_abund, all_prevalence, all_mean_abund, niche_flag] = get_prevalence_abundance(centroids_data_matrix, centroids_list, metadata)

	if niche_flag:
		imp_centroids = get_important_niche_centroids(centroid_prev_abund, all_prevalence, all_mean_abund, tshld, output_folder)
	else:
		imp_centroids = get_important_centroids(centroid_prev_abund, all_prevalence, all_mean_abund, tshld, output_folder)
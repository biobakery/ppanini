import os
import sys
import pdb
import re
import argparse
import numpy

from src import create_fastas
from src import create_annotations

basename = ''

def read_gene_table(gene_table_fname):
	'''Returns the different elements from the gene table

	Input: gene_table_fname = Filename of the gene_table

	Output: metadata = [metadata strings]; Rows with # as first character in table
			uniref_gis = {UniRef_XYZ: [list of gene ids]}
			gis_unannotated = {sample_id: [gene ids]}
			gene_ids = [List of all gene ids]
			data_matrix = The abundance table [[0, 0,.],[0, 0,.],...]'''

	gene_table = open(gene_table_fname) 
	gtab_lines = gene_table.readlines()
	metadata, gene_ids, data_matrix = [], [], []
	gis_unannotated, uniref_gis = {}, {} # Cluster IDs: [Gene IDs]
	
	for line in gtab_lines:
		#Metadata containing Sample names, Niche specifics and optionally fasta file locations
		if line.startswith('#'):
			metadata += [line]
		samples = [re.sub('[\r\t\n]','',i) for i in metadata[-1].split('\t')[1:]]
		
		if not line.startswith('#'):
			split_i = line.split('\t')
			annot = split_i[0].split('|') #geneID column split 
			u90_annot = [i for i in annot if 'UniRef90' in i][0]
			gene_ids += [annot[0]]
		
			data_row = [float(i) for i in split_i[1:]]
			sample_inds = [i for i, val in enumerate(data_row) if val > 0]
			data_matrix += [data_row]
			#Add unknown UniRef gene ids to a dict to be processed for clustering later {Sample: [gids]}
			if 'UniRef90_unknown' == u90_annot:
				for i in sample_inds:
					if samples[i] not in gis_unannotated:
						gis_unannotated[samples[i]] = [annot[0]]
					else:
						gis_unannotated[samples[i]] += [annot[0]]
			else:

				try:
					#Add gene id to the list of uniref90 cluster id genes
					uniref_gis[u90_annot] += [annot[0]] 
				except KeyError:
					#Initiate a list of gene ids that belong to a specific UniRef90 ID
					uniref_gis[u90_annot] = [annot[0]] 
					
	return [metadata, uniref_gis, gis_unannotated, gene_ids, data_matrix] 

def extract_fasta_names(metadata):
	'''Returns the dict of fasta files corresp. to each sample from metadata
	
	Input:	metadata = [metadata strings]; Rows with # as first character in table
	Output: location = {sample_id : path_to_fasta_file}'''

	location = {}
	
	[line, ind] = is_present(metadata, '#FASTAS')

	samples = metadata[-1].split('\t')[1:]
	
	if line:
		split_i = line.split('\t')
		for i, val in enumerate(split_i[1:]):
			location[re.sub('[\t\r\n]', '', samples[i])] = re.sub('[\t\n\r]', '', val)
	else:
		raise Exception("Missing #FASTAS metadata for sample names!")

	return location


def get_centroids(gis_unannotated, metadata, usearch_folder, uniref_gis, uclust_file):
	'''Returns the dict of all centroids containing clusters of gene IDs

	Input:	gis_unannotated = {sample_id: [gene ids]}
			metadata = [metadata strings]; Rows with # as first character in table
			uniref_gis = {UniRef_XYZ: [list of gene ids]}
			gis_unannotated = {sample_id: [gene ids]}
			usearch_folder = Location of the USEARCH program

	Output: all_centroids = {gene_centroid : [List of gene ids]}'''


	centroids_fasta = {}

	location = extract_fasta_names(metadata)
	if not uclust_file:
		for sample in gis_unannotated:
			genes = gis_unannotated[sample]
			file_fasta = create_fastas.read_fasta(location[sample])
			for gid in genes:
				try:
					centroids_fasta[gid] = file_fasta[gid]
				except KeyError:
					raise Exception('Unannotated gene '+gene+' not found in sample FASTA: '+location[sample])

		centroid_gis = get_clusters(centroids_fasta, usearch_folder) #all UniRef90_unknowns are clustered across samples
	else:
		pre_centroid_gis = create_annotations.get_clusters_dict(gene_centroid_clusters_file_path)
		uniref_ids = []
		centroids_gis = {}
		for uid in uniref_gis:
			uniref_ids += uniref_gis[uid]
		
		for i in pre_centroid_gis:
			#assuming that if cluster not in UniRef, none of the cluster members are either?!!!!!
			if i not in uniref_ids: 
				centroids_gis[i] = pre_centroid_gis[i]

	for sample in gis_unannotated:
		for gi in gis_unannotated[sample]:
			if gi not in centroid_gis:
				centroid_gis[gi] = [gi]
	
	for gi in uniref_gis:
		centroid_gis[gi] = uniref_gis[gi]
	
	return centroid_gis


def get_clusters(centroids_fasta, usearch_folder): #ONLY FOR THE UNIREF UNANNOTATED
	'''Returns the dict of unannotated gene centroids containing clusters of genes at 90% similarity
	Input:	centroids_fasta = {UniRef90_unknown: [List of gene ids]}
			usearch_folder = Location of the USEARCH program
	Output: centroids_gis = {gene_centroid: [List of genes in the cluster]}
			'''
	allgenes_file_path = 'data_files/'+basename+'_centroids_for_clustering.fasta'
	gene_centroids_file_path = 'data_files/'+basename+'_centroids.fasta'
	gene_centroid_clusters_file_path = 'data_files/'+basename+'_clusters.uc'

	create_fastas.write_fasta(centroids_fasta, allgenes_file_path) #ensures all clustering happens to FAAs

	create_annotations.run_uclust(usearch_folder, allgenes_file_path, gene_centroids_file_path, gene_centroid_clusters_file_path, 0.9)
	
	centroid_gis = create_annotations.get_clusters_dict(gene_centroid_clusters_file_path)

	return centroid_gis


def get_centroids_table(gene_ids, all_centroids, data_matrix, metadata):
	'''Returns data matrix containing gene centroids and abundance per sample
	Input:	metadata = [metadata strings]; Rows with # as first character in table
			gene_ids = [List of all gene ids]
			data_matrix = The abundance table [[0, 0,.],[0, 0,.],...]
			all_centroids = {gene_centroid : [List of gene ids]}
	Output: centroids_data_matrix = {gene_centroid: [Gene centroid abundance across samples]}
			'''
	
	centroids_data_matrix = {}
	
	for centroid in all_centroids:
		centroid_row = [[data_matrix[gene_ids.index(gene)]] for gene in all_centroids[centroid]]
		centroids_data_matrix[centroid] = sum(numpy.array(centroid_row))
		
		if len(centroids_data_matrix[centroid]) == 1: 
			#raise Exception('length of data row is one?')#if only one member in cluster
			centroids_data_matrix[centroid] = [i for i in centroids_data_matrix[centroid][0]]
		else:
			centroids_data_matrix[centroid] = list(centroids_data_matrix[centroid])#[0])
	
	gene_centroids_table_file_path = 'data_files/'+basename+'_gene_centroids_table.txt'
	
	foo = open(gene_centroids_table_file_path,'w')	

	foo.writelines(metadata[:-1])
	foo.writelines([str.join('\t', ['Centroids'] + metadata[-1].split('\t')[1:]) +  '\n'])
	for centroid in centroids_data_matrix:
		foo.writelines([str.join('\t', [centroid] + [str(i) for i in centroids_data_matrix[centroid]]) + '\n'])
	foo.close()
	
	return centroids_data_matrix 

def is_present(metadata, meta_type):
	'''Returns True if meta_type is present in metadata extracted from mappert_file'''
	line = []
	ind = []
	for i, val in enumerate(metadata):
		if val.upper().startswith(meta_type):
			line = val
			ind = i
			break
	return [line, ind]

def get_prevalence_abundance(centroids_data_matrix, metadata):
	'''Returns the dict of centroids with their prevalence and abundance
	Input:	centroids_data_matrix = {gene_centroid: [Gene centroid abundance across samples]}
			metadata = [metadata strings]; Rows with # as first character in table
	Output: centroid_prev_abund = {centroid: {'abund': mean abundance, 'prev': prevalence}}
			all_prevalence = [List of all observed gene centroid prevalence values (>0) across samples]
			all_mean_abund = [List of all calculated mean gene centroid abundance across samples]
			flag = True (if NICHE PRESENT) or False(if NICHE ABSENT)'''

	centroid_prev_abund_file_path = 'data_files/'+basename+'_centroid_prev_abund.txt'
	
	[line, ind] = is_present(metadata, '#NICHE')
	
	if line:
		[centroid_prev_abund, all_alpha_prev, all_mean_abund] = get_niche_prevalence_abundance (centroids_data_matrix, metadata, line, ind)
		return [centroid_prev_abund, all_alpha_prev, all_mean_abund, True]
	else:
		centroid_prev_abund = {}
		all_prevalence = [] 
		all_mean_abund = []
		
		for centroid in centroids_data_matrix:
			#abund only where the gene is present in sample
			centroid_prev_abund[centroid] = {'abund': numpy.mean(numpy.array([i for i in  centroids_data_matrix[centroid] if i > 0])), \
											 'prev': sum(numpy.array(centroids_data_matrix[centroid]) > 0)}
			
			all_prevalence += [centroid_prev_abund[centroid]['prev']]
			all_mean_abund += [centroid_prev_abund[centroid]['abund']]
		
		write_prev_abund_matrix(centroid_prev_abund, centroid_prev_abund_file_path)
		
		return [centroid_prev_abund, all_prevalence, all_mean_abund, False]

def get_niche_prevalence_abundance(centroids_data_matrix, metadata, line, ind):
	'''Returns the dict of centroids with their prevalence and abundance
	Input:	centroids_data_matrix = {gene_centroid: [Gene centroid abundance across samples]}
			metadata = [metadata strings]; Rows with # as first character in table
	Output: centroid_prev_abund = {centroid: {'abund': float_mean abundance, 
											  'a_prev': {'NICHE**': prevalence of centroid across samples within niche, ...}},
											  'b_prev': float_median of alpha prevalence within niches for centroid}}
			all_alpha_prev = {niche: List of all observed gene centroid alpha prevalence values (>0) across samples within each niche}
			all_mean_abund = [List of all calculated mean gene centroid abundance across samples]'''

	centroid_prev_abund_file_path = 'data_files/'+basename+'_centroid_prev_abund.txt'
	
	niches = {}
	
	split_i = [re.sub('[\r\t\n]', '', i) for i in line.split('\t')[1:]]
	for i, val in enumerate(split_i):
		if val not in niches:
			niches[val] = [i]
		else:
			niches[val] += [i]
	
	niches_label = niches.keys()
	centroid_prev_abund = {}

	all_mean_abund = []
	all_alpha_prev = {}

	for niche in niches_label:
		all_alpha_prev[niche] = []

	for centroid in centroids_data_matrix:
		centroid_prev_abund[centroid] = {'abund': numpy.mean(numpy.array([i for i in centroids_data_matrix[centroid] if i > 0]))}
		
		a_prev = {}
		for niche in niches_label:
			a_prev[niche] = 0
 			for i in niches[niche]:
				a_prev[niche] += int(centroids_data_matrix[centroid][i] > 0)
			a_prev[niche] = float(a_prev[niche])/float(len(niches[niche]))
			all_alpha_prev[niche] += [a_prev[niche]]
		centroid_prev_abund[centroid]['a_prev'] = a_prev

		b_prev = numpy.mean(a_prev.values())
		centroid_prev_abund[centroid]['b_prev'] = b_prev
		all_mean_abund += [centroid_prev_abund[centroid]['abund']]

	write_prev_abund_matrix(centroid_prev_abund, centroid_prev_abund_file_path)
	
	return [centroid_prev_abund, all_alpha_prev, all_mean_abund]


def get_important_niche_centroids(centroid_prev_abund, all_alpha_prev, all_mean_abund, output_folder):
	'''Returns the dict of important gene centroids [>= 10th percentile of alpha_prevalence and mean abundance]
	Input:	centroid_prev_abund = {centroid: {'abund': mean abundance, 'prev': prevalence}}
			all_alpha_preva = {niche: List of all observed gene centroid alpha prevalence values (>0) across samples within each niche}
			all_mean_abund = [List of all calculated mean gene centroid abundance across samples]
			output_folder = Location of the results folder
	Output: imp_centroids = {centroid: {'mean_abundance': mean abundance, 
										'beta_prevalence': median of alpha prevalences observed in each niche, 
										'alpha_prevalence_NICHEX': alpha_prevalence for niche X, ...}}'''
	imp_centroid_prev_abund_file_path = 'imp_centroid_prev_abund.txt'
	tshld_prev = {}
	#not good for small sample space? 10th percentile; need a better metric? z-score?
	for niche in all_alpha_prev:
		tshld_prev[niche] = numpy.percentile(all_alpha_prev[niche], 10) 
	tshld_abund = numpy.percentile(all_mean_abund, 10)

	imp_centroids = {}

	for centroid in centroid_prev_abund:
		abund_check = centroid_prev_abund[centroid]['abund'] >= tshld_abund
		#If Alpha-prevalence of the centroid is higher than the niche-specific threshold in ANY of the niches;
		#Note: SUM is used to implement the OR functionality
		prev_check = sum([centroid_prev_abund[centroid]['a_prev'][niche]>=tshld_prev[niche] for niche in centroid_prev_abund[centroid]['a_prev']])
		
		if abund_check and prev_check:
			imp_centroids[centroid]={'mean_abundance': centroid_prev_abund[centroid]['abund'], \
									 'beta_prevalence': centroid_prev_abund[centroid]['b_prev']}
			for niche in centroid_prev_abund[centroid]['a_prev']:
				imp_centroids[centroid]['alpha_prevalence_' + niche] = centroid_prev_abund[centroid]['a_prev'][niche]
	
	write_prev_abund_matrix(imp_centroids, output_folder + '/' + imp_centroid_prev_abund_file_path)
	
	return imp_centroids


def get_important_centroids(centroid_prev_abund, all_prevalence, all_mean_abund, output_folder):
	'''Returns the dict of important gene centroids [>= 10th percentile of prevalence and abundance]
	Input:	centroid_prev_abund = {centroid: {'abund': mean abundance, 'prev': prevalence}}
			all_prevalence = [List of all observed gene centroid prevalence values (>0) across samples]
			all_mean_abund = [List of all calculated mean gene centroid abundance across samples]
			output_folder = Location of the results folder
	Output: imp_centroids = {centroid: {'abund': mean abundance, 'prev': prevalence}}'''
	
	imp_centroid_prev_abund_file_path = 'imp_centroid_prev_abund.txt'
	
	tshld_prev = numpy.percentile(all_prevalence, 10)
	tshld_abund = numpy.percentile(all_mean_abund, 10)

	imp_centroids = {}

	for centroid in centroid_prev_abund:
		abund_check = centroid_prev_abund[centroid]['abund'] >= tshld_abund
		prev_check = centroid_prev_abund[centroid]['prev'] >= tshld_prev
		if abund_check and prev_check:
			imp_centroids[centroid]={'mean_abundance': centroid_prev_abund[centroid]['abund'], \
									 'prevalence': centroid_prev_abund[centroid]['prev']}
	
	write_prev_abund_matrix(imp_centroids, output_folder + '/' + imp_centroid_prev_abund_file_path)
	
	return imp_centroids

def write_prev_abund_matrix(centroid_prev_abund, out_file):
	'''Writes the centroids prevalence and abundance information in text file'''

	foo = open(out_file,'w')
	#Assumes the dictionary structure to be consistent and fixed
	#dict_X = {key1: {subkey1: [value]}}
	keys = centroid_prev_abund.values()[0].keys()
	
	foo.writelines(['Centroids\t' + str.join('\t', keys) + '\n'])
	for centroid in centroid_prev_abund:
		foo.writelines([str.join('\t', [centroid] + [str(centroid_prev_abund[centroid][key]) for key in keys]) + '\n'])
	foo.close()

def create_folders(list_folders):

	for fname in list_folders:
		try:
			os.mkdir(fname)
		except:
			pass

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--input_table', help='Gene abundance table with metadata', required=True)
	parser.add_argument('-u','--usearch_folder', nargs = '?' , help='Path for USEARCH program; if not provided, assumed to be in path')
	parser.add_argument('-o','--output_folder', help='Folder containing results', default='.')
	parser.add_argument('--uclust_file', help='File containing UCLUST results')
	parser.add_argument('--basename', help='BASENAME for all the output files')

	args = parser.parse_args()

	create_folders(['data_files', args.output_folder])

	basename = args.basename
	input_table = args.input_table
	uclust_file = args.uclust_file

	if not basename:
		basename = input_table.split('.')[0].split('/')[-1]
	
	[metadata, uniref_gis, gis_unannotated, gene_ids, data_matrix] = read_gene_table(input_table)
	all_centroids = get_centroids(gis_unannotated, metadata, args.usearch_folder, uniref_gis, uclust_file)
	centroids_data_matrix = get_centroids_table(gene_ids, all_centroids, data_matrix, metadata)
	[centroid_prev_abund, all_prevalence, all_mean_abund, flag] = get_prevalence_abundance(centroids_data_matrix, metadata)

	if flag:
		#Niche-classification present
		imp_centroids = get_important_niche_centroids(centroid_prev_abund, all_prevalence, all_mean_abund, args.output_folder)
	else:
		imp_centroids = get_important_centroids(centroid_prev_abund, all_prevalence, all_mean_abund, args.output_folder)

	
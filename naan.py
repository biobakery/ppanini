import os
import sys
import pdb
import re
import argparse
import numpy

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
			u90_annot = [i for i in annot if 'UniRef90' in i]
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
		
				if u90_annot in uniref_gis: ###CHECK AND EDIT
					#Add gene id to the list of uniref90 cluster id genes
					uniref_gis[u90_annot] += [annot[0]] 
				else:
					#Initiate a list of gene ids that belong to a specific UniRef90 ID
					uniref_gis[u90_annot] = [annot[0]] 
	
	return [metadata, uniref_gis, gis_unannotated, gene_ids, data_matrix] 

def extract_fasta_names(metadata, fasta_folder):
	'''Returns the dict of fasta files corresp. to each sample from metadata
	
	Input:	metadata = [metadata strings]; Rows with # as first character in table
			fasta_folder = Location of the fasta files
	Output: location = {sample_id : path_to_fasta_file}'''

	location = {}
	
	[line, ind] = is_present(metadata, '#FASTA')

	samples = metadata[-1].split('\t')[1:]
	
	if line:
		split_i = line.split('\t')
		for i, val in enumerate(split_i[1:]):
			location[re.sub('[\t\r\n]', '', samples[i])] = fasta_folder + '/' + re.sub('[\t\n\r]', '', val)
	else:
		raise Exception("Missing #FASTA metadata for sample names!")

	return location


def get_centroids(gis_unannotated, fasta_folder, metadata, usearch_folder, uniref_gis):
	'''Returns the dict of all centroids containing clusters of gene IDs

	Input:	gis_unannotated = {sample_id: [gene ids]}
			metadata = [metadata strings]; Rows with # as first character in table
			uniref_gis = {UniRef_XYZ: [list of gene ids]}
			gis_unannotated = {sample_id: [gene ids]}
			fasta_folder = Location of the fasta files
			usearch_folder = Location of the USEARCH program

	Output: all_centroids = {gene_centroid : [List of gene ids]}'''

	
	centroids_fasta = []
	
	location = extract_fasta_names(metadata, fasta_folder)

	for sample in gis_unannotated:
		genes = gis_unannotated[sample]
		foo = open(location[sample],'r')
		check = False
		genes_done = 0
		for line in foo.readlines():
			if '>' in line and re.sub('[\t\r\n]', '', line[1:]).strip() in genes:
				check = True
				genes_done += 1
				centroids_fasta += [line]
			elif '>' not in line and check:
				centroids_fasta += [line]
			elif '>' in line and re.sub('[\t\r\n]', '', line[1:]).strip() not in genes:
				check = False
			if not check and genes_done == len(genes):
				break
		foo.close()
	centroid_gis = get_clusters(centroids_fasta, usearch_folder)
	
	all_centroids = {}

	for sample in gis_unannotated:
		for gi in gis_unannotated[sample]:
			if gi not in centroid_gis:
				centroid_gis[gi] = [gi]
		all_centroids = dict(centroid_gis.items() + uniref_gis.items())

	
	return all_centroids 


def get_clusters(centroids_fasta, usearch_folder): #ONLY FOR THE UNIREF UNANNOTATED
	'''Returns the dict of unannotated gene centroids containing clusters of genes at 90% similarity
	Input:	centroids_fasta = {UniRef90_unknown: [List of gene ids]}
			usearch_folder = Location of the USEARCH program
	Output: centroids_gis = {gene_centroid: [List of genes in the cluster]}
			'''
	allgenes_file_path = 'data_files/centroids_for_clustering.fasta'
	gene_centroids_file_path = 'data_files/centroids.fasta'
	gene_centroid_clusters_file_path = 'data_files/clusters.uc'

	foo = open(allgenes_file_path,'w')
	foo.writelines(centroids_fasta)
	foo.close()

	out_clust = os.system(usearch_folder + \
		'/usearch -cluster_fast ' + allgenes_file_path +' \
		  		  -id 0.9 \
		          -centroids '+ gene_centroids_file_path + ' \
		          -uc ' + gene_centroid_clusters_file_path)

	cluster_txt = os.popen('grep -w H ' + gene_centroid_clusters_file_path)
	centroid_gis = {}
	for line in cluster_txt.xreadlines():
		split_i = [re.sub('[\r\t\n]', '', i) for i in line.split('\t')]
		if split_i[-1] in centroid_gis:
			centroid_gis[split_i[-1]] += [split_i[-2]]
		else:
			centroid_gis[split_i[-1]] = [split_i[-2], split_i[-1]]
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
		
		if len(centroids_data_matrix[centroid]) == 1: #if only one member in cluster
			centroids_data_matrix[centroid] = [i for i in centroids_data_matrix[centroid][0]]
		else:
			centroids_data_matrix[centroid] = list(centroids_data_matrix[centroid][0])
	
	gene_centroids_table_file_path = 'data_files/gene_centroids_table.txt'
	
	foo = open(gene_centroids_table_file_path,'w')	

	foo.writelines(metadata[:-1])
	foo.writelines([str.join('\t', ['Centroids'] + metadata[-1].split('\t')[1:]) +  '\n'])
	for centroid in centroids_data_matrix:
		foo.writelines([str.join('\t', [centroid] + [str(i) for i in centroids_data_matrix[centroid]]) + '\n'])
	foo.close()
	
	return centroids_data_matrix 

def is_present(metadata, meta_type):
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

	centroid_prev_abund_file_path = 'data_files/centroid_prev_abund.txt'
	
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

	centroid_prev_abund_file_path = 'data_files/centroid_prev_abund.txt'
	
	niches = {}
	
	split_i = [re.sub('[\r\t\n]', '', i) for i in line.split('\t')[2:]]
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
										'alpha_prevalence_NICHE***': alpha_prevalence for niche X, ...}}'''
	imp_centroid_prev_abund_file_path = 'imp_centroid_prev_abund.txt'
	tshld_prev = {}
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

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--input_table', help='Gene abundance table with metadata')
	parser.add_argument('-f','--fasta_folder', help='Folder containing fasta files')
	parser.add_argument('-u','--usearch_folder', nargs = '?' , help='Path for USEARCH program')
	parser.add_argument('-o','--output_folder', help='Folder containing results')
	args = parser.parse_args()

	try:
		os.mkdir('data_files')
	except:
		pass

	try:
		os.mkdir(args.output_folder)
	except:
		pass
	
	[metadata, uniref_gis, gis_unannotated, gene_ids, data_matrix] = read_gene_table(args.input_table)
	all_centroids = get_centroids(gis_unannotated, args.fasta_folder, metadata, args.usearch_folder, uniref_gis)
	centroids_data_matrix = get_centroids_table(gene_ids, all_centroids, data_matrix, metadata)
	[centroid_prev_abund, all_prevalence, all_mean_abund, flag] = get_prevalence_abundance(centroids_data_matrix, metadata)

	if flag:
		#Niche-classification present
		imp_centroids = get_important_niche_centroids(centroid_prev_abund, all_prevalence, all_mean_abund, args.output_folder)
	else:
		imp_centroids = get_important_centroids(centroid_prev_abund, all_prevalence, all_mean_abund, args.output_folder)

	
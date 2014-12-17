import os
import sys
import pdb
import re
import argparse
import subprocess
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
		if line[0] == '#':
			metadata += [line]
		samples = [re.sub('[\r\t\n]','',i) for i in metadata[-1].split('\t')[2:]]
		
		if not line[0] == '#':
			split_i = line.split('\t')
			gene_ids += [split_i[0]]
		
			data_row = [float(i) for i in split_i[2:]]
			sample_inds = [i for i , val in enumerate(data_row) if val>0]
			data_matrix += [data_row]
		
			#Add unknown UniRef gene ids to a dict to be processed for clustering later {Sample: [gids]}
			if 'UniRef_unknown' == split_i[1]:
		
				for i in sample_inds:
					if samples[i] not in gis_unannotated:
						gis_unannotated[samples[i]] = [split_i[0]]
					else:
						gis_unannotated[samples[i]] += [split_i[0]]
			else:
		
				if split_i[1] in uniref_gis:
					#Add gene id to the list of uniref90 cluster id genes
					uniref_gis[split_i[1]] += [split_i[0]] 
				else:
					#Initiate a list of gene ids that belong to a specific UniRef90 ID
					uniref_gis[split_i[1]] = [split_i[0]] 
	
	return [metadata, uniref_gis, gis_unannotated, gene_ids, data_matrix] 

def extract_fasta_names(metadata, fasta_folder):
	'''Returns the dict of fasta files corresp. to each sample from metadata
	
	Input:	metadata = [metadata strings]; Rows with # as first character in table
			fasta_folder = Location of the fasta files
	Output: location = {sample_id : path_to_fasta_file}'''

	location = {}
	line = []
	for i in metadata:
		if '#FASTA' in i:
			line = i
			break
	samples = metadata[-1].split('\t')[2:]
	
	if line:
		split_i = line.split('\t')
		for i in range(len(split_i[2:])):
			location[re.sub('[\t\r\n]', '', samples[i])] = fasta_folder+'/'+re.sub('[\t\n\r]','',split_i[2:][i])
	else:
		raise Exception("Missing #FASTA metadata for sample names!")

	return location


def get_centroids(gis_unannotated, fasta_folder, metadata, uclust_folder, uniref_gis):
	'''Returns the dict of all centroids containing clusters of gene IDs

	Input:	gis_unannotated = {sample_id: [gene ids]}
			metadata = [metadata strings]; Rows with # as first character in table
			uniref_gis = {UniRef_XYZ: [list of gene ids]}
			gis_unannotated = {sample_id: [gene ids]}
			fasta_folder = Location of the fasta files
			uclust_folder = Location of the USEARCH program

	Output: all_centroids = {gene_centroid : [List of gene ids]}'''

	#samples = metadata[-1][3:] Pending functionality--ability to point to fasta files??
	centroids_fasta = []
	
	location = extract_fasta_names(metadata, fasta_folder)

	for sample in gis_unannotated:
		genes = gis_unannotated[sample]
		foo = open(location[sample],'r')
		check = False
		genes_done = 0
		for line in foo.readlines():
			if '>' in line and re.sub('[\t\r\n]','', line[1:]).strip() in genes:
				check = True
				genes_done += 1
				centroids_fasta += [line]
			elif '>' not in line and check:
				centroids_fasta += [line]
			elif '>' in line and not re.sub('[\t\r\n]','', line[1:]).strip() in genes:
				check = False
			if not check and genes_done == len(genes):
				break
		foo.close()
	centroid_gis = get_clusters(centroids_fasta, uclust_folder)
	for sample in gis_unannotated:
		for gi in gis_unannotated[sample]:
			if not gi in centroid_gis:
				centroid_gis[gi] = [gi]
		all_centroids = dict(centroid_gis.items()+ uniref_gis.items())

	return all_centroids 


def get_clusters(centroids_fasta, uclust_folder):
	'''Returns the dict of unannotated gene centroids containing clusters of genes at 90% similarity
	Input:	centroids_fasta = {UniRef90_unknown: [List of gene ids]}
			uclust_folder = Location of the USEARCH program
	Output: centroids_gis = {gene_centroid: [List of genes in the cluster]}
			'''

	foo = open('tmp/centroids_for_clustering.fasta','w')
	foo.writelines(centroids_fasta)
	foo.close()

	out_clust = os.system(uclust_folder + \
		'/usearch -cluster_fast tmp/centroids_for_clustering.fasta -id 0.9 \
		          -centroids tmp/centroids.fasta -uc tmp/clusters.uc')

	cluster_txt = os.popen('grep -w H tmp/clusters.uc')
	centroid_gis = {}
	for line in cluster_txt.xreadlines():
		split_i = [re.sub('[\r\t\n]', '',i) for i in line.split('\t')]
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
		centroids_data_matrix[centroid] = sum(numpy.array([[data_matrix[gene_ids.index(gene)]] \
												for gene in all_centroids[centroid]]))
		
		if len(centroids_data_matrix[centroid]) == 1:
			centroids_data_matrix[centroid] = [i for i in centroids_data_matrix[centroid][0]]
		else:
			centroids_data_matrix[centroid] = list(centroids_data_matrix[centroid][0])
	
	foo = open('tmp/gene_centroids_table.txt','w')	
	foo.writelines(metadata[:-1])
	foo.writelines([str.join('\t', ['Centroids']+metadata[-1].split('\t')[2:])+'\n'])
	for centroid in centroids_data_matrix:
		foo.writelines([str.join('\t', [centroid]+[str(i) for i in centroids_data_matrix[centroid]])+'\n'])
	foo.close()
	
	return centroids_data_matrix 

def get_prevalence_abundance(centroids_data_matrix, metadata):
	'''Returns the dict of centroids with their prevalence and abundance
	Input:	centroids_data_matrix = {gene_centroid: [Gene centroid abundance across samples]}
			metadata = [metadata strings]; Rows with # as first character in table
	Output: centroid_prev_abund = {centroid: {'abund': mean abundance, 'prev': prevalence}}
			all_prevalence = [List of all observed gene centroid prevalence values (>0) across samples]
			all_mean_abund = [List of all calculated mean gene centroid abundance across samples]'''
	
	#Niche specific Beta-prevalence?
	centroid_prev_abund = {}
	all_prevalence = []
	all_mean_abund = []
	
	for centroid in centroids_data_matrix:
		centroid_prev_abund[centroid] = {'abund':numpy.mean(numpy.array(centroids_data_matrix[centroid])), \
										 'prev': sum(numpy.array(centroids_data_matrix[centroid])>0)}
		
		all_prevalence += [centroid_prev_abund[centroid]['prev']]
		all_mean_abund += [centroid_prev_abund[centroid]['abund']]
	
	write_prev_abund_matrix(centroid_prev_abund, 'tmp/centroid_prev_abund.txt')
	
	return [centroid_prev_abund, all_prevalence, all_mean_abund]

def get_important_centroids(centroid_prev_abund, all_prevalence, all_mean_abund, output_folder):
	'''Returns the dict of important gene centroids [>= 10th percentile of prevalence and abundance]
	Input:	centroid_prev_abund = {centroid: {'abund': mean abundance, 'prev': prevalence}}
			all_prevalence = [List of all observed gene centroid prevalence values (>0) across samples]
			all_mean_abund = [List of all calculated mean gene centroid abundance across samples]
			output_folder = Location of the results folder
	Output: imp_centroids = {centroid: {'abund': mean abundance, 'prev': prevalence}}'''
	
	tshld_prev = numpy.percentile(all_prevalence, 10)
	tshld_abund = numpy.percentile(all_mean_abund, 10)

	imp_centroids = {}

	for centroid in centroid_prev_abund:
		if centroid_prev_abund[centroid]['abund'] >=tshld_abund and centroid_prev_abund[centroid]['prev'] >= tshld_prev:
			imp_centroids[centroid]={'abund': centroid_prev_abund[centroid]['abund'], \
									 'prev': centroid_prev_abund[centroid]['prev']}
	
	write_prev_abund_matrix(centroid_prev_abund, output_folder+'/imp_centroid_prev_abund.txt')
	
	return imp_centroids

def write_prev_abund_matrix(centroid_prev_abund, out_file):
	'''Writes the centroids prevalence and abundance information in text file'''

	foo = open(out_file,'w')
	foo.writelines(['Centroids\tAbundance\tPrevalence\n'])
	for centroid in centroid_prev_abund:
		foo.writelines([str.join('\t', [centroid, str(centroid_prev_abund[centroid]['abund']), \
												str(centroid_prev_abund[centroid]['prev'])])+'\n'])
	foo.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--input_table', help='Gene abundance table with metadata')
	parser.add_argument('-f','--fasta_folder', help='Folder containing fasta files')
	parser.add_argument('-u','--uclust_folder', nargs = '?' , help='Path for UCLUST program')
	parser.add_argument('-t','--type', help='Type of analysis Choose: [gene_table, reads_assemblies]')
	parser.add_argument('-o','--output_folder', help='Folder containing results')
	args = parser.parse_args()

	try:
		os.mkdir('tmp')
	except:
		pass

	try:
		os.mkdir(args.output_folder)
	except:
		pass
	
	[metadata, uniref_gis, gis_unannotated, gene_ids, data_matrix] = read_gene_table(args.input_table)
	all_centroids = get_centroids(gis_unannotated, args.fasta_folder, metadata, args.uclust_folder, uniref_gis)
	centroids_data_matrix = get_centroids_table(gene_ids, all_centroids, data_matrix, metadata)
	[centroid_prev_abund, all_prevalence, all_mean_abund] = get_prevalence_abundance(centroids_data_matrix, metadata)
	imp_centroids = get_important_centroids(centroid_prev_abund, all_prevalence, all_mean_abund, args.output_folder)

	
import os
import sys
import pdb
import re
import argparse
import subprocess
import numpy

def read_gene_table(gene_table_fname):
	gene_table = open(gene_table_fname) # Gene table containing UniRef annotations and gene abundance information.
	gtab_lines = gene_table.readlines()
	metadata = []
	uniref_gis = {} # Cluster IDs: [Gene IDs]
	gis_unannotated = {}
	gene_ids = []
	data_matrix = []
	for line in gtab_lines:
		if line[0] == '#':
			#Metadata containing Sample names, Niche specifics and optionally fasta file locations
			metadata += [line]
		samples = [re.sub('[\r\t\n]','',i) for i in metadata[-1].split('\t')[2:]]
		if not line[0] == '#':
			split_i = line.split('\t')
			gene_ids += [split_i[0]]
			#data_matrix contains the numeric matrix of genes vs. sample abundance
			data_row = [float(i) for i in split_i[3:]]
			sample_inds = [i for i , val in enumerate(data_row) if val>0]
			data_matrix += [data_row]
			if 'UniRef_unknown' == split_i[1]:
				#Add unknown UniRef gene ids to a dict to be processed for clustering later {Sample: [gids]}
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
	pdb.set_trace()
	return [metadata, uniref_gis, gis_unannotated, gene_ids, data_matrix] 

def get_centroids(gis_unannotated, fasta_folder, metadata, uclust_folder, uniref_gis):
	#samples = metadata[-1][3:] Pending functionality--ability to point to fasta files??
	fasta_folder = fasta_folder
	gis_unannotated = gis_unannotated
	centroids_fasta = []
	for sample in gis_unannotated:
		genes = gis_unannotated[sample]
		foo = open(fasta_folder+'/'+sample+'.with_fasta.gff3.faa','r')
		check = False
		genes_done = 0
		for line in foo.readlines():
			print 'enters'
			if '>' in line and re.sub('[\t\r\n]','', line[1:]).strip() in genes:
				print 'comes here'
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
	pdb.set_trace()
	centroid_gis = get_clusters(centroids_fasta, uclust_folder)

	for sample in gis_unannotated:
		for gi in gis_unannotated[sample]:
			if not gi in centroid_gis:
				centroid_gis[gi] = [gi]
		all_centroids = dict(centroid_gis.items()+ uniref_gis.items())

	return all_centroids 


def get_clusters(centroids_fasta, uclust_folder):
	try:
		os.mkdir('tmp')
	except:
		pdb.set_trace()
	foo = open('tmp/centroids_for_clustering.fasta','w')
	foo.writelines(centroids_fasta)
	foo.close()

	out_clust = os.system(uclust_folder+'/usearch -cluster_fast tmp/centroids_for_clustering.fasta -id 0.9 -centroids tmp/centroids.fasta -uc tmp/clusters.uc')
	cluster_txt = os.popen('grep -w H clusters.uc')
	centroid_gis = {}
	for line in cluster_txt.xreadlines():
		split_i = [re.sub('[\r\t\n]', '',i) for i in line.split('\t')]
		if split_i[-1] in centroid_gis:
			centroid_gis[split_i[-1]] += [split_i[-2]]
		else:
			centroid_gis[split_i[-1]] = [split_i[-2], split_i[-1]]
	return centroid_gis



def get_centroids_table(gene_ids, all_centroids, data_matrix, metadata):
	gene_ids = gene_ids #List of all genes
	centroids_data_matrix = {}
	for centroid in all_centroids:
		try:
			centroids_data_matrix[centroid] = sum(numpy.array([[data_matrix[gene_ids.index(gene)]] for gene in all_centroids[centroid]]))
			if len(centroids_data_matrix[centroid]) == 1:
				centroids_data_matrix[centroid] = [i for i in centroids_data_matrix[centroid][0]]
			else:
				centroids_data_matrix[centroid] = list(centroids_data_matrix[centroid][0])
		except:
			pdb.set_trace()
	foo = open('tmp/gene_centroids_table.txt','w')	
	foo.writelines(metadata[:-1])
	try:
		foo.writelines([str.join('\t', ['Centroids', metadata[-1].split('\t')[2:]])+'\n'])
	except:
		pdb.set_trace()
	for centroid in centroids_data_matrix:
		foo.writelines([str.join('\t', [centroid, [str(i) for i in centroids_data_matrix[centroid]]])+'\n'])
	foo.close()
	return centroids_data_matrix 

def get_prevalence_abundance(centroids_data_matrix, metadata):
	#Niche specific Beta-prevalence?
	centroid_prev_abund = {}
	all_prevalence = []
	all_mean_abund = []
	
	for centroid in centroids_data_matrix:
		centroid_prev_abund[centroid] = {'abund':numpy.mean(numpy.array(centroids_data_matrix[centroid])), 'prev': sum(numpy.array(centroids_data_matrix[centroid])>0)}
		all_prevalence += [centroid_prev_abund[centroid]['prev']]
		all_mean_abund += [centroid_prev_abund[centroid]['abund']]
	
	foo = open('tmp/centroid_prev_abund.txt','w')
	foo.writelines(['Centroids\tAbundance\tPrevalence\n'])
	for centroid in centroid_prev_abund:
		foo.writelines([str.join('\t', [centroid, str(centroid_prev_abund[centroid]['abund']), str(centroid_prev_abund[centroid]['prev'])])+'\n'])
	foo.close()

	return [centroid_prev_abund, all_prevalence, all_mean_abund]

#def get_important_centroids(centroid_prev_abund, all_prevalence, all_mean_abund):


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--input_table', help='Gene abundance table with metadata')
	parser.add_argument('-f','--fasta_folder', help='Folder containing fasta files')
	parser.add_argument('-u','--uclust_folder', nargs = '?' , help='Path for UCLUST program')
	parser.add_argument('-t','--type', help='Type of analysis Choose: [gene_table, reads_assemblies]')
	parser.add_argument('-o','--output_folder', help='Folder containing results')
	args = parser.parse_args()

	[metadata, uniref_gis, gis_unannotated, gene_ids, data_matrix] = read_gene_table(args.input_table)
	all_centroids = get_centroids(gis_unannotated, args.fasta_folder, metadata, args.uclust_folder, uniref_gis)
	centroids_data_matrix = get_centroids_table(gene_ids, all_centroids, data_matrix, metadata)
	[centroid_prev_abund, all_prevalence, all_mean_abund] = get_prevalence_abundance(centroids_data_matrix, metadata)

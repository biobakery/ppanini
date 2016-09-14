import os
import sys
import numpy
import re
import argparse
import logging
import pdb
from ppanini import utilities


logger = logging.getLogger(__name__)

def parse_annotation_table(annotations_file, fasta_sequences, thld_ref):
	'''Parses annotations result from RAPSEARCH/DIAMOND to ensure they pass threshold

	Input: annotations_file = path_to_RAPSEARCH/DIAMOND_results_file
		   fasta_sequences = {sequence_header: sequence}
		   thld_ref = int, threshold for sequence acceptance
	
	Output: [sample_annotations, search50_seqs]
	where, sample_annotations = {sequence: UniRef annotation}
		   diamond50_seqs = {sequence_header:  sequence} for all sequences that need to be redone'''

	logger.debug('parse_annotation_table '+annotations_file)

	foo = open(annotations_file)
	sample_annotations = {}
	for line in foo:
		if not line.startswith('#'):
			#FILE_FORMAT: Query_Seq<\t>Target_Seq<\t>Identity<\t>Aligned length
			split_line = [re.sub('[\t\r\n]', '', i) for i in line.split('\t')]

			query_length = float(len(fasta_sequences[split_line[0]]))
			aln_len = float(split_line[3])
			identity = float(split_line[2])
			thld = aln_len*identity/query_length
			
			if thld >= thld_ref:
				sample_annotations[split_line[0].strip()] = split_line[1].strip()
	
	search50_seqs = {}
	for seq in fasta_sequences:
		if not seq in sample_annotations:
			search50_seqs[seq] = fasta_sequences[seq]

	return [sample_annotations, search50_seqs]

def run_diamond(query_file, db, out_fname, all_paths, nprocesses):
	'''Runs DIAMOND on query_file to produce results in out_fname
	Input: query_file = path to query_fasta_file
		   all_paths = {'uniref90': path_to_uniref90_index, 
						'uniref50': path_to_uniref50_index, 
						'umap90_50': path_to_uniref90_uniref50_mapping}
		   out_fname = path of output_file to put the results in
		   nprocesses = Number of processes
		   db = DIAMOND preprocessed database'''

	logger.debug('run_diamond '+query_file)
	os.system(all_paths['diamond']+' blastp -q ' + query_file + ' \
													-d ' + db + ' \
													-k 1 \
													-a ' + out_fname + ' \
													-p ' + str(nprocesses)) #check if command is correct.

	os.system(all_paths['diamond']+' view -a ' + out_fname + '.daa \
												  -o ' + out_fname + '.m8 \
												  -p ' + str(nprocesses))

def run_rapsearch(query_file, db, out_fname, all_paths, nprocesses):
	'''Runs RAPSEARCH2 on query_file to produce results in out_fname

	Input: query_file = path to query_fasta_file
		   all_paths = {'uniref90': path_to_uniref90_index, 
						'uniref50': path_to_uniref50_index, 
						'umap90_50': path_to_uniref90_uniref50_mapping}
		   out_fname = path of output_file to put the results in
		   nprocesses = Number of processes
		   db = RAPSEARCH2 preprocessed database'''

	logger.debug('run_rapsearch '+query_file)
	os.system(all_paths['rapsearch']+' -q ' + query_file + ' \
												 -d ' + db + ' \
												 -o ' + out_fname + '.m8 \
												 -u 2 \
												 -b 0 \
												 -v 1 \
												 -z' + str(nprocesses))	

def run_uclust(usearch_folder, allgenes_file_path, gene_centroids_file_path, gene_centroid_clusters_file_path, perc_id, nprocesses):
	'''Runs USEARCH UCLUST on query_file to produce results in out_fname

	Input: usearch_folder = path to folder containing USEARCH
		   allgenes_file_path = path to input fasta file
		   gene_centroids_file_path = path to place all the centroids produced in
		   gene_centroid_clusters_file_path = path to file containing clustering results in
		   id = %ID to cluster sequences at
		   nprocesses= number of threads'''

	logger.debug('run_uclust '+allgenes_file_path)

	os.system(usearch_folder+ ' -cluster_fast ' + allgenes_file_path +' \
								 -id '+str(perc_id)+' \
								 -centroids '+ gene_centroids_file_path + ' \
								 -uc ' + gene_centroid_clusters_file_path+ '\
								 -threads '+str(nprocesses))

def run_vclust(vsearch_folder, allgenes_file_path, gene_centroids_file_path, gene_centroid_clusters_file_path, perc_id, nprocesses):
	'''Runs USEARCH UCLUST on query_file to produce results in out_fname

		Input: usearch_folder = path to folder containing USEARCH
				   allgenes_file_path = path to input fasta file
				   gene_centroids_file_path = path to place all the centroids produced in
				   gene_centroid_clusters_file_path = path to file containing clustering results in
				   id = %ID to cluster sequences at
				   nprocesses= number of threads'''

	logger.debug('run_vclust '+allgenes_file_path) 


	os.system(vsearch_folder+ ' --cluster_fast ' + allgenes_file_path +' \
								 --id '+str(perc_id)+' \
								 --centroids '+ gene_centroids_file_path + ' \
								 --uc ' + gene_centroid_clusters_file_path+ '\
								 --threads '+str(nprocesses))

def get_clusters_dict(gene_centroid_clusters_file_path):
	'''Return dict containing clusters
	Input: filepath to centroid file
	Output: centroid_gis (dict) {gene_centroid: [List of genes], }'''

	logger.debug('get_clusters_dict '+gene_centroid_clusters_file_path)

	cluster_txt = open(gene_centroid_clusters_file_path)
	centroid_gis = {}
	for line in cluster_txt:
		if line.startswith('H'):
			split_i = [re.sub('[\r\t\n]', '', i) for i in line.split('\t')[-2:]]
			try:
				centroid_gis[split_i[1]] += [split_i[0]]
			except KeyError:
				centroid_gis[split_i[1]] = [split_i[0], split_i[1]]
	return centroid_gis


def get_annotations_dict(centroid_annotations, centroid_gis):
	'''Returns annotations for all the genes and centroids
	Input: centroid_annotations (dict) {centroid: annotation, ...}
		   centroid_gis (dict) {centroid: [List of genes], ...}
	Output: annotation_dict (dict) {gene: annotation, ...}'''
	
	logger.debug('get_annotations_dict')

	annotation_dict = {}
	for centroid in centroid_annotations:
		try:
			gis = centroid_gis[centroid]
			for gi in gis:
				annotation_dict[gi] = centroid_annotations[centroid]
		except:
			annotation_dict[centroid] = centroid_annotations[centroid]
	return annotation_dict

	
if __name__ == '__main__':
	pass

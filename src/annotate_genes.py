import os
import sys
import numpy
import re
import argparse

from . import utilities

def parse_annotation_table(annotations_file, fasta_sequences, thld_ref):
	'''Parses annotations result from RAPSEARCH/DIAMOND to ensure they pass threshold

	Input: annotations_file = path_to_RAPSEARCH/DIAMOND_results_file
		   fasta_sequences = {sequence_header: sequence}
		   thld_ref = int, threshold for sequence acceptance
	
	Output: [sample_annotations, search50_seqs]
	where, sample_annotations = {sequence: UniRef annotation}
		   diamond50_seqs = {sequence_header:  sequence} for all sequences that need to be redone'''

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

def run_diamond(query_file, db, out_fname, nprocesses):
	'''Runs DIAMOND on query_file to produce results in out_fname
	Input: query_file = path to query_fasta_file
		   all_paths = {'uniref90': path_to_uniref90_index, 
                                            'uniref50': path_to_uniref50_index, 
                                            'umap90_50': path_to_uniref90_uniref50_mapping}
           out_fname = path of output_file to put the results in
           nprocesses = Number of processes
           db = DIAMOND preprocessed database'''

	os.system(all_paths['diamond']+'/diamond blastp -q ' + query_file + ' \
													-d ' + db + ' \
													-k 1 \
													-o ' + out_fname + ' \
													-p ' + str(nprocesses)) #check if command is correct.

def run_rapsearch(query_file, db, out_fname, nprocesses):
	'''Runs RAPSEARCH2 on query_file to produce results in out_fname

	Input: query_file = path to query_fasta_file
		   all_paths = {'uniref90': path_to_uniref90_index, 
                                            'uniref50': path_to_uniref50_index, 
                                            'umap90_50': path_to_uniref90_uniref50_mapping}
           out_fname = path of output_file to put the results in
           nprocesses = Number of processes
           db = RAPSEARCH2 preprocessed database'''
    
	os.system(all_paths['rapsearch']+'/rapsearch -q ' + query_file + ' \
												 -d ' + db + ' \
												 -o ' + out_fname + ' \
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
		   id = %ID to cluster sequences at'''

	if not usearch_folder:
		usearch_folder = 'usearch'
	else:
		usearch_folder = usearch_folder

	os.system(usearch_folder + ' -cluster_fast ' + allgenes_file_path +' \
								  		-id '+str(perc_id)+' \
								        -centroids '+ gene_centroids_file_path + ' \
								        -uc ' + gene_centroid_clusters_file_path+ '\
								        -threads '+str(nprocesses))

def run_vclust(usearch_folder, allgenes_file_path, gene_centroids_file_path, gene_centroid_clusters_file_path, perc_id, nprocesses):
        '''Runs USEARCH UCLUST on query_file to produce results in out_fname

        Input: usearch_folder = path to folder containing USEARCH
                   allgenes_file_path = path to input fasta file
                   gene_centroids_file_path = path to place all the centroids produced in
                   gene_centroid_clusters_file_path = path to file containing clustering results in
                   id = %ID to cluster sequences at'''

        if not usearch_folder:
                usearch_folder = 'vsearch'
        else:
                usearch_folder = usearch_folder

        os.system(usearch_folder + ' --cluster_fast ' + allgenes_file_path +' \
                                     --id '+str(perc_id)+' \
                                     --centroids '+ gene_centroids_file_path + ' \
                                     --uc ' + gene_centroid_clusters_file_path+ '\
                                     --threads '+str(nprocesses))

def get_clusters_dict(gene_centroid_clusters_file_path):
	'''Return dict containing clusters
	Input: filepath to centroid file
	Output: centroid_gis (dict) {gene_centroid: [List of genes], }'''

	cluster_txt = open(gene_centroid_clusters_file_path)
	centroid_gis = {}
	for line in cluster_txt:
		if line.startswith('H'):
			split_i = [re.sub('[\r\t\n]', '', i) for i in line.split('\t')]
			try:
				centroid_gis[split_i[-1]] += [split_i[-2]]
			except KeyError:
				centroid_gis[split_i[-1]] = [split_i[-2], split_i[-1]]
	return centroid_gis


def get_annotations_dict(centroid_annotations, centroid_gis):
	'''Returns annotations for all the genes and centroids
	Input: centroid_annotations (dict) {centroid: annotation, ...}
		   centroid_gis (dict) {centroid: [List of genes], ...}
	Output: annotation_dict (dict) {gene: annotation, ...}'''

	annotation_dict = {}
	for centroid in centroid_annotations:
		gis = centroid_gis[centroid]
		for gi in gis:
			annotation_dict[gi] = centroid_annotations[centroid]
	return annotation_dict

	
if __name__ == '__main__':
	pass
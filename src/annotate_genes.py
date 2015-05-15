import os
import sys
import numpy
import re
import argparse

from src import utilities

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
			#FILE_FORMAT: Seq_Name<\t>Seq_Length<\t>Mapped_Reads<\t>Unmapped_Reads
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
	cluster_txt = open(gene_centroid_clusters_file_path)#os.popen('grep -w H ' + gene_centroid_clusters_file_path)
	centroid_gis = {}
	for line in cluster_txt:#.xreadlines():
		if line.startswith('H'):
			split_i = [re.sub('[\r\t\n]', '', i) for i in line.split('\t')]
			try:
				centroid_gis[split_i[-1]] += [split_i[-2]]
			except KeyError:
				centroid_gis[split_i[-1]] = [split_i[-2], split_i[-1]]
	return centroid_gis
	
def get_genes_samples(mapper):
	genes_samples = {} #Dicitonary of genes and their samples

	for sample in mapper:
		seqs_i = utilities.read_fasta(mapper[sample]['FAAS'])
		seqs_i = seqs_i.keys()

		for gene in seqs_i:
			if not gene in genes_samples:
				genes_samples[gene] = sample
			else:
				raise Exception('Non-unique gene identifiers across samples; please rename gene ids')
	return genes_samples

def get_annotations_dict(centroid_annotations, centroid_gis):
	annotations_dict = {}
	for centroid in centroid_annotations:
		gis = centroid_gis[centroid]
		for gi in gis:
			annotations_dict[gi] = centroid_annotations[centroid]

	return annotations_dict
# def get_annotations_dict(mapper, genes_samples, centroid_annotations, centroid_gis):
# 	annotations_dict = {}
# 	for sample in mapper:
# 		annotations_dict[sample] = {}

# 	for centroid in centroid_annotations:
# 		annot_i = centroid_annotations[centroid]
		
# 		if not centroid in centroid_gis:
# 			centroid_gis[centroid] = [centroid] #If no cluster with centroid, then one member cluster

# 		genes_i = centroid_gis[centroid]
# 		for gene in genes_i:
# 			if not gene in annotations_dict[genes_samples[gene]]:
# 				annotations_dict[genes_samples[gene]][gene] = annot_i
# 			else:
# 				raise Exception('Multiple annotation assignment for the same gene')

# 	return annotations_dict
######################
##WITHPANGENOME
######################
def generate_annotation(mapper, all_paths, nprocesses, basename):
	'''Returns an annotation dictionary {sample: {gene:annotation}, ...}
	'''
	allgenes_file_path = 'tmp/'+basename+'_preppanini_allsequences.fasta'
	gene_centroids_file_path = 'tmp/'+basename+'_preppanini_centroids.fasta'
	gene_centroid_clusters_file_path = 'tmp/'+basename+'_preppanini_centroid_clusters.uc'
	out_u90_fname = 'tmp/'+basename+'_preppanini_centroids_u90'
	out_u50_fname = 'tmp/'+basename+'_preppanini_centroids_u50'
	u50_gene_input = 'tmp/'+basename+'_preppanini_centroids_u50input.fasta'

	os.system('cat '+str.join(' ', [mapper[i]['FAAS'] for i in mapper[i]]) + ' > ' + all_fastas_seqs)
	run_uclust(all_paths['usearch'], allgenes_file_path, gene_centroids_file_path, gene_centroid_clusters_file_path, 0.9, nprocesses)

	run_diamond(gene_centroids_file_path, all_paths, out_fname, nprocesses, 'uniref90')

	centroid_sequences = utilities.read_fasta(gene_centroids_file_path)
	[centroid_annotations90, diamond50_seqs] = parse_annotation_table(u90_out_fname+'.m8', centroid_sequences, 90.0)
	
	centroid_annotations = centroid_annotations90

	if not diamond50_seqs == {}:
		utilities.write_fasta(diamond50_seqs, u50_gene_input, True)
		run_diamond(u50_gene_input, all_paths, u50_out_fname, nprocesses, 'uniref50')
		[centroid_annotations50, diamondukn_seqs] = parse_annotation_table(u50_out_fname+'.m8', diamond50_seqs, 50.0)
		
		for gid in centroid_annotations50:
			if gid not in centroid_annotations:
				centroid_annotations[gid] = sample_annotations50[gid]
			else:
				raise Exception('GeneID for UniRef50 exists in UniRef90;\nSearching twice for the same gene;\nError in parse_annotation_table\n')

	genes_samples = get_genes_samples(mapper)
	centroid_gis = get_clusters_dict(gene_centroid_clusters_file_path)

	return get_annotations_dict(mapper, genes_samples, centroid_annotations, centroid_gis)

	
if __name__ == '__main__':
	pass
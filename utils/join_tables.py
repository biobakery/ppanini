import os
import sys
import re
import pdb
import time

def read_table(gene_table_fname):
	'''Returns the different elements from the gene table

	Input: gene_table_fname = Filename of the gene_table

	Output: metadata = [metadata strings]; Rows with # as first character in table
			uniref_gis = {UniRef_XYZ: [list of gene ids]}
			gis_unannotated = {sample_id: [gene ids]}
			gene_ids = [List of all gene ids]
			data_matrix = The abundance table [[0, 0,.],[0, 0,.],...]'''

	gene_table = open(gene_table_fname) 
	gtab_lines = gene_table.readlines()
	samples = []
	genes, metadata = {}, {}

	
	for line in gtab_lines:
		#Metadata containing Sample names, Niche specifics and optionally fasta file locations
	#	samples = []
		if line.startswith('#'):
			split_i = line.split('\t')
			metadata[split_i[0]]= [re.sub('[\r\t\n]','',i) for i in split_i[1:]]
			samples = metadata[split_i[0]]	
		if not line.startswith('#'):
			split_i = line.split('\t')
			genes[split_i[0]] = [re.sub('[\r\t\n]','',i) for i in split_i[1:]]
	return [metadata, samples, genes]

def join(orig_table, add_table):
	[metadata, samples, genes] = orig_table
	[metadata_i, samples_i, genes_i] = add_table

	metadata_final = []
	for key in metadata:
		metadata_final += ['\t'.join([key]+metadata[key]+metadata_i[key])]

	samples_final = samples + samples_i
	zeros_tba_o = len(samples_i)
	zeros_tba_i = len(samples)
	genes_final = {}
	
	for gene in genes:
		genes_final[gene] = genes[gene] + ["0" for i in range(zeros_tba_o)]

	for gene in genes_i:
		tmp =["0" for i in range(zeros_tba_i)] + genes_i[gene]
		genes_final[gene] = tmp
	#pdb.set_trace()
	return [metadata_final, samples_final, genes_final] 

def print_table(metadata_final, genes_final):
	for line in metadata_final:
		print line.strip()
	for gene in genes_final:
		row = [gene]+genes_final[gene]
		#print len(row)
		print '\t'.join([gene]+genes_final[gene]).strip()

if __name__ == '__main__':
	list_tables = []
	t1 = time.time()
	for i in range(1, len(sys.argv)):
		#print sys.argv[i]
		list_tables +=[read_table(sys.argv[i])]
	t2 = time.time()
	print 'Read tables '+str((t2-t1)/60)+' mins elapsed'
	[metadata_final, samples_final, genes_final] = list_tables[0]

	for i in range(1, len(list_tables)):
		[metadata_final, samples_final, genes_final] = join([metadata_final, samples_final, genes_final], list_tables[i])
	t3 = time.time()
	print 'Joined tables, now printing '+str((t3-t2)/60)+' mins elapsed'
	print_table(metadata_final, genes_final)








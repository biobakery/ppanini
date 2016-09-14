import os
import sys
import re
import time

def read_table(gene_table_fname):
	'''Returns the different elements from the gene table

	Input: gene_table_fname = Filename of the gene_table

	Output: metadata = {'#KEY': [row split ]} 
			genes = [] order of genes as they appear in the table
			data_matrix = [['0','0','1', ...],...] matrix of the gene abundance in strings'''

	gene_table = open(gene_table_fname) 
	genes = []
	metadata = {}
	data_matrix = []
	
	for line in gene_table:
		if line.startswith('#'):
			split_i = line.split('\t')
			metadata[split_i[0]]= [re.sub('[\r\t\n]','',i) for i in split_i[1:]]
		else:
			split_i = line.split('\t')
			genes += [split_i[0]] 
			data_matrix += [[re.sub('[\r\t\n]','',i) for i in split_i[1:]]]
	return [metadata, genes, data_matrix]
def main():
	cmd_h = ['-h', '--help']
	if sys.argv[1] in cmd_h:
		print 'usage: python utils/join_tables.py <table1> <table2> ... > merged_table.txt'
		sys.exit(0)

	original_table = sys.argv[1]
	[metadata_o, genes_o, data_matrix_o] = read_table(original_table)

	for i in range(2, len(sys.argv)):
		new_table = sys.argv[i]
		print 'Adding '+new_table+' to original table'

		[metadata, genes, data_matrix] = read_table(new_table)
		
		no_added_o = [len(metadata['#SAMPLES']), len(genes_o)]
		o_row = ['0' for j in range(len(metadata['#SAMPLES']))]

		no_added_i = [len(metadata_o['#SAMPLES']), len(genes)]
		i_row = ['0' for j in range(len(metadata_o['#SAMPLES']))]
		
		for j, val in enumerate(data_matrix_o): #Adding cols to old table
			data_matrix_o[j] = val+o_row

		for j, val in enumerate(data_matrix): #Adding cols to new table for old table
			data_matrix[j] = i_row+val
		
		data_matrix_o = data_matrix_o+data_matrix
		genes_o = genes_o + genes
		
		for key in metadata: #Assumes all keys in metadata exist in metadata_o
			metadata_o[key] += metadata[key]

	for key in metadata_o:
		print '\t'.join([key]+metadata_o[key])
	for i, val in enumerate(genes_o):
		print '\t'.join([val]+data_matrix_o[i])

if __name__ == '__main__':
	main()
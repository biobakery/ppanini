import os
import sys
import re

def write_dict(dictX, keys, keys, gene_annotations_file):
	'''Writes dictionary of genes and their annotations into text file
	Input: dictX = {geneID: annotation}
		   gene_annotations_file = path_to_output_gene_annotations_table'''

	with open(gene_annotations_file, 'w') as foo:
		if not keys:
			foo.writelines('#ID\tCOLUMN2\n')
			for i in dictX:
				foo.writelines(['\t'.join([str(i), str(dictX[i])])+'\n'])
		else:
			foo.writelines('#ID'+'\t'.join(keys)+'\n')
			for i in dictX:
				foo.writelines(['\t'.join([str(i)] + [str(dictX[i][j]) for j in keys])+'\n'])
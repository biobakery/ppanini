
import os
import sys
import pdb
import re
import numpy
import logging
import argparse
import subprocess
import multiprocessing

logger = logging.getLogger(__name__)

def generate_gene_table(abundance_dict, annotations_dict, niche_flag, mapper, output_table, samples):
	'''Input: abundance_dict: {sample:{gene:abundance,...},...}
			  annotation_dict: {gene: annotation},
			  niche_flag: True if NICHE exists
			  mapper'''
	logger.debug('generate_gene_table '+output_table)
#	samples = abundance_dict.keys()
	# fasta_row = [mapper[i]['FAAS'] for i in samples]

	with open(output_table, 'w') as foo:
		if niche_flag:
			niche_row = [mapper[i]['NICHE'] for i in samples]
			foo.writelines([str.join('\t', ['#NICHE']+niche_row)+'\n']) #header
		# foo.writelines([str.join('\t', ['#FAAS']+fasta_row)+'\n']) #header
		foo.writelines([str.join('\t', ['#SAMPLES']+samples)+'\n']) #header
		
		#for i, sample in enumerate(samples):
		for gene in abundance_dict:
			#abund_x_i = abundance_dict[sample][gene]
			#data_row = numpy.zeros(len(samples))
			#data_row[i] = abund_x_i
			str_data_row = [str(ele) for ele in abundance_dict[gene]]
			if gene in annotations_dict: 
				annot_x_i = annotations_dict[gene]
				if annot_x_i.startswith('UniRef90'):
					umap_i = 'UniRef50_unknown'
					annot_x = annot_x_i + '|' + umap_i
				else:
					annot_x = 'UniRef90_unknown|' + annot_x_i 
			else:
				annot_x = 'UniRef90_unknown|UniRef50_unknown'
			foo.writelines([str.join('\t', [gene+'|'+annot_x]+str_data_row)+'\n'])

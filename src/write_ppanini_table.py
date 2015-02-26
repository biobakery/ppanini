
import os
import sys
import pdb
import re
import argparse
import numpy
import subprocess
import multiprocessing



def generate_gene_table(abundance_dict, annotations_dict, all_paths, niche_flag, mapper, output_table):

	#annotations_dict = mapper_with_annotations_dict['annotations_dict']	#sample:{gene:annotation}
	#abundance_dict = mapper_with_annotations_dict['abundance_dict'] #sample:{gene:abundance}
	umap90_50 = read_id_mapping(all_paths['uniref_map'])
	
	samples = abundance_dict.keys()
	fasta_row = [mapper[i]['FASTAS'] for i in samples]

	if niche_flag:
		niche_row = [mapper[i]['NICHE'] for i in samples]

	with open(output_table, 'w') as foo:
		if niche_flag:
			foo.writelines([str.join('\t', ['#NICHE']+niche_row)+'\n']) #header

		foo.writelines([str.join('\t', ['#FASTAS']+fasta_row)+'\n']) #header
		foo.writelines([str.join('\t', ['#GENES']+samples)+'\n']) #header
		
		for i, sample in enumerate(samples):
			for gene in abundance_dict[sample]:
		
				abund_x_i = abundance_dict[sample][gene]
				data_row = numpy.zeros(len(samples))
				data_row[i] = abund_x_i
				str_data_row = [str(ele) for ele in data_row]

				if gene in annotations_dict[sample][gene]: 
					annot_x_i = annotations_dict[sample][gene]
					if annot_x.startswith('UniRef90'):
						annot_x = annot_x_i + '|' + umap90_50[annot_x_i]
					else:
						annot_x = 'UniRef90_unknown|' + annot_x_i 
				else:
					annot_x = 'UniRef90_unknown|UniRef50_unknown'
				
				foo.writelines([str.join('\t', [gene+'|'+annot_x]+str_data_row)+'\n'])
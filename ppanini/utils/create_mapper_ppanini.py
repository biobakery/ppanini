import pdb
import os
import sys
import matplotlib
import re
import numpy
import time
import numpy
import argparse
import write_dictionary
from . import utilities
from matplotlib import pyplot

'''To write the mapper for plot2vals'''

numpy.seterr(divide='ignore', invalid='ignore')

def read_go_table(map_go_fname):
	map = open(map_go_fname,'r')
	mapper = {}
	for line in map:
		if not line.startswith('#'):
			split_line = line.split('\t')
			for i in split_line[1:]:
				mapper[re.sub('[\r\t\n]','',i)] = split_line[0]
	return mapper

def get_go_mapping(genes, mapper):
	go_genes = []
	for gene in genes:
		if gene in mapper:
			go_genes +=[gene]
	print len(go_genes)
	return go_genes

def read_genes(input_table):
	genes = []
	with open(input_table) as foo:
		for line in foo:
			if not line.startswith('#'):
				genes += [re.sub('[\r\t\n]','',line.split('\t')[0]).strip()]
	print len(genes)
	return genes

def get_final_mapper(genes, go_mapper, labels):
	go_genes = get_go_mapping(genes, go_mapper)
	final_map = {}
	for gene in genes:
		final_map[gene] = {}
		if gene.startswith('UniRef90'):
			final_map[gene]= {'color': labels['color_uniref90'],\
							  'zorder': labels['zorder_uniref90']}
		else:
			final_map[gene]= {'color': labels['color_unknown'],\
							  'zorder': labels['zorder_unknown']}
	for gene in go_genes:
		final_map[gene]= {'color': labels['color_go'],\
							  'zorder': labels['zorder_go']}
	print(len(final_map))
	return final_map

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--input_table', help='Gene abundance table with metadata', required=True)
	parser.add_argument('--go-mapper', dest = 'go_mapper', help='', required=True)
	parser.add_argument('--basename', default=False, help='')
	parser.add_argument('--color-go', dest = 'color_go', default='blue', help='')
	parser.add_argument('--color-uniref90', dest ='color_uniref90', default='red', help='')
	parser.add_argument('--color-unknown', dest ='color_unknown', default='yellow', help='')
	parser.add_argument('--zorder-go', dest ='zorder_go', default=2, help='')
	parser.add_argument('--zorder-uniref90', dest='zorder_uniref90', default=1, help='')
	parser.add_argument('--zorder-unknown', dest = 'zorder_unknown', default=3, help='')
	parser.add_argument('-o','--output-filename', dest ='output_filename',default=False, help='filename for figure')

	args = parser.parse_args()
	go_mapper = read_go_table(args.go_mapper)
	genes = read_genes(args.input_table)
	
	basename = args.input_table.split('/')[-1].split('.')[0]
	
	if args.basename:
		basename = basename

	labels = {'color_go': args.color_go,\
	'color_uniref90': args.color_uniref90,\
	'color_unknown': args.color_unknown,\
	'zorder_go': int(args.zorder_go),\
	'zorder_uniref90': int(args.zorder_uniref90),\
	'zorder_unknown': int(args.zorder_unknown)}

	out = basename+'_mapper.txt'
	if args.output_filename:
		out = args.output_filename

	final_map = get_final_mapper(genes, go_mapper, labels)
	write_dictionary.write_dict(final_map, ['color', 'zorder'], out)



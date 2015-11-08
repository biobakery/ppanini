import pdb
import os
import sys
import matplotlib
import re
import numpy
import time
import numpy
import logging
import argparse

from src import utilities
from matplotlib import pyplot

numpy.seterr(divide='ignore', invalid='ignore')

basename = ''
logger = logging.getLogger(__name__)

def read_map(table_name):
	table_obj = open(table_name)
	keys_list = []
	keys ={}
	genes_cmap = {}
	for line in table_obj:
		if line.startswith('#'):
			keys_list = [str.lower(re.sub('[#\r\t\n]','',i).strip()) for i in line.split('\t')] #Assuming #Centroids\tCol1\tCol2
			try:
				keys['color'] = keys_list.index('color')
			except:
				pass
			try:
				keys['zorder'] = keys_list.index('zorder')
			except:
				pass
			if keys == {}:
				raise Exception('Specify header for mapper file; mark it with # at the beginning of the line; and specify zorder and/or color')
		else:
			split_i = [re.sub('[\r\t\n]','',i).strip() for i in line.split('\t')]
			genes_cmap[split_i[0]] = {}
			for i in keys:
				genes_cmap[split_i[0]][i]= split_i[keys[i]]	
	return genes_cmap


def parse_table(table_name, cols):
	table_obj = open(table_name)
	genes_prab = {}
	keys = []
	for line in table_obj:
		split_i = [re.sub('[\r\t\n]','',i).strip() for i in line.split('\t')]
		if line.startswith('#'):
			keys = [split_i[cols[0]], split_i[1]]
		else:
			genes_prab[split_i[0]] = tuple([float(split_i[cols[0]]), float(split_i[cols[1]])]) #(prevalence, abundance) always <--

	return [genes_prab, keys]

def populate_cmap(genes_cmap, genes_prab):
	for gene in genes_prab:
		if not gene in genes_cmap:
			genes_cmap[gene]={}
			genes_cmap[gene]['zorder'] = 1
			genes_cmap[gene]['color'] = 'gray'
		else:
			if not 'zorder' in genes_cmap[gene]:
				genes_cmap[gene]['zorder'] = 1
			if not 'color' in genes_cmap[gene]:
				genes_cmap[gene]['color'] = 'gray'
	return genes_cmap

def split_cmap(genes_cmap, genes_prab):
	zorder_dict = {}
	color_dict = {}
	prevalence_dict = {}
	abundance_dict = {}
	ind = 0 
	zorder_dict[ind] = {}
	color_dict[ind] = {}
	prevalence_dict[ind] = {}
	abundance_dict[ind] = {}
	for gene in genes_cmap:
		z = genes_cmap[gene]['zorder']
		try:
			zorder_dict[z] +=[gene]
		except:
			zorder_dict[z] =[gene]
			color_dict[z] = {}
			prevalence_dict[z] = {}
			abundance_dict[z] = {}
		color_dict[z][gene] = genes_cmap[gene]['color']
		try:
			prevalence_dict[z][gene] = genes_prab[gene][0]
		except:
			pdb.set_trace()
		abundance_dict[z][gene] = genes_prab[gene][1]
	return [zorder_dict, color_dict, prevalence_dict, abundance_dict]


def plot_prev_abund(zorder_dict, color_dict, prevalence_dict, abundance_dict, labels, margins):
	pyplot.figure()
	pyplot.xlabel(labels['xlabel'])#'Alpha Prevalence_'+keys[i])
	pyplot.ylabel(labels['ylabel'])
	pyplot.title(labels['title']) #'Mean abundance')
	
	for i in zorder_dict:
		print i
		pyplot.scatter(numpy.log(prevalence_dict[i].values(), dtype='float64'), \
				   numpy.log(abundance_dict[i].values(), dtype='float64'), \
				   c=color_dict[i].values(), \
				   alpha=labels['alpha'], \
				   linewidths=0.0, \
				   zorder=i, \
				   marker='o')
	if margins:
		[prev_tshld, abund_tshld] = [float(i) for i in margins.split(',')]
		pyplot.axvline(y=numpy.log(float(margins[0])), alpha=0.5, color='gray')
		pyplot.axhline(x=numpy.log(float(margins[1])), alpha=0.5, color='gray')

	pyplot.savefig(labels['filename'])
	pyplot.savefig(labels['filename']+'.png')

def populate_vars(labels, args):
	if args.xlabel:
		labels['xlabel'] = args.xlabel
	if args.xlabel:
		labels['ylabel'] = args.ylabel
	if args.xlabel:
		labels['title'] = args.title
	if args.xlabel:
		labels['alpha'] = float(args.alpha)
	if args.xlabel:
		labels['filename'] = args.output_filename+'.'+args.format

	return labels
		
if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--input_table', help='Genes with prevalence and abundance values', required=True)
	parser.add_argument('-m', '--map_file', default=False, help='Mapping file of genes to color and zorder values', required=True)
	parser.add_argument('--margins', default=False, help='Add threshold margins: prevalence_threshold,abundance_threshold: e.g. 0.001,0.001')
	parser.add_argument('--alpha', default=False, help='Transparency')
	parser.add_argument('--cols', default='1,2', help='Columns to plot: x,y')
	parser.add_argument('--xlabel', default=False, help='Custom xlabel')
	parser.add_argument('--ylabel', default=False, help='Custom ylabel')
	parser.add_argument('--title', default=False, help='Custom title')
	parser.add_argument('--basename', default=False, help='Basename')
	parser.add_argument('--hexplot', default=False, help='Plot hexplots instead of')
	parser.add_argument('-o','--output_filename', default=False, help='filename for figure')
	parser.add_argument('--format', default='pdf', help='Format of file')

	args = parser.parse_args()
	basename = args.input_table.split('/')[-1].split('.')[0]
	
	if args.basename:
		basename = basename
	cols = [int(i) for i in args.cols.split(',')]
	[genes_prab, ids] = parse_table(args.input_table, cols)

	if args.map_file:
		genes_cmap = read_map(args.map_file)
	else:
		genes_cmap = {}

	genes_cmap = populate_cmap(genes_cmap, genes_prab)

	labels = {'xlabel': ids[0],\
			  'ylabel': ids[1],\
			  'title': 'Scatterplot',\
			  'alpha': 0.1,\
			  'filename': basename+'.'+args.format}

	labels = populate_vars(labels, args)

	[zorder_dict, color_dict, prevalence_dict, abundance_dict] = split_cmap(genes_cmap, genes_prab)
	plot_prev_abund(zorder_dict, color_dict, prevalence_dict, abundance_dict, labels, args.margins)
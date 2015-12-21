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

from . import utilities
from matplotlib import pyplot

'''Tmp file to parse results'''

numpy.seterr(divide='ignore', invalid='ignore')

basename = ''
logger = logging.getLogger(__name__)

def read_map(table_name):
	table_obj = open(table_name)
	keys = []
	genes_cmap = {}
	for line in table_obj:
		if line.startswith('#'):
			keys = [re.sub('[\r\t\n]','',i).strip() for i in line.split('\t')[1:]] #Assuming #Centroids\tCol1\tCol2
			if not keys:
				raise Exception('Specify header for mapper file; mark it with # at the beginning of the line; and specify zorder and/or color')
		else:
			split_i = [re.sub('[\r\t\n]','',i).strip() for i in line.split('\t')]
			genes_cmap[split_i[0]] = {}
			for i in keys:
				genes_cmap[split_i[0]][keys[i]]= float(split_i[i])
	return genes_cmap


def parse_table(table_name, beta_flag):
	table_obj = open(table_name)
	genes_prab = {}
	inds = []
	for line in table_obj:
		if line.startswith('#'):
			keys = [re.sub('[\r\t\n]','',i).strip() for i in line.split('\t')[1:]]
			try:
				if not beta_flag:
					alpha_i = [i for i, val in enumerate(keys) if 'alpha' in val][0] #need to add multiple alpha prevalence functionality
					inds+=[alpha_i]
				else:
					beta_i = [i for i, val in enumerate(keys) if 'beta' in val][0]
					inds+=[beta_i]
				
				abund_i = [i for i, val in enumerate(keys) if 'abund' in val][0]
				inds+=[abund_i]
			except:
				raise Exception('No Metadata provided: Please insert Metadata row and insert # at the beginning\n \
						The keywords for alpha_prevalence, beta_prevalence and mean_abundance are: alpha, beta, abund\n \
						Atleast one of alpha or beta prevalence **AND** mean abundance should be present. See demo file for example.')
		else:
			split_i = [re.sub('[\r\t\n]','',i).strip() for i in line.split('\t')]
			genes_prab[split_i[0]] = tuple([float(split_i[inds[0]]), float(split_i[inds[1]])]) #(prevalence, abundance) always <--

	return genes_prab

def populate_cmap(genes_cmap, genes_prab):
	for gene in genes_prab:
		if not gene in genes_cmap:
			genes_cmap[gene]['zorder'] = 1
			genes_cmap[gene]['color'] = 'gray'
		else:
			if not 'zorder' in genes_cmap[gene]:
				genes_cmap[gene]['zorder'] = 1
			if not 'color' in genes_cmap[gene]:
				genes_cmap[gene]['color'] = 'gray'
	return genes_cmap

def plot_prev_abund(genes_prab, genes_cmap, labels, margins):
	genes = []
	prev = []
	abund = []
	color = []
	zorder = []

	for gene in genes_prab:
		genes +=[gene]
		prev +=[genes_prab[gene][0]]
		abund +=[genes_prab[gene][1]]
		zorder +=[genes_cmap[gene]['zorder']]
		color +=[genes_cmap[gene]['color']]

	pyplot.figure()
	pyplot.xlabel(labels['xlabel'])#'Alpha Prevalence_'+keys[i])
	pyplot.ylabel(labels['ylabel'])
	pyplot.title(labels['title']) #'Mean abundance')
	
	pyplot.scatter(numpy.log(prev, dtype='float64'), \
				   numpy.log(abund, dtype='float64'), \
				   c=color, \
				   alpha=labels['alpha'], \
				   linewidths=0.0, \
				   zorder=zorder, \
				   marker='o')
	if margins:
		[prev_tshld, abund_tshld] = [ float(i) for i in margins.split(',')]
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
	parser.add_argument('--beta_flag', default=False, action="store_true",help='Beta_Prevalence against Abundance')
	parser.add_argument('--margins', default=False, help='Add threshold margins: prevalence_threshold,abundance_threshold: e.g. 0.001,0.001')
	parser.add_argument('--alpha', default=False, help='Transparency')
	parser.add_argument('--xlabel', default=False, help='Custom xlabel')
	parser.add_argument('--ylabel', default=False, help='Custom ylabel')
	parser.add_argument('--title', default=False, help='Custom title')
	parser.add_argument('-o','--output_filename', default=False, help='filename for figure')
	parser.add_argument('--format', default='pdf', help='Format of file')

	args = parser.parse_args()
	basename = args.input_table.split('/')[-1].split('.')[0]
	
	if args.basename:
		basename = basename
	
	[genes_prab, ids] = parse_table(args.input_table)
	
	if args.map_file:
		[genes_cmap, keys] = read_map(args.map_file)
	else:
		genes_cmap = {}

	genes_cmap = populate_cmap(genes_cmap, genes_prab)

	labels = {'xlabel': 'Prevalence',\
			  'ylabel': 'Mean Abundance',\
			  'title': 'Scatterplot: Prevalence and abundance',\
			  'alpha': 0.1,\
			  'filename': basename+'.'+args.format}
	labels = populate_vars(labels, args)

	plot_prev_abund(genes_prab, genes_cmap, labels, args.margins)



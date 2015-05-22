import os
import sys
import matplotlib
import re
import numpy
import pdb
import time
import numpy
import argparse
import ppanini
import src

#from ppanini import src
from matplotlib import pyplot
from src import utilities

'''Analysis of Genome Hits per gene showing how many genomes each gene is found in.'''

def parse_table(m8_filename, fasta_filename):
	'''Parse the BLAST results to give gene hits to genomes
	Input: filename of blast results
	Output: table = {gene: [List of genomes]}'''

	fasta_dict = utilities.read_fasta(fasta_filename)

	for seq in fasta_dict:
		fasta_dict[seq] = float(len(fasta_dict[seq]))
	table = {}
	foo = open(m8_filename)
	for line in foo:
		split_i = line.split('\t')
		try:
			threshold = float(split_i[2])*float(split_i[3])/fasta_dict[split_i[0]]
		except:
			raise Exception('Gene '+split_i[0]+' not found in fasta: '+fasta_filename)
		if threshold > 90.0:
			sp = split_i[1].split('|')[-3]
			if split_i[0] not in table:
				table[split_i[0]] = [sp]
			elif sp not in table[split_i[0]]:
				table[split_i[0]] += [sp]
	return table

def read_parsed(m8_filename):
	'''Read parsed table for {gene: genomes}
	Input: filename of blast results
	Output: table = {gene: [List of genomes]}'''
	table = {}
	foo = open(m8_filename)

	for line in foo:
		split_i = [i.strip() for i in line.split('\t')]
		try:
			table[split_i[0]] += [split_i[1]]
		except:
			table[split_i[0]] = [split_i[1]]
	return table

def plot_scatter(table, m8_filename, no_uniq_genomes):
	'''Plot Scatter plot for genome hits per gene'''
	labels = {'xlabel': 'Prioritized Centroids',\
			  'ylabel':'No. of Genomes (Log10)', \
			  'title':'Metagenome vs. Genome Prioritization',\
			  'filename': m8_filename+'_prioritizationScatter.pdf'}	

	all_genes = []
	for gene in table:
		all_genes +=[len(table[gene])/float(no_uniq_genomes)]

	all_genes.sort()
	pyplot.figure()
	pyplot.xlabel(labels['xlabel'])
	pyplot.ylabel(labels['ylabel'])
	pyplot.title(labels['title']) 
	pyplot.scatter(numpy.arange(len(all_genes)), \
				   numpy.log(all_genes), \
				   c='grey', \
				   alpha=0.1, \
				   linewidths=0.0, \
				   zorder=1, \
				   marker='o',\
				   label='All Centroids')
	pyplot.legend( loc=4, \
				   fontsize='x-small', \
				   framealpha=0.4, )	
	pyplot.savefig(labels['filename'])
	pyplot.savefig(labels['filename']+'.png')

def plot_hexbin(table, m8_filename):
	'''Plots HexBin plots for the genome hits per gene'''
	labels = {'xlabel': 'Prioritized Centroids',\
			  'ylabel':'No. of Genomes', \
			  'title':'Metagenome vs. Genome Prioritization',\
			  'filename': m8_filename+'_prioritizationHEXBIN.pdf'}
	pyplot.figure()
	pyplot.xlabel(labels['xlabel'])
	pyplot.ylabel(labels['ylabel'])
	pyplot.title(labels['title']) 

	all_genes = []
	for gene in table:
		all_genes +=[len(table[gene])/4189.0]
	all_genes.sort()
	pyplot.hold(True)
	image = pyplot.hexbin(numpy.arange(len(all_genes)), \
				   all_genes,\
				   bins='log', gridsize=30, mincnt=1, cmap=pyplot.cm.Spectral_r, zorder=1)
	cb = pyplot.colorbar(image, spacing='uniform', extend='max')
	
	
	pyplot.savefig(labels['filename'])

def plot_hist(table, m8_filename, no_uniq_genomes):
	'''Plots histogram for the genome hits per gene'''
	labels = {'ylabel': 'Centroids',\
			  'xlabel':'Genomes', \
			  'title':'Metagenome vs. Genome Prioritization',\
			  'filename': m8_filename+'_prioritizationHIST.pdf'}
	
	all_genes = []
	for gene in table:
		all_genes +=[len(table[gene])/float(no_uniq_genomes)]
	
	print 'No. of unique genomes: '+str(no_uniq_genomes)

	pyplot.figure()
	pyplot.xlabel(labels['xlabel'])
	pyplot.ylabel(labels['ylabel'])
	pyplot.title(labels['title'])
	pyplot.legend('All Centroids')
	pyplot.hist(all_genes, log=True, bins=30, edgecolor='white', color='gray')
	pyplot.savefig(labels['filename'])
	pyplot.savefig(labels['filename']+'.png')

	pyplot.figure()
	pyplot.xlabel('Genomes Fraction [Genome hits/Total no. of unique genomes in niche]')
	pyplot.ylabel('Centroids coverage')
	pyplot.title(labels['title'])
	pyplot.legend('All Centroids')
	hist, bins = numpy.histogram(all_genes, bins=200)
	offset = bins[1:]-bins[:-1]
	pyplot.plot(bins[:-1]+offset, numpy.cumsum(hist), color='gray')
	pyplot.savefig(labels['filename']+'_cumsum.pdf')
	pyplot.savefig(labels['filename']+'_cumsum.png')

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input_file', help='Gene Genomes blast results', required=True)
	parser.add_argument('--parsed', default=False, action='store_true', help='Input file is parsed')
	parser.add_argument('--metagenome_fasta', help='Metagenome FASTA file')
	parser.add_argument('--bypass_hist', default=False, action='store_true', help='Generates Histogram')
	parser.add_argument('--bypass_scatter', default=False, action='store_true', help='Generates Scatterplot')

	args = parser.parse_args()

	m8_filename = args.input_file

	if not args.parsed:
		try:
			fasta_filename = args.metagenome_fasta
		except:
			raise Exception('Metagenome fasta not entered')
		table = parse_table(m8_filename, fasta_filename)
		with open(m8_filename+'_parsed.m8','w') as foo:
			for i in table:
				for j in table[i]:
					foo.writelines('\t'.join([i, j])+'\n')
	else:
		table = read_parsed(m8_filename)

	uniq_genomes = []
	for gene in table:
		for genome in table[gene]:
			if genome not in uniq_genomes:
				uniq_genomes +=[genome]

	no_uniq_genomes = len(uniq_genomes)

	if not args.bypass_scatter:
		plot_scatter(table, m8_filename, no_uniq_genomes)
	if not args.bypass_hist:
		plot_hist(table, m8_filename, no_uniq_genomes)

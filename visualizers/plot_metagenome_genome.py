import os
import sys
import matplotlib
import re
import numpy
import pdb
import time
import numpy
import argparse

from matplotlib import pyplot
from src import create_fastas

def parse_table(m8_filename, fasta_filename):
	fasta_dict = create_fastas.read_fasta(fasta_filename)
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
	table = {}
	foo = open(m8_filename)

	for line in foo:
		split_i = [i.strip() for i in line.split('\t')]
		try:
			table[split_i[0]] += [split_i[1]]
		except:
			table[split_i[0]] = [split_i[1]]
	return table

def plot_scatter(table, m8_filename):
	labels = {'xlabel': 'Prioritized Centroids',\
			  'ylabel':'No. of Genomes (Log10)', \
			  'title':'Metagenome vs. Genome Prioritization',\
			  'filename': m8_filename+'_prioritizationScatter.pdf',\
			  'fname_split': m8_filename+'_prioritizationScatter_split.pdf'}	
	# uniref = []
	# unannotated = []
	all_genes = []
	for gene in table:
		all_genes +=[len(table[gene])/4189.0]
		# if gene.startswith('UniRef90'):
		# 	uniref +=[len(table[gene])/4189.0]
		# else:
		# 	unannotated += [len(table[gene])/4189.0]
	# uniref.sort()
	# unannotated.sort()
	# all_genes = uniref+unannotated
	all_genes.sort()
	pyplot.figure()
	pyplot.xlabel(labels['xlabel'])#'Alpha Prevalence_'+keys[i])
	pyplot.ylabel(labels['ylabel'])
	pyplot.title(labels['title']) #'Mean abundance')
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
	#pyplot.axhline(y=1, alpha=0.5, color='gray', zorder=5, hold=True, label='y=1')
	pyplot.savefig(labels['filename'])

	# pyplot.figure()
	# pyplot.xlabel(labels['xlabel'])#'Alpha Prevalence_'+keys[i])
	# pyplot.ylabel(labels['ylabel'])
	# pyplot.title(labels['title']) #'Mean abundance')

	# pyplot.scatter(numpy.arange(len(uniref)), \
	# 			   numpy.log(uniref), \
	# 			   c='grey', \
	# 			   alpha=0.1, \
	# 			   linewidths=0.0, \
	# 			   zorder=1, \
	# 			   marker='o',\
	# 			   label='UniRef90 Centroids')
	# pyplot.hold(True)
	# pyplot.scatter(numpy.arange(len(unannotated)), \
	# 			   numpy.log(unannotated), \
	# 			   c='red', \
	# 			   alpha=0.1, \
	# 			   linewidths=0.0, \
	# 			   zorder=2, \
	# 			   marker='o',\
	# 			   label='Unannotated Centroids')
	# pyplot.legend( loc=4, \
	# 			   fontsize='x-small', \
	# 			   framealpha=0.4, )	
	# #pyplot.axhline(y=1, alpha=0.5, color='gray', zorder=5, hold=True, label='y=1')
	# pyplot.savefig(labels['fname_split'])

def plot_hexbin(table, m8_filename):
	labels = {'xlabel': 'Prioritized Centroids',\
			  'ylabel':'No. of Genomes', \
			  'title':'Metagenome vs. Genome Prioritization',\
			  'filename': m8_filename+'_prioritizationHEXBIN.pdf'}
	pyplot.figure()
	pyplot.xlabel(labels['xlabel'])#'Alpha Prevalence_'+keys[i])
	pyplot.ylabel(labels['ylabel'])
	pyplot.title(labels['title']) #'Mean abundance')

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

def plot_hist(table, m8_filename):
	labels = {'ylabel': 'Centroids',\
			  'xlabel':'Genomes', \
			  'title':'Metagenome vs. Genome Prioritization',\
			  'filename': m8_filename+'_prioritizationHIST.pdf',
			  'fname_split':m8_filename+'_prioritizationHIST_split.pdf'}
	all_genes = []
	for gene in table:
		all_genes +=[len(table[gene])/4189.0]
	
	# uniref = []
	# unannotated = []
	# for gene in table:
	# 	if gene.startswith('UniRef90'):
	# 		uniref +=[len(table[gene])/4189.0]
	# 	else:
	# 		unannotated += [len(table[gene])/4189.0]
	# uniref.sort()
	# unannotated.sort()
	# all_genes = uniref+unannotated
	# all_genes.sort()

	# pyplot.figure()
	# pyplot.xlabel(labels['xlabel'])#'Alpha Prevalence_'+keys[i])
	# pyplot.ylabel(labels['ylabel'])
	# pyplot.title(labels['title']) #'Mean abundance')

	# pyplot.hist([numpy.array(uniref), numpy.array(unannotated)], log=True, bins=30, edgecolor='white')
	# pyplot.legend(['UniRef90 Centroids','Unannotated Centroids'])
	# pyplot.savefig(labels['fname_split'])

	pyplot.figure()
	pyplot.xlabel(labels['xlabel'])#'Alpha Prevalence_'+keys[i])
	pyplot.ylabel(labels['ylabel'])
	pyplot.title(labels['title']) #'Mean abundance')
	pyplot.legend('All Centroids')
	pyplot.hist(all_genes, log=True, bins=30, edgecolor='white')
	pyplot.savefig(labels['filename'])

	pyplot.figure()
	pyplot.xlabel('Genomes Fraction')#'Alpha Prevalence_'+keys[i])
	pyplot.ylabel('Centroids coverage')
	pyplot.title(labels['title']) #'Mean abundance')
	pyplot.legend('All Centroids')
	hist, bins = numpy.histogram(all_genes, bins=200)
	offset = bins[1:]-bins[:-1]
	pyplot.plot(bins[:-1]+offset, numpy.cumsum(hist))
	pyplot.savefig(labels['filename']+'_cumsum.pdf')

if __name__ == '__main__':

	m8_filename = sys.argv[1]

	if sys.argv[2]=='--to_parse':
		fasta_filename = sys.argv[3]
		table = parse_table(m8_filename, fasta_filename)
		with open(m8_filename+'_parsed.m8','w') as foo:
			for i in table:
				for j in table[i]:
					foo.writelines('\t'.join([i, j])+'\n')
	else:
		table = read_parsed(m8_filename)


	plot_scatter(table, m8_filename)
	# plot_hexbin(table, m8_filename)
	plot_hist(table, m8_filename)

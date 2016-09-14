import os
import sys
import matplotlib
import re
import numpy
import time
import numpy
import argparse
import pdb
from matplotlib import pyplot
from . import utilities

'''Tmp file to parse results'''
'''Analysis of Genome Hits per gene showing how many genomes each gene is found in.'''

def parse_table(m8_filename, fasta_filename):
	'''Parse the BLAST results to give gene hits to genomes
	Input: 
	m8_filename = filename of blast results
	fasta_filename = filename of corresponding fasta file
	
	Output: 
	table = {gene: [List of genomes]}'''

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
			split_sp = split_i[1].split('|')
			sp = [i for i in split_sp if 'g__' in i and '.s__' in i]
			if split_i[0] not in table:
				table[split_i[0]] = sp
			elif sp not in table[split_i[0]]:
				table[split_i[0]] += sp
	return table

def read_parsed(m8_filename):
	'''Read parsed table for {gene: genomes}
	Input: 
	m8_filename = filename of blast results

	Output: 
	table = {gene: [List of genomes]}'''

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
	'''Plots Scatter plot for genome hits per gene
	Input:
	m8_filename = filename of blast results
	table = {gene: [List of genomes]}
	no_uniq_genomes = Number of Unique Genomes in Metagenomic niche'''

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
	'''Plots HexBin plots for the genome hits per gene
	Input:
	m8_filename = filename of blast results
	table = {gene: [List of genomes]}'''

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
	'''Plots histogram for the genome hits per gene
	Input:
	m8_filename = filename of blast results
	table = {gene: [List of genomes]}
	no_uniq_genomes = Number of Unique Genomes in Metagenomic niche'''

	labels = {'ylabel': 'Centroids',\
			  'xlabel':'Genomes', \
			  'title':'Metagenome vs. Genome Prioritization',\
			  'filename': m8_filename+'_prioritizationHIST.pdf'}
	
	all_genes = []
	for gene in table:
		all_genes +=[len(table[gene])/float(no_uniq_genomes)]

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

def read_abund_prev(filename):
	keys = {'abund':0, 'alpha':0,'beta':0}
	abund = []
	alpha = []
	genes = []
	with open(filename) as foo:
		for line in foo:
			split_line = [re.sub('[\r\t\n]','', i) for i in line.split('\t')]
			if line.startswith('Centroids'):
				for i, val in enumerate(split_line):
					if 'abund' in val:
						keys['abund'] = i
					elif 'alpha' in val:
						keys['alpha'] = i
					elif 'beta' in val:
						keys['beta'] = i
			else:
				# pdb.set_trace()
				genes +=[split_line[0]]
				abund +=[float(split_line[keys['abund']])]
				alpha +=[float(split_line[keys['alpha']])]
	abund_prev = {'genes': genes, 'abundance': abund, 'prevalence': alpha}
	return abund_prev

def plot_metagenomic_priority(abund_prev, table, no_uniq_genomes, filename):
	mp_gp = {}
	genes = abund_prev['genes']
	abund = abund_prev['abundance']
	alpha = abund_prev['prevalence']

	abund = numpy.array(abund)/max(abund)
	alpha = numpy.array(alpha)/max(alpha)
	
	gp = []
	mp = []
	for i in range(len(genes)):
		gene = genes[i]
		try:
			gp += [len(table[gene])]
		except:
			gp += [0]
		mp += [min((abund[i], alpha[i]))]
	gp = numpy.array(gp)/float(no_uniq_genomes)
	pyplot.figure()
	pyplot.xlabel('Genomic Priority')
	pyplot.ylabel('Metagenomic Priority')
	pyplot.title('Metagenomic vs. Genomic Priority') 
	pyplot.scatter(numpy.log(gp), \
				   numpy.log(mp), \
				   # c='grey', \
				   alpha=0.1, \
				   linewidths=0.0, \
				   zorder=1, \
				   marker='o',\
				   label='All Centroids')
	pyplot.savefig(filename+'_mp_gp.pdf')
	pyplot.savefig(filename+'_mp_gp.png')
	
	pyplot.figure()
	pyplot.xlabel('Genomic Priority')
	pyplot.ylabel('Metagenomic Priority')
	pyplot.title('Metagenomic vs. Genomic Priority') 
	pyplot.scatter(gp, \
				   mp, \
				   # c='grey', \
				   alpha=0.1, \
				   linewidths=0.0, \
				   zorder=1, \
				   marker='o',\
				   label='All Centroids')
	pyplot.savefig(filename+'_nonlog_mp_gp.pdf')
	pyplot.savefig(filename+'_nonlog_mp_gp.png')

	pyplot.figure()
	pyplot.xlabel('Genomic Priority')
	pyplot.ylabel('Metagenomic Priority')
	pyplot.title('Metagenomic vs. Genomic Priority') 
	pyplot.hexbin(numpy.log(gp), \
				   numpy.log(mp), \
				   cmap='Blues',
				   gridsize=10)
	pyplot.colorbar()
	pyplot.savefig(filename+'_hexplot_mp_gp.pdf')
	pyplot.savefig(filename+'_hexplot_mp_gp.png')
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input-file', dest='input_file',  help='Gene Genomes blast results', required=True)
	parser.add_argument('--bypass-parse', dest= 'bypass_parse' , default=False, action='store_true', help='Input file is parsed')
	parser.add_argument('--parse-only', dest= 'parse_only' , default=False, action='store_true', help='To only parse')
	parser.add_argument('--metagenome-fasta',dest='metagenome_fasta'  ,help='Metagenome FASTA file')
	parser.add_argument('--write-no-genomes',dest= 'write_no_genomes' ,default=False, action='store_true', help='Write Gene to No. of genomes')
	parser.add_argument('--bypass-hist', dest= 'bypass_hist', default=False, action='store_true', help='Generates Histogram')
	parser.add_argument('--bypass-scatter',dest= 'bypass_scatter' ,  default=False, action='store_true', help='Generates Scatterplot')
	parser.add_argument('--abund-prev', dest= 'abund_prev' , default=False, help='Centroid prevalence and abundance')
	parser.add_argument('--bypass-priority-scatter',dest= 'bypass_priority_scatter',  default=False, action='store_true', help='Generates Scatterplot')

	args = parser.parse_args()

	m8_filename = args.input_file

	if not args.bypass_parse:
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
	
	if args.write_no_genomes:
		with open(m8_filename+'_no_genomes.m8','w') as foo:
			for i in table:
				foo.writelines('\t'.join([i, str(len(table[i]))])+'\n')
	
	uniq_genomes = []
	for gene in table:
		for genome in table[gene]:
			if genome not in uniq_genomes:
				uniq_genomes +=[genome]

	no_uniq_genomes = len(uniq_genomes)
	if not args.bypass_priority_scatter:
		abund_prev= read_abund_prev(args.abund_prev)
		plot_metagenomic_priority(abund_prev, table, no_uniq_genomes, args.input_file)
	print 'No. of unique genomes: '+str(no_uniq_genomes)

	if args.parse_only:
		sys.exit('Input files parsed: '+args.input_file)
	
	if not args.bypass_scatter:
		plot_scatter(table, m8_filename, no_uniq_genomes)
	if not args.bypass_hist:
		plot_hist(table, m8_filename, no_uniq_genomes)
	
if __name__ == '__main__':
	main()
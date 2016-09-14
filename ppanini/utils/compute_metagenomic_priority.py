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

'''Analysis of Genome Hits per gene showing how many genomes each gene is found in.'''

def parse_chocophlan_table(m8_filename):
	'''Parse the BLAST results to give gene hits to genomes
	Input: 
	m8_filename = filename of blast results
	fasta_filename = filename of corresponding fasta file
	
	Output: 
	table = {gene: [List of genomes]}'''

	table = {}
	foo = open(m8_filename)
	all_species = []
	for line in foo:
		if not line.startswith('#'):
			split_i = line.split('\t')
			sp = [re.sub('[\r\t\n]','',i) for i in split_i[1:]]
			try:
				table[split_i[0]] = sp
			except:
				table[split_i[0]] += sp
			all_species += sp
	no_uniq_genomes = len(set(all_species))
	return [table, no_uniq_genomes]

def parse_abundance_table(table_name, cols):
	table_obj = open(table_name)
	genes = []
	prevalence = []
	abundance = []

	for line in table_obj:
		split_i = [re.sub('[\r\t\n]','',i).strip() for i in line.split('\t')]
		if line.startswith('#'):
			keys = [split_i[cols[0]], split_i[1]]
		else:
			genes += [split_i[0]]
			prevalence +=[float(split_i[cols[0]])]
			abundance +=[float(split_i[cols[1]])]
	return [keys, genes, prevalence, abundance]

def compute_priority(genes, prevalence, abundance, genes_genomes, no_uniq_genomes):
	mp = []
	gp = []
 	max_prev = max(prevalence)
 	max_abund = max(abundance)
	for i, gene in enumerate(genes):
		mp += [min((prevalence[i]/max_prev, abundance[i]/max_abund))]
		try:
			gp += [len(genes_genomes[gene])/float(no_uniq_genomes)]
		except:
			gp += [0]
	return [mp, gp]
	
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i1', '--input_prab', help='Tabulated Genes Prevalence/Abundance Table', required=True)
	parser.add_argument('-i2', '--input_choco', help='Tabulated genes to genomes table', required=True)
	parser.add_argument('--cols', default='1,2', help='Columns to plot: x,y')
	
	args = parser.parse_args()
	
	cols = [int(i) for i in args.cols.split(',')]
	[genes_genomes, no_uniq_genomes] = parse_chocophlan_table(args.input_choco)
	[keys, genes, prevalence, abundance] = parse_abundance_table(args.input_prab, cols)

	[mp,gp] = compute_priority(genes, prevalence, abundance, genes_genomes, no_uniq_genomes)

	print '\t'.join(['#ID','Metagenomic_Priority', 'Genomic_Priority'])
	for i, gene in enumerate(genes):
		print '\t'.join([gene, str(mp[i]), str(gp[i])])

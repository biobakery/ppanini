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
from src import utilities

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
			if '|' in split_i[1]:
				split_sp = split_i[1].split('|')
				sp = [i for i in split_sp if 'g__' in i and '.s__' in i]
			else:
				sp = [re.sub('[\r\t\n]','',i) for i in split_i[1].split('\t')]
			try:
				if split_i[0] not in table:
					table[split_i[0]] = sp
					all_species += sp
				elif sp[0] not in table[split_i[0]]:
					table[split_i[0]] += sp
					all_species += sp
			except:
				pdb.set_trace()
	no_uniq_genomes = len(set(all_species))
	print str(no_uniq_genomes)+'here'
	return [table, no_uniq_genomes]

def parse_abundance_table(table_name, cols):
	table_obj = open(table_name)
	genes_prab = {}
	keys = []
	max1 = []
	max2 = []
	for line in table_obj:
		split_i = [re.sub('[\r\t\n]','',i).strip() for i in line.split('\t')]
		if line.startswith('#'):
			keys = [split_i[cols[0]], split_i[1]]
		else:
			genes_prab[split_i[0]] = tuple([float(split_i[cols[0]]), float(split_i[cols[1]])]) #(prevalence, abundance) always <--
			max1 += [float(split_i[cols[0]])]
			max2 += [float(split_i[cols[1]])]
	max_s = [max(max1), max(max2)]
	return [genes_prab, keys, max_s]

def compute_priority(genes_prab, genes_genomes, no_uniq_genomes, max_s):
	mp_gp = {}
	# max_genomes = 0
	# for gene in genes_prab:
	# 	i =len(genes_genomes[gene])
	# 	if i > max_genomes:
	# 		max_genomes = i
	# print i
	for gene in genes_prab:
		mp = min((genes_prab[gene][0]/max_s[0], genes_prab[gene][1]/max_s[1]))
		try:
			gp = len(genes_genomes[gene])/float(no_uniq_genomes)
			mp_gp[gene] = (mp, gp)
		except:
			mp_gp[gene] = (mp, 0)
	return mp_gp
	
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i1', '--input_prab', help='Prevalence/Abundance Table', required=True)
	parser.add_argument('-i2', '--input_choco', help='Chocophlan parsed Table', required=True)
	parser.add_argument('--cols', default='1,2', help='Columns to plot: x,y')
	
	args = parser.parse_args()
	
	cols = [int(i) for i in args.cols.split(',')]
	[genes_genomes, no_uniq_genomes] = parse_chocophlan_table(args.input_choco)
	[genes_prab, keys, max_s] = parse_abundance_table(args.input_prab, cols)

	mp_gp = compute_priority(genes_prab, genes_genomes, no_uniq_genomes, max_s)

	print '\t'.join(['#ID','Metagenomic_Priority', 'Genomic_Priority'])
	for gene in mp_gp:
		print '\t'.join([gene]+[str(i) for i in mp_gp[gene]])

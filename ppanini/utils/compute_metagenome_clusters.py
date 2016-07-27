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

def get_clusters_dict(gene_centroid_clusters_file_path):
	'''Return dict containing clusters
	Input: filepath to centroid file
	Output: centroid_gis (dict) {gene_centroid: [List of genes], }'''

	cluster_txt = open(gene_centroid_clusters_file_path)
	centroid_gis = {}
	for line in cluster_txt:
		if line.startswith('H'):
			split_i = [re.sub('[\r\t\n]', '', i) for i in line.split('\t')[-2:]]
			try:
				centroid_gis[split_i[1]] += [split_i[0]]
			except KeyError:
				centroid_gis[split_i[1]] = [split_i[0], split_i[1]]
	return centroid_gis

def get_uniref_dict(gene_centroid_clusters_file_path):
	'''Return dict containing clusters; Supply ppanini file here
	Input: filepath to centroid file
	Output: centroid_gis (dict) {gene_centroid: [List of genes], }'''

	cluster_txt = open(gene_centroid_clusters_file_path)
	centroid_gis = {}
	for line in cluster_txt:
		if not line.startswith('#'):
			split_i = [re.sub('[\r\t\n]', '', i).strip() for i in line.split('\t')[0].split('|')]
			if not split_i[1]=='UniRef90_unknown':
				try:
					centroid_gis[split_i[1]] += [split_i[0]]
				except KeyError:
					centroid_gis[split_i[1]] = [split_i[0]]
	return centroid_gis

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
			split_sp = split_i[1].split('|')
			sp = [i for i in split_sp if 'g__' in i and '.s__' in i]
			try:
				if split_i[0] not in table:
					table[split_i[0]] = sp
					all_species += sp
				elif sp[0] not in table[split_i[0]]:
					table[split_i[0]] += sp
					all_species += sp
			except:
				pdb.set_trace()
	no_uniq_genomes = len(all_species)
	# print str(no_uniq_genomes)+'here'
	return [table, no_uniq_genomes]

foo1 = sys.argv[1] #chocophlan
foo2 = sys.argv[2] #gc_gi
foo3 = sys.argv[3] #ppanini_table

[table, no_uniq_genomes] = parse_chocophlan_table(foo1)
uniref_gis = get_uniref_dict(foo3)
centroid_gis = get_clusters_dict(foo2)

sub_chocophlan = {}
for gene in uniref_gis:
	genes = uniref_gis[gene]
	for gi in genes:
		try:
			sub_chocophlan[gene] += table[gi]
		except:
			sub_chocophlan[gene] = table[gi]
	try:
		sub_chocophlan[gene] = list(set(sub_chocophlan[gene]))
	except:
		pdb.set_trace()

for gene in centroid_gis:
	genes = centroid_gis[gene]
	for gi in genes:
		try:
			sub_chocophlan[gene] += table[gi]
		except:
			sub_chocophlan[gene] = table[gi]
	try:
		sub_chocophlan[gene] = list(set(sub_chocophlan[gene]))
	except:
		pdb.set_trace()

for gene in sub_chocophlan:
	if not sub_chocophlan[gene]:
		print '\t'.join([gene,'UniRef90_unknown'])
	else:
		print '\t'.join([gene]+sub_chocophlan[gene])



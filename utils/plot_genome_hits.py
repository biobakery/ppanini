import os
import sys
import matplotlib
import re
import numpy
import pdb
import time
import numpy
import argparse
import random

from matplotlib import pyplot
from matplotlib import colors

'''Analysis of genomes in the niche giving a sense of what genomes are in the community'''

numpy.seterr(divide='ignore', invalid='ignore')

def read_parsed(m8_filename, go_table):
	'''Returns genomes dict {genome: {'UniRef90': #hits, 'UniRef90_NA': #hits, 'NA': #hits}}'''
	table = {}
	foo = open(m8_filename)
	for line in foo:
		split_i = [i.strip() for i in line.split('\t')]
		if 'UniRef90' in split_i[0]:
			if 'NA' == go_table[split_i[0]]:
				key = 'UniRef90'
			else:
				key = 'UniRef90_NA'
		else:
			key = 'NA'
		try:
			table[split_i[1]][key] += 1
		except:
			table[split_i[1]] = {'UniRef90':0, 'UniRef90_NA':0, 'NA':0}
			table[split_i[1]][key] += 1
	return table

def read_go_map(m8_filename):
	'''Returns dict {gene: GO annotation}'''
	table = {}
	foo = open(m8_filename)

	for line in foo:
		split_i = [i.strip() for i in line.split('\t')]
		table[split_i[0]] = split_i[1]
	return table

def plot_scatter(table, m8_filename):
	'''Plots a scatter plot for genome hits over prioritized genes'''
	labels = {'xlabel': 'No. of Prioritized genes',\
			  'ylabel':'Frequency(No. of Genomes)', \
			  'title':'Metagenome vs. Genome Prioritization',\
			  'filename': m8_filename+'_genomehits.pdf'}
	genomes = []
	hits = []

	for genome in table:
		hits += [sum(table[genome].values())]
	hits.sort()
	inds = numpy.arange(len(hits))

	pyplot.figure()
	pyplot.hist(hits, color='b', bins=1000)
	pyplot.xlabel(labels['xlabel'])
	pyplot.ylabel(labels['ylabel'])
	pyplot.title(labels['title'])
	pyplot.savefig(labels['filename'])

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input_file', help='Gene Genomes blast results parsed**', required=True)
	parser.add_argument('--map', help='Gene to GO mapper', required=True)
	parser.add_argument('--bypass_scatter', default=False, action='store_true', help='Sctter plot for genomes')
	parser.add_argument('--bypass_stats', default=False, action='store_true', help='Write stats for genome gene hits')
	parser.add_argument('--bypass_graphlan_rings', default=False, action='store_true', help='Generates graphlan rings file')
	args = parser.parse_args()

	m8_filename = args.input_file
	go_map = read_go_map(args.map)
	genomes = read_parsed(m8_filename, go_map)
	all_values = []

	if not args.bypass_scatter:
		plot_scatter(genomes, m8_filename)

	if not args.bypass_stats:
		with open(m8_filename+'_allstats.txt', 'w') as foo:
			for gene in genomes:
				x = genomes[gene].values()
				all_values += x
				foo.writelines([gene+'\t'+str(sum(x))+'\n'])
	max_val = float(max(all_values))
	
	if not args.bypass_graphlan_rings:
		##GRAPHLAN RINGS FILE
		ref = numpy.arange(int(max_val)+1)
		print ref #Reference
		ref_i = 0
		with open(sys.argv[1]+'_allrings.txt','w') as foo_rings:
			foo_rings.writelines(['ring_internal_separator_thickness\t1\t1.0\n'])# #UniRe90+GO
			foo_rings.writelines(['ring_internal_separator_thickness\t2\t1.0\n'])# #UniRef90NA
			foo_rings.writelines(['ring_internal_separator_thickness\t3\t1.0\n']) #SRS
			foo_rings.writelines(['ring_internal_separator_thickness\t4\t1.0\n'])#REFERENCE

			for genome in genomes:
				while ref_i <= len(refs):
					foo_rings.writelines([genome+'\tring_height\t4\t'+str(inds[check])+'\n'])	
					foo_rings.writelines([genome+'\tring_color\t4\t#000000\n'])
					ref_i += 1
				foo_rings.writelines([genome+'\tring_height\t1\t'+str(numpy.log(genomes[genome]['UniRef90']))+'\n'])
				foo_rings.writelines([genome+'\tring_height\t2\t'+str(numpy.log(genomes[genome]['UniRef90_NA']))+'\n'])
				foo_rings.writelines([genome+'\tring_height\t3\t'+str(numpy.log(genomes[genome]['NA']))+'\n'])

				foo_rings.writelines([genome+'\tring_color\t1\t#0000FF\n'])
				foo_rings.writelines([genome+'\tring_color\t2\t#FF0000\n'])
				foo_rings.writelines([genome+'\tring_color\t3\t#FFFF00\n'])

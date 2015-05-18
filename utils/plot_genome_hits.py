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
from matplotlib import colors
from src import create_fastas

numpy.seterr(divide='ignore', invalid='ignore')

def read_parsed(m8_filename, go_table):
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
	table = {}
	foo = open(m8_filename)

	for line in foo:
		split_i = [i.strip() for i in line.split('\t')]
		table[split_i[0]] = split_i[1]
	return table

def plot_scatter(table, m8_filename):
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

	m8_filename = sys.argv[1]
	go_map = read_go_map(sys.argv[2])
	genomes = read_parsed(m8_filename, go_map)
	all_values = []
	plot_scatter(genomes, m8_filename)
	if sys.argv[3]=='--write_all':
		with open(m8_filename+'_allstats.txt', 'w') as foo:
			for gene in genomes:
				x = genomes[gene].values()
				all_values += x
				foo.writelines([gene+'\t'+str(sum(x))+'\n'])
	max_val = float(max(all_values))
	
	with open(sys.argv[1]+'_allrings.txt','w') as foo_rings:
		foo_rings.writelines(['ring_internal_separator_thickness\t1\t1.0\n'])#ring_width\t5\t0.5\n']) #UniRe90
		foo_rings.writelines(['ring_internal_separator_thickness\t2\t1.0\n'])#ring_width\t5\t0.5\n']) #UniRef90NA
		foo_rings.writelines(['ring_internal_separator_thickness\t3\t1.0\n'])#ring_width\t5\t0.5\n']) #SRS

		for genome in genomes:
			foo_rings.writelines([genome+'\tring_height\t1\t'+str(numpy.log(genomes[genome]['UniRef90']))+'\n'])
			foo_rings.writelines([genome+'\tring_height\t2\t'+str(numpy.log(genomes[genome]['UniRef90_NA']))+'\n'])
			foo_rings.writelines([genome+'\tring_height\t3\t'+str(numpy.log(genomes[genome]['NA']))+'\n'])
			
			foo_rings.writelines([genome+'\tring_color\t1\t#AAAA00\n'])#+str(genomes[genome]['UniRef90']/max_val)+'\n'])
			foo_rings.writelines([genome+'\tring_color\t2\t#00AAAA\n'])#+str(genomes[genome]['UniRef90_NA']/max_val)+'\n'])
			foo_rings.writelines([genome+'\tring_color\t3\t#AA00AA\n'])#+str(genomes[genome]['NA']/max_val)+'\n'])
	
	import matplotlib.patches as patches
	# cdict = {'white':[(0.0,1.0,1.0),(1.0,1.0,1.0),(1.0,1.0,1.0)], '#AAAA00':[(0.0,1.0,1.0),(1.0,0.0,0.0),(1.0,0.0,0.0)]}
	fig5 = pyplot.figure()
	ax5 = fig5.add_subplot(111, aspect='equal')
	all_values = numpy.array(all_values)
	x=numpy.arange(100)
	x=x/100.0
	x_i = 0.01
	for p in range(len(x)):
		# print x[p]
		ax5.add_patch(patches.Rectangle((x_i,0.01),0.01,0.05,alpha=x[p], facecolor='#AAAA00', edgecolor='none'))
		ax5.add_patch(patches.Rectangle((x_i,0.10),0.01,0.05,alpha=x[p], facecolor='#00AAAA', edgecolor='none'))
		ax5.add_patch(patches.Rectangle((x_i,0.20),0.01,0.05,alpha=x[p], facecolor='#AA00AA', edgecolor='none'))
		x_i +=0.01
	# pyplot.imshow([numpy.array(all_values)/max_val, numpy.array(all_values)/max_val], cmap=matplotlib.colors.LinearSegmentedColormap('cutom_XYZ', cdict, 256))
	# pyplot.colorbar()
	fig5.savefig(sys.argv[1]+'_cbar.pdf', dpi=300, bbox_inches='tight')

import pdb
import os
import sys
import matplotlib
import re
import numpy
import time
import numpy
import argparse
import random

from . import utils
from utils import plot_metagenome_genome
from . import utilities
from matplotlib import pyplot
from matplotlib import colors

'''Tmp file to parse results'''

'''Analysis of genomes in the niche giving a sense of what genomes are in the community'''

numpy.seterr(divide='ignore', invalid='ignore')

def read_parsed_genome(m8_filename, go_table):
	'''Returns genomes dict {genome: {'UniRef90': #hits, 'UniRef90_NA': #hits, 'NA': #hits}}
	Input:
	m8_filename = filename for parsed blast results
	go_table = {uniref_id : GO_ID, ...}

	Output:
	table = {genome: {UniRef90: #hits, UniRef90_NA: #hits, NA: #hits}, ...}'''

	table = {}
	foo = open(m8_filename)
	for line in foo:
		split_i = [i.strip() for i in line.split('\t')]
		if 'UniRef90' in split_i[0]:
			try:
				if 'NA' == go_table[split_i[0]]:
					key = 'UniRef90'
				else:
					key = 'UniRef90_NA'
			except:
				pdb.set_trace()
		else:
			key = 'NA'
		try:
			table[split_i[1]][key] += 1
		except:
			table[split_i[1]] = {'UniRef90':0, 'UniRef90_NA':0, 'NA':0}
			table[split_i[1]][key] += 1
	return table

def read_go_map(filename):
	'''Returns dict {gene: GO annotation}
	Input:
	filename = filename for Uniref90 to GO mapping

	Output:
	go_table = {uniref_id : GO_ID, ...}
	'''
	go_table = {}
	foo = open(filename)

	for line in foo:
		split_i = [i.strip() for i in line.split('\t')]
		go_table[split_i[0]] = split_i[1]

	return go_table

def plot_scatter(table, m8_filename):
	'''Plots a scatter plot for genome hits over prioritized genes
	Input:
	table = {genome: {UniRef90: #hits, UniRef90_NA: #hits, NA: #hits}, ...}
	m8_filename = filename for the output figure'''

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
	pyplot.hist(hits, color='gray', bins=30, edgecolor='none')
	pyplot.xlabel(labels['xlabel'])
	pyplot.ylabel(labels['ylabel'])
	pyplot.title(labels['title'])
	pyplot.savefig(labels['filename'])
	pyplot.savefig(labels['filename']+'.png')

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input_file', help='Gene Genomes blast results parsed**', required=True)
	parser.add_argument('--map', help='Gene to GO mapper from ppanini_visualizer', required=True)
	parser.add_argument('--bypass-scatter', dest ='bypass_scatter' , default=False, action='store_true', help='Scatter plot for genomes')
	parser.add_argument('--bypass-write-stats', dest = 'bypass_write_stats', default=False, action='store_true', help='Write stats for genome gene hits')
	parser.add_argument('--bypass-graphlan-rings',dest = 'bypass_graphlan_rings', default=False, action='store_true', help='Generates graphlan rings file')
	parser.add_argument('--pangenome-size', dest = 'pangenome_size', default=False, help='Pangenome size mapping file')
	parser.add_argument('--metagenome-fasta', dest= 'metagenome_fasta', help='Metagenome FASTA file')
	parser.add_argument('--bypass-parse', dest= 'bypass_parse', default=False, action='store_true', help='Input file is parsed')

	args = parser.parse_args()

	m8_filename = args.input_file
	go_map = read_go_map(args.map)
	
	if not args.bypass_parse:
		try:
			fasta_filename = args.metagenome_fasta
		except:
			raise Exception('Metagenome fasta not entered')
		table = plot_metagenome_genome.parse_table(m8_filename, fasta_filename)
		with open(m8_filename+'_parsed.m8','w') as foo:
			for i in table:
				for j in table[i]:
					foo.writelines('\t'.join([i, j])+'\n')
		m8_filename = m8_filename+'_parsed.m8'
		genomes = read_parsed_genome(m8_filename, go_map)
	else:
		genomes = read_parsed_genome(m8_filename, go_map) 

	all_values = []
	
	pgsize = utilities.read_dict_num(args.pangenome_size) 
	pangenome = {}
	for genome in pgsize:
		name = '.'.join(genome.split('.')[:2])
		pangenome[name] = pgsize[genome]	

	if not args.bypass_scatter:
		plot_scatter(genomes, m8_filename)


	if not args.bypass_write_stats:
		with open(m8_filename+'_allstats.txt', 'w') as foo:
			for gene in genomes:
				x = genomes[gene].values()
				all_values += x
				foo.writelines([gene+'\t'+str(sum(x))+'\n'])

		with open(m8_filename+'_allstats_stratified.txt', 'w') as foo:
			for gene in genomes:
				x = genomes[gene].values()
				all_values += x
				foo.writelines([gene+'\t'+\
								str(genomes[gene]['UniRef90'])+'\t'+\
								str(genomes[gene]['UniRef90_NA'])+'\t'+\
								str(genomes[gene]['NA'])+'\n'])

	# all_values = numpy.log(numpy.array(all_values))
	# max_val = int(max(all_values))
	# print max_val

	##GRAPHLAN RINGS FILE

	if not args.bypass_graphlan_rings:
		frac = {'UniRef90': {}, 'UniRef90_NA': {}, 'NA':{}}
		for genome in genomes:
				split_genome = genome.split('.')
				g_ind = [i for i in range(len(split_genome)) if 'g__' in split_genome[i]][0]
				tmp_genome = '.'.join(split_genome[g_ind:g_ind+2])
				if tmp_genome not in pangenome:
					pdb.set_trace()
					raise Exception('Genome '+genome+' not found in PANGENOME!!!')
				pg_i = pangenome[tmp_genome]
				frac['UniRef90'][genome] = genomes[genome]['UniRef90']/float(pg_i)
				frac['UniRef90_NA'][genome] = genomes[genome]['UniRef90_NA']/float(pg_i)
				frac['NA'][genome] = genomes[genome]['NA']/float(pg_i)

		u90_vals = [i for i in frac['UniRef90'].values() if not i == 0.0]
		u90na_vals = [i for i in frac['UniRef90_NA'].values() if not i == 0.0]
		na_vals = [i for i in frac['NA'].values() if not i == 0.0]

		factor = {}
		if u90_vals:
			factor['UniRef90'] = numpy.average(u90_vals)
		else:
			factor['UniRef90'] = 1 
		
		if na_vals:
			factor['UniRef90_NA'] = numpy.average(u90na_vals)
		else:
			factor['UniRef90_NA'] = 1 
			
		if na_vals:
			factor['NA'] = numpy.average(na_vals)
		else:
			factor['NA'] = 1 
		

		genomes_frac = {}
		for i in frac:
			for genome in frac[i]:
				if not genome in genomes_frac:
					genomes_frac[genome] = {i: frac[i][genome]/factor[i]}
				else:
					genomes_frac[genome][i] = frac[i][genome]/factor[i]

		# ref = numpy.arange(max_val+1)
		# print ref #Reference
		ref_i = 0
		with open(m8_filename+'_allrings.txt','w') as foo_rings:
			foo_rings.writelines(['ring_internal_separator_thickness\t1\t1.0\n'])# #UniRe90+GO
			foo_rings.writelines(['ring_internal_separator_thickness\t2\t1.0\n'])# #UniRef90NA
			foo_rings.writelines(['ring_internal_separator_thickness\t3\t1.0\n']) #SRS
			foo_rings.writelines(['ring_internal_separator_thickness\t4\t1.0\n'])#REFERENCE		
			
			ref = []
			for genome in genomes_frac:
				# x = [genomes_frac[genome]['UniRef90'], \
				# 	 genomes_frac[genome]['UniRef90_NA'],\
				# 	 genomes_frac[genome]['NA']]

				x = [numpy.log(genomes_frac[genome]['UniRef90']), \
					 numpy.log(genomes_frac[genome]['UniRef90_NA']),\
					 numpy.log(genomes_frac[genome]['NA'])]
				for i, val in enumerate(x):
					if -1*numpy.inf==val or numpy.inf==val:
						x[i] = 0
					else:
						x[i] = abs(x[i])
				ref +=x
				foo_rings.writelines([genome+'\tring_height\t1\t'+str(x[0])+'\n'])
				foo_rings.writelines([genome+'\tring_height\t2\t'+str(x[1])+'\n'])
				foo_rings.writelines([genome+'\tring_height\t3\t'+str(x[2])+'\n'])

				foo_rings.writelines([genome+'\tring_color\t1\t#0000FF\n'])
				foo_rings.writelines([genome+'\tring_color\t2\t#FF0000\n'])
				foo_rings.writelines([genome+'\tring_color\t3\t#FFFF00\n'])
			max_val = max(ref)
			# pdb.set_trace()
			ref = numpy.arange(max_val+1)
			print ref
			for genome in genomes_frac:
				foo_rings.writelines([genome+'\tring_height\t4\t'+str(ref[ref_i])+'\n'])	
				foo_rings.writelines([genome+'\tring_color\t4\t#000000\n'])
				ref_i +=1
				if ref_i >= len(ref):
					break
	
if __name__ == '__main__':
	main()

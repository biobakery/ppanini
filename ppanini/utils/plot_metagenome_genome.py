import os
import sys
import matplotlib
import re
import numpy
import time
import numpy
import argparse
import pdb
import matplotlib.pyplot as plt
from matplotlib import colors
#from . import utils

from .. import utilities

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
	print"Total No. of genes:", len(table) 
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
	plt.figure()
	plt.xlabel(labels['xlabel'])
	plt.ylabel(labels['ylabel'])
	plt.title(labels['title']) 
	plt.scatter(numpy.arange(len(all_genes)), \
				   numpy.log(all_genes), \
				   c='grey', \
				   alpha=0.1, \
				   linewidths=0.0, \
				   zorder=1, \
				   marker='o',\
				   label='All Centroids')
	'''plt.legend( loc=4, \
				   fontsize='x-small', \
				   framealpha=0.4, )
				   '''	
	plt.savefig(labels['filename'])
	plt.savefig(labels['filename']+'.png')

def plot_hexbin(table, m8_filename):
	'''Plots HexBin plots for the genome hits per gene
	Input:
	m8_filename = filename of blast results
	table = {gene: [List of genomes]}'''

	labels = {'xlabel': 'Prioritized Centroids',\
			  'ylabel':'No. of Genomes', \
			  'title':'Metagenome vs. Genome Prioritization',\
			  'filename': m8_filename+'_prioritizationHEXBIN.pdf'}
	plt.figure()
	plt.xlabel(labels['xlabel'])
	plt.ylabel(labels['ylabel'])
	plt.title(labels['title']) 

	all_genes = []
	for gene in table:
		all_genes +=[len(table[gene])/4189.0]
	all_genes.sort()
	plt.hold(True)
	image = plt.hexbin(numpy.arange(len(all_genes)), \
				   all_genes,\
				   bins='log', gridsize=30, mincnt=1, cmap= 'darkgoldenrod', zorder=1) 
	cb = plt.colorbar(image, spacing='uniform', extend='max')
	
	plt.savefig(labels['filename'])

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

	plt.figure()
	plt.xlabel(labels['xlabel'])
	plt.ylabel(labels['ylabel'])
	plt.title(labels['title'])
	plt.legend('All Centroids')
	plt.hist(all_genes, log=True, bins=30, edgecolor='white', color='gray')
	plt.savefig(labels['filename'])
	plt.savefig(labels['filename']+'.png')
	fig, axarr = plt.subplots()
	#plt.figure()
	#plt.xlabel('Genomes Fraction [Genome hits/Total no. of unique genomes in niche]')
	#plt.ylabel('Centroids coverage')
	#plt.title(labels['title'])
	#plt.legend('All Centroids')
	axarr.get_xaxis().set_tick_params(which='both', labelsize=8,top='off',  direction='out')
	axarr.get_yaxis().set_tick_params(which='both', labelsize=8, right='off', direction='out')
	axarr.yaxis.set_label_position('left') 
	axarr.set_ylabel('Centroids coverage')
	axarr.set_xlabel('Genomes Fraction [Genome hits/Total no. of unique genomes in niche]')
	plt.title(labels['title'])
	plt.legend('All Centroids') 
	
	
	hist, bins = numpy.histogram(all_genes, bins=200)
	offset = bins[1:]-bins[:-1]
	axarr.plot(bins[:-1]+offset, numpy.cumsum(hist), color='darkgoldenrod')
	plt.savefig(labels['filename']+'_cumsum.pdf')
	plt.savefig(labels['filename']+'_cumsum.png')

def read_abund_prev(filename):
	keys = {'abund':0, 'alpha':0,'beta':0}
	abund = []
	prev = []
	genes = []
	ppanini_score = []
	with open(filename) as foo:
		for line in foo:
			split_line = [re.sub('[\r\t\n]','', i) for i in line.split('\t')]
			if line.startswith('#'):
				for i, val in enumerate(split_line):
					if 'abundance' in val:
						keys['abundance'] = i
					elif 'prevalence' in val:
						keys['prevalence'] = i
					#elif 'beta' in val:
					#	keys['beta'] = i
					elif 'ppanini_score' in val:
						keys['ppanini_score'] = i	
			else:
				# pdb.set_trace()
				genes +=[split_line[0]]
				ppanini_score +=[float(split_line[keys['ppanini_score']])]
				abund +=[float(split_line[keys['abundance']])]
				prev +=[float(split_line[keys['prevalence']])]
	abund_prev = {'genes': genes, 'ppanini_score':ppanini_score, 'abundance': abund, 'prevalence': prev}
	return abund_prev
def read_data(mg_file, ppanini_output_file):
	metagenomic_table  = read_parsed(mg_file)
	uniq_genomes = []
	for gene in metagenomic_table:
		for genome in metagenomic_table[gene]:
			if genome not in uniq_genomes:
				uniq_genomes +=[genome]
	no_uniq_genomes = len(uniq_genomes)
	print 'No. of unique genomes: '+str(no_uniq_genomes)
	ppanini_output = read_abund_prev(ppanini_output_file)
	return metagenomic_table, ppanini_output, no_uniq_genomes 
def master_plot():
    data_scale = 'none'
    fig, axarr = plt.subplots(nrows=2, ncols=2, dpi=300)#, sharex=False, sharey=False)
    fig.set_size_inches(10, 10)
    metagenomic_table1, ppanini_output1, no_uniq_genomes1  = read_data('./PARSED_BLAST_RESULTS/stool_mg.m8', '/Users/rah/Documents/Hutlab/ppanini/output_tables/stool_table.txt')
    scatter_plot_metagenomic_priority(axarr[0, 0], ppanini_output1, metagenomic_table1, no_uniq_genomes1, title = 'Stool', scale = data_scale)
    
    metagenomic_table2, ppanini_output2, no_uniq_genomes2 = read_data('./PARSED_BLAST_RESULTS/AN_mg.m8', '/Users/rah/Documents/Hutlab/ppanini/output_tables/AN_table.txt')
    scatter_plot_metagenomic_priority(axarr[1, 1], ppanini_output2, metagenomic_table2, no_uniq_genomes2, title = 'Anterior nares', scale = data_scale)
    
    metagenomic_table3, ppanini_output3, no_uniq_genomes3 = read_data('./PARSED_BLAST_RESULTS/BM_mg.m8', '/Users/rah/Documents/Hutlab/ppanini/output_tables/BM_table.txt')
    scatter_plot_metagenomic_priority(axarr[1, 0], ppanini_output3, metagenomic_table3, no_uniq_genomes3, title = 'Buccal mucosa', scale = data_scale)
    
    metagenomic_table4, ppanini_output4, no_uniq_genomes4 = read_data('./PARSED_BLAST_RESULTS/PF_mg.m8', '/Users/rah/Documents/Hutlab/ppanini/output_tables/PF_table.txt')
    scatter_plot_metagenomic_priority(axarr[0, 1], ppanini_output4, metagenomic_table4, no_uniq_genomes4, title = 'Posterior fornix', scale = data_scale)
    
    # Add a colorbar
    #fig.subplots_adjust(right=1.25)
    '''cbar_ax = fig.add_axes([0.9, 0.037, 0.02, 0.863])
    fig.colorbar(im, cax = cbar_ax)'''
    #plt.axis('tight')
    
    #plt.margins(.75)
    #plt.xlim(xmin=-.15, xmax =5)
    #plt.ylim(xmin =1,ymin=1)
    #plt.subplots_adjust(top=0.19, right=0.99)
    
    #axarr.set_autoscale_on(True)
    #axarr.set_adjustable('box-forced')
    plt.tight_layout()#pad=0.4, w_pad=0.5, h_pad=1.0
    # Fine-tune figure; make subplots close to each other and hide x ticks for
    # all but bottom plot.
    #fig.subplots_adjust(hspace=0)
    #plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    #plt.setp([axarr], visible=False)
    st = fig.suptitle('Metagenomic vs. Genomic Priority', fontsize=12, fontweight='bold', va='top')
    
    # shift subplots up:
    st.set_y(0.95)
    fig.subplots_adjust(top=.9)
    fig.subplots_adjust()
    
    
    plt.savefig('metagenomic_genomic_priority_sctterplot_'+'sclae_'+str(data_scale)+'.pdf', pad_inches = .05, dpi=300) 
    fig, axarr = plt.subplots(nrows=2, ncols=2, dpi=300)#, sharex=False, sharey=False)
    fig.set_size_inches(10, 10)
    hexbin_plot_metagenomic_priority(axarr[0, 0], ppanini_output1, metagenomic_table1, no_uniq_genomes1, title = 'Stool', scale = data_scale)
    
    hexbin_plot_metagenomic_priority(axarr[1, 1], ppanini_output2, metagenomic_table2, no_uniq_genomes2, title = 'Anterior nares', scale = data_scale)
    
    hexbin_plot_metagenomic_priority(axarr[1, 0], ppanini_output3, metagenomic_table3, no_uniq_genomes3, title = 'Buccal mucosa', scale = data_scale)
    
    im= hexbin_plot_metagenomic_priority(axarr[0, 1], ppanini_output4, metagenomic_table4, no_uniq_genomes4, title = 'Posterior fornix', scale = data_scale)
    
    # Add a colorbar
    #fig.subplots_adjust(right=1.25)
    '''cbar_ax = fig.add_axes([0.9, 0.037, 0.02, 0.863])
    fig.colorbar(im, cax = cbar_ax)'''
    #plt.axis('tight')
    
    #plt.margins(.75)
    #plt.xlim(xmin=-.15, xmax =5)
    #plt.ylim(xmin =1,ymin=1)
    #plt.subplots_adjust(top=0.19, right=0.99)
    
    #axarr.set_autoscale_on(True)
    #axarr.set_adjustable('box-forced')
    plt.tight_layout()#pad=0.4, w_pad=0.5, h_pad=1.0
    # Fine-tune figure; make subplots close to each other and hide x ticks for
    # all but bottom plot.
    #fig.subplots_adjust(hspace=0)
    #plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
    #plt.setp([axarr], visible=False)
    st = fig.suptitle('Metagenomic vs. Genomic Priority', fontsize=12, fontweight='bold', va='top')
    
    # shift subplots up:
    st.set_y(0.95)
    fig.subplots_adjust(top=.9)
    fig.subplots_adjust()
    
    
    plt.savefig('metagenomic_genomic_priority_hexbitplot'+'sclae_'+str(data_scale)+'.pdf', pad_inches = .05, dpi=300) 

def scatter_plot_metagenomic_priority(axe, abund_prev, table, no_uniq_genomes, title,scale = None):
	if axe is None:
		axe = plt.gca()
	mp_gp = {}
	genes = abund_prev['genes']
	abund = abund_prev['abundance']
	prev = abund_prev['prevalence']
	ppanini_score = abund_prev['ppanini_score']
	abund = numpy.array(abund)/max(abund)
	prev = numpy.array(prev)/max(prev)
	
	gp = []
	mp = []
	for i in range(len(genes)):
		gene = genes[i]
		try:
			gp += [len(table[gene])]
		except:
			gp += [0]
		#mp += [min((abund[i], alpha[i]))]
		mp +=[ppanini_score[i]]
	gp = numpy.array(gp)/float(no_uniq_genomes)
	if scale == 'log':
		gp = numpy.log(gp)
		mp = numpy.log(mp)
	#fig, axe = plt.subplots()
	def ncolors( n, colormap="jet" ):
		"""utility for defining N evenly spaced colors across a color map"""
		cmap = plt.get_cmap( colormap )
		cmap_max = cmap.N
		return [cmap( int( k * cmap_max / (n - 1) ) ) for k in range( n )]
	my_color = ncolors(4)
	#x_dic= {'Gut': '#b87333', 'Skin':'#ffff00', 'Oral':'#009fff', 'Vaginal':'#ff4d00'}
	color_dic= {'Stool': my_color[0], 'Anterior nares':my_color[2], 'Buccal mucosa':my_color[1], 'Posterior fornix':my_color[3]}
	axe.scatter(gp, \
			   mp, \
			   c= color_dic[title],\
			   #'darkgoldenrod', \
			   ##'slategray'
			   alpha=0.1, \
			   linewidths=0.01, \
			   zorder=0, \
			   marker='o',\
			   label='All Centroids')
	if scale == 'log':
		axe.set_ylabel('Metagenomic Priority (log)', fontsize=10)
		axe.set_xlabel('Genomic Priority (log)', fontsize=10)
	else:
		axe.set_ylabel('Metagenomic Priority', fontsize=10)
		axe.set_xlabel('Genomic Priority')
	axe.get_xaxis().set_tick_params(which='both', labelsize=8,top='off',  direction='out')
	axe.get_yaxis().set_tick_params(which='both', labelsize=8, right='off', direction='out')
	axe.yaxis.set_label_position('left') 
	axe.set_title(title, fontsize=10, fontweight='bold')

def hexbin_plot_metagenomic_priority(axe, abund_prev, table, no_uniq_genomes, title, scale = None):
	if axe is None:
		axe = plt.gca()
	mp_gp = {}
	genes = abund_prev['genes']
	abund = abund_prev['abundance']
	prev = abund_prev['prevalence']
	ppanini_score = abund_prev['ppanini_score']
	abund = numpy.array(abund)/max(abund)
	prev = numpy.array(prev)/max(prev)
	
	gp = []
	mp = []
	for i in range(len(genes)):
		gene = genes[i]
		try:
			gp += [len(table[gene])]
		except:
			gp += [0]
		#mp += [min((abund[i], alpha[i]))]
		mp +=[ppanini_score[i]]
	gp = numpy.array(gp)/float(no_uniq_genomes)
	if scale == 'log':
		gp = numpy.log(gp)
		mp = numpy.log(mp)
	def ncolors( n, colormap="jet" ):
		"""utility for defining N evenly spaced colors across a color map"""
		cmap = plt.get_cmap( colormap )
		cmap_max = cmap.N
		return [cmap( int( k * cmap_max / (n - 1) ) ) for k in range( n )]
	my_color = ncolors(4)
	#x_dic= {'Gut': '#b87333', 'Skin':'#ffff00', 'Oral':'#009fff', 'Vaginal':'#ff4d00'}
	color_dic= {'Stool': my_color[0], 'Anterior nares':my_color[2], 'Buccal mucosa':my_color[1], 'Posterior fornix':my_color[3]}
	im  = axe.hexbin(gp, \
				   mp, \
				   cmap='YlOrBr',#'Blues',
				   gridsize=10)
	if scale == 'log':
		axe.set_ylabel('Metagenomic Priority (log)',  fontsize=10)
		axe.set_xlabel('Genomic Priority (log)', fontsize=10)
	else:
		axe.set_ylabel('Metagenomic Priority', fontsize=10)
		axe.set_xlabel('Genomic Priority', fontsize=10)
	axe.set_title(title, fontsize=10, fontweight='bold')#, va='top')
	axe.get_xaxis().set_tick_params(which='both', labelsize=8,top='off',  direction='out')
	axe.get_yaxis().set_tick_params(which='both', labelsize=8, right='off', direction='out')
	axe.yaxis.set_label_position('left') 
	return im
def main():
	master_plot()
	return
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input-file', dest='input_file',  help='Gene Genomes blast results', required=False)
	parser.add_argument('--bypass-parse', dest= 'bypass_parse' , default=True, action='store_true', help='Input file is parsed')
	parser.add_argument('--parse-only', dest= 'parse_only' , default=False, action='store_true', help='To only parse')
	parser.add_argument('--metagenome-fasta',dest='metagenome_fasta'  ,help='Metagenome FASTA file')
	parser.add_argument('--write-no-genomes',dest= 'write_no_genomes' ,default=False, action='store_true', help='Write Gene to No. of genomes')
	parser.add_argument('--bypass-hist', dest= 'bypass_hist', default=False, action='store_true', help='Generates Histogram')
	parser.add_argument('--bypass-scatter',dest= 'bypass_scatter' ,  default=True, action='store_true', help='Generates Scatterplot')
	parser.add_argument('--abund-prev', dest= 'abund_prev' , default=False, help='Centroid prevalence and abundance')
	parser.add_argument('--bypass-priority-scatter',dest= 'bypass_priority_scatter',  default=True, action='store_true', help='Generates Scatterplot')

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
	print"Total No. of genes:", len(table)
	uniq_genomes = []
	for gene in table:
		for genome in table[gene]:
			if genome not in uniq_genomes:
				uniq_genomes +=[genome]
	no_uniq_genomes = len(uniq_genomes)
	print 'No. of unique genomes: '+str(no_uniq_genomes)
	
	if not args.bypass_priority_scatter:
		abund_prev= read_abund_prev(args.abund_prev)
		master_plot(abund_prev, table, filename= args.input_file)#plot_metagenomic_priority(abund_prev, table, no_uniq_genomes, args.input_file)
	
	if args.parse_only:
		sys.exit('Input files parsed: '+args.input_file)
	
	if not args.bypass_scatter:
		plot_scatter(table, m8_filename, no_uniq_genomes)
	if not args.bypass_hist:
		plot_hist(table, m8_filename, no_uniq_genomes)
	
if __name__ == '__main__':
	main()
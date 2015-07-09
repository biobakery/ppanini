import os
import sys
import matplotlib
import re
import numpy
import time
import numpy
import argparse

from matplotlib import pyplot

numpy.seterr(divide='ignore', invalid='ignore')

def parse_table(table_name):
	'''Returns genes and uniref dicts from PPANINI output
	Input: 
	table_name =  filename of PPANINI output
	
	Output: 
	genes_table = {gene_id: {'alpha': 0.002, 'beta': 0.001, 'abund': 100.00} , ...}
	uniref_table = {uniref_id: {'alpha': 0.002, 'beta': 0.001, abund: 100.00}, ...} 
	inds = {'alpha': index, 'beta': index, abund: index}
	keys = ['alpha_xx', 'beta_xx', 'abund_xx']'''

	table_obj = open(table_name)
	table = table_obj.readlines()
	genes_table, uniref_table = {}, {}

	keys = [re.sub('[\r\t\n]','',i) for i in table[0].split('\t')[1:]]
	## only applicable for one alpha prevalence for NICHE <--should add multiple niche adaptability
	alpha_is = [i for i in range(len(keys)) if 'alpha' in keys[i]][0]
	beta_i = [i for i in range(len(keys)) if 'beta' in keys[i]][0]
	abund_i = [i for i in range(len(keys)) if 'abund' in keys[i]][0]
	
	for line in table[1:]:
		split_i = [re.sub('[\r\t\n]','',i) for i in line.split('\t')]
		if 'UniRef90' in split_i[0]: 
			uniref_table[split_i[0]] = [float(i) for i in split_i[1:]]
		else:
			genes_table[split_i[0]] = [float(i) for i in split_i[1:]]
	inds = {'alpha': alpha_is, 'beta': beta_i, 'abund': abund_i}
	return [genes_table, uniref_table, inds, keys]

def split_table(genes_table, inds):
	'''Splits gene_table into x, y coordinates of genes
	Input:
	genes_table = {gene_id : {'alpha': 0.001, 'beta': 0.001, 'abund': 100.00}, ...}
	inds = {'alpha': index, 'beta': index, abund: index}

	Output:
	dict = {'x': [list of x coordinates], 'y': [list of y coordinates]}'''

	x = inds[0]
	y = inds[1]

	x_i, y_i = [], []
	for gene in genes_table:
		x_i += [genes_table[gene][x]]
		y_i += [genes_table[gene][y]]
	return {'x': x_i, 'y': y_i}

def split_x_y(table, inds):
	'''Returns x and y coordinates for the uniref_ids and genes_ids separately
	Input:
	table = [genes_table = {gene_id : {'alpha': 0.001, 'beta': 0.001, 'abund': 100.00}, ...}, \
			 uniref_table = {uniref_id: {'alpha': 0.002, 'beta': 0.001, abund: 100.00}, ...}]
	inds = {'alpha': index, 'beta': index, abund: index}		 
	
	Output:
	genes_x_y =  {'x': [list of x coordinates], 'y': [list of y coordinates]}
	uniref_x_y = {'x': [list of x coordinates], 'y': [list of y coordinates]}'''

	[genes_table, uniref_table] = table

	genes_x_y = split_table(genes_table, inds)
	uniref_x_y = split_table(uniref_table, inds)
	return [genes_x_y, uniref_x_y]

def draw_cloud(cloud_points, data_points, labels, margins, zord):
	'''Saves a figure of Mean abundance vs. Prevalence for all the gene centroids 
	with important centroids highlighted

	Input:
	cloud_points = [c_genes_x_y = {'x': [list of x coordinates], 'y': [list of y coordinates]},\
					c_uniref_x_y = {'x': [list of x coordinates], 'y': [list of y coordinates]}]
	data_points = [d_genes_x_y = {'x': [list of x coordinates], 'y': [list of y coordinates]},\
				   d_uniref_x_y = {'x': [list of x coordinates], 'y': [list of y coordinates]},\
				   d_uniref_go_x_y = {'x': [list of x coordinates], 'y': [list of y coordinates]}]
	labels = {'xlabel': '', 'ylabel': '', 'title': '', 'filename': ''}
	margins = [abundance_threshold, prevalence_threshold]
	zord = Order of plotting the UniRef_GO, UniRef_NA, NA and Unprioritized'''

	pyplot.figure()
	pyplot.xlabel(labels['xlabel'])#'Alpha Prevalence_'+keys[i])
	pyplot.ylabel(labels['ylabel'])
	pyplot.title(labels['title']) #'Mean abundance')

	[c_genes_x_y, c_uniref_x_y] = cloud_points
	[d_genes_x_y, d_uniref_x_y, d_uniref_go_x_y] = data_points
	z = numpy.log(c_genes_x_y['y']+c_uniref_x_y['y'], dtype='float64')
	pyplot.scatter(numpy.log(c_genes_x_y['x'], dtype='float64'), \
				   numpy.log(c_genes_x_y['y'], dtype='float64'), \
				   c='grey', \
				   alpha=0.1, \
				   linewidths=0.0, \
				   zorder=zord[0], \
				   marker='o',\
				   label='Low Priority')
	pyplot.hold(True)
	pyplot.scatter(numpy.log(c_uniref_x_y['x'], dtype='float64'), \
				   numpy.log(c_uniref_x_y['y'], dtype='float64'), \
				   c='grey', \
				   alpha=0.1,\
				   linewidths =0.0, \
				   zorder=zord[0], \
				   marker='o',
				   label='Low Priority')
	##WithoutGO/UniRef90
	pyplot.scatter(numpy.log(d_genes_x_y['x'], dtype='float64'), \
				   numpy.log(d_genes_x_y['y'], dtype='float64'), \
				   c='yellow', \
				   alpha=0.5,\
				   linewidths =0.0, \
				   zorder=zord[3], \
				   marker= 'o', \
				   edgecolors='none',\
				   label='Prioritized unannotated')
	##WithUniRef90/WithoutGO
	pyplot.scatter(numpy.log(d_uniref_x_y['x'], dtype='float64'), \
				   numpy.log(d_uniref_x_y['y'], dtype='float64'), \
				   c='Red', \
				   alpha=0.5,\
				   linewidths =0.0, \
				   zorder=zord[1], \
				   marker='o', \
				   edgecolors='none',\
				   label='Prioritized UniRef90 annotated')
	##WITH UniRef AND GO
	pyplot.scatter(numpy.log(d_uniref_go_x_y['x'], dtype='float64'), \
	 			   numpy.log(d_uniref_go_x_y['y'], dtype='float64'), \
	 			   c='Blue', \
	 			   alpha=0.5,\
	 			   linewidths =0.0, \
	 			   zorder=zord[2],
	 			   edgecolors='none',\
	 			   label='Prioritized UniRef90 and GO annotated')

	pyplot.axhline(y=numpy.log(float(margins[0])), alpha=0.5, color='gray', label='abund>0.1+2SE')


	pyplot.axvline(x=numpy.log(float(margins[1])), alpha=0.5, color='gray', label='prev>0.1+2SE')

	pyplot.legend( loc=4, \
				   fontsize='x-small', \
				   framealpha=0.4, )	

	pyplot.savefig(labels['filename'])
	pyplot.savefig(labels['filename']+'.png')

def get_go_mapping(uniref_table, inds, map_go_fname, out_fname):
	'''Writes the uniref to GO mapping to file AND 
	Returns the x,y coordinates of UniRef_GO and UniRef_NA centroids

	Input:
	uniref_table = {uniref_id: {'alpha': 0.002, 'beta': 0.001, abund: 100.00}, ...}
	inds = {'alpha': index, 'beta': index, abund: index}		 
	map_go_fname = UniRef90 to GO mapping filename
	out_fname = output_filename

	Output:
	uniref_x_y = {'x': [list of x coordinates], 'y': [list of y coordinates]}
	uniref_go_x_y = {'x': [list of x coordinates], 'y': [list of y coordinates]}'''

	uniref, uniref_go = {}, {}
	map = open(map_go_fname,'r')
	mapper = {}
	for line in map:
		split_line = line.split('\t')
		for i in split_line[1:]:
			mapper[re.sub('[\r\t\n]','',i)] = split_line[0]
	with open(out_fname, 'w') as foo:
		for gene in uniref_table:
			if gene in mapper:
				uniref_go[gene] = uniref_table[gene]
				foo.writelines(gene+'\t'+mapper[gene]+'\n')
			else:
				foo.writelines(gene+'\tNA\n')
				uniref[gene] = uniref_table[gene]
	uniref_x_y = split_table(uniref, inds)
	uniref_go_x_y = split_table(uniref_go, inds)
	return [uniref_x_y, uniref_go_x_y]

def draw_prev_plot(genes_table, uniref_table, inds, labels, check):
	'''Draws and saves the sorted Prevalence or Abundance across gene centroids
	Input:
	genes_table = {gene_id: {'alpha': 0.002, 'beta': 0.001, abund: 100.00}, ...}
	uniref_table = {uniref_id: {'alpha': 0.002, 'beta': 0.001, abund: 100.00}, ...}
	inds = {'alpha': index, 'beta': index, abund: index}
	labels = {'xlabel': '', 'ylabel': '', 'title': '', 'filename': ''}
	check = DUBIOUS CHECK ??? '''
	alpha_is = inds
	prev = []
	
	for gene in genes_table:
		prev +=[genes_table[gene][alpha_is]]
	
	for gene in uniref_table:
		prev +=[uniref_table[gene][alpha_is]]
	arr_prev = numpy.array(prev)
	
	if check:
		prev_tmp = [i for i in prev if i>0]
	else:
		prev_tmp = prev
	
	tshld = 0.1+(2.0*float(numpy.std(numpy.array(prev_tmp)))/float(numpy.sqrt(len(prev_tmp))))
	
	arr_prev.sort()
	prev = arr_prev
	prev.sort()
	
	x= range(1, len(prev)+1)

	pyplot.figure()
	
	pyplot.scatter(x, numpy.log(prev), marker='.')
	pyplot.hold(True)
	pyplot.axhline(numpy.log(float(check)), alpha=0.5, color='red', label='>0.1+2SE')
	pyplot.legend(loc='upper left', fontsize='small')
	pyplot.xlabel(labels['xlabel'])
	pyplot.ylabel(labels['ylabel'])
	pyplot.title(labels['title'])
	pyplot.savefig(labels['filename'])
	pyplot.savefig(labels['filename']+'.png')

def draw_hexbin(genes_table, uniref_x_y, go_x_y, labels, inds):
	'''Draws and saves the HEXBIN plot for Prevalence and Abundance across gene centroids
	Input:
	genes_table = {gene_id: {'alpha': 0.002, 'beta': 0.001, abund: 100.00}, ...}
	uniref_x_y = {'x': [list of x coordinates], 'y': [list of y coordinates]}
	uniref_go_x_y = {'x': [list of x coordinates], 'y': [list of y coordinates]}

	inds = {'alpha': index, 'beta': index, abund: index}
	labels = {'xlabel': '', 'ylabel': '', 'title': '', 'filename': ''}'''

	pyplot.subplot(2,2,1)
	genes_x_y = split_table(genes_table, [inds['alpha'], inds['abund']])
	
	pyplot.hexbin(numpy.log(genes_x_y['x']), numpy.log(genes_x_y['y']), cmap='Blues', gridsize=10)
	pyplot.title(labels['title']+'_Genes')
	pyplot.xlabel(labels['xlabel'])
	pyplot.ylabel(labels['ylabel'])
	pyplot.colorbar()

	pyplot.subplot(2,2,2)
	
	pyplot.hexbin(numpy.log(uniref_x_y['x']), numpy.log(uniref_x_y['y']), cmap='Blues', gridsize=10)
	pyplot.title(labels['title']+'_UniRef')
	pyplot.xlabel(labels['xlabel'])
	pyplot.ylabel(labels['ylabel'])
	pyplot.colorbar()

	pyplot.subplot(2,2,3)

	pyplot.hexbin(numpy.log(go_x_y['x']), numpy.log(go_x_y['y']), cmap='Blues', gridsize=10)
	pyplot.title(labels['title']+'_GO')
	pyplot.xlabel(labels['xlabel'])
	pyplot.ylabel(labels['ylabel'])
	pyplot.colorbar()
	
	pyplot.subplot(2,2,4)
	pyplot.hold(True)
	pyplot.hexbin(numpy.log(go_x_y['x']), numpy.log(go_x_y['y']), cmap='Blues', gridsize=10, alpha=0.5)
	pyplot.hexbin(numpy.log(uniref_x_y['x']), numpy.log(uniref_x_y['y']), cmap='Reds', gridsize=10, alpha=0.5)
	pyplot.hexbin(numpy.log(genes_x_y['x']), numpy.log(genes_x_y['y']), cmap='Greys', gridsize=10, alpha=0.5)

	pyplot.savefig(labels['filename']+'.png')

def calculate_priority(array_x):
	'''Returns Metagenomic and Genomic priority for the gene lists on a scale of 0 to 1
	Input:
	array_x_y = {'x': [list of x coordinates], 'y': [list of y coordinates]}
	
	Output:
	gp = Genomic Priority [0.5,0.1, 0.3, ...]
	mp = Metagenomic Priority [0.9, 0.1, 0.2, ...]'''

	abund = numpy.array(array_x['y'])/max(array_x['y'])
	prev  = numpy.array(array_x['x'])/max(array_x['x'])

	# abund_order = [float(i) for i in numpy.argsort(abund)]
	# prev_order = [float(i) for i in numpy.argsort(prev)]

	# abund_order = numpy.array(abund_order)
	# prev_order = numpy.array(prev_order)

	gp = abund
	mp = [numpy.mean([abund[i], prev[i]]) for i in range(len(abund))]
	# mp = [(abund_order[i]*prev_order[i])/(abund_order[i]+prev_order[i]) for i in range(len(abund_order))]
	# gp = [float(i) for i in numpy.argsort(abund_order)]
	# mp = [float(i) for i in numpy.argsort(mp)]
	
	return [gp, mp]


def plot_priority(cloud_points, data_points, zord, labels):
	'''Plots Metagenome vs. Genome Priority plots
	
	cloud_points = [c_genes_x_y = {'x': [list of x coordinates], 'y': [list of y coordinates]},\
					c_uniref_x_y = {'x': [list of x coordinates], 'y': [list of y coordinates]}]
	data_points = [d_genes_x_y = {'x': [list of x coordinates], 'y': [list of y coordinates]},\
				   d_uniref_x_y = {'x': [list of x coordinates], 'y': [list of y coordinates]},\
				   d_uniref_go_x_y = {'x': [list of x coordinates], 'y': [list of y coordinates]}]
	labels = {'xlabel': '', 'ylabel': '', 'title': '', 'filename': ''}
	zord = Order of plotting the UniRef_GO, UniRef_NA, NA and Unprioritized'''

	[c_genes_x_y, c_uniref_x_y] = cloud_points
	[d_genes_x_y, d_uniref_x_y, d_uniref_go_x_y] = data_points

	master_x_y = {'x':[], 'y':[]}
	master_x_y['x'] += c_genes_x_y['x']
	master_x_y['y'] += c_genes_x_y['y']
	inds = [0, len(c_genes_x_y['x'])]
	
	master_x_y['x'] += c_uniref_x_y['x']
	master_x_y['y'] += c_uniref_x_y['y']
	inds += [len(inds), len(inds)+len(c_uniref_x_y['x'])]

	master_x_y['x'] += d_genes_x_y['x']
	master_x_y['y'] += d_genes_x_y['y']
	inds += [len(inds), len(inds)+len(d_genes_x_y['x'])]

	master_x_y['x'] += d_uniref_x_y['x']
	master_x_y['y'] += d_uniref_x_y['y']
	inds += [len(inds), len(inds)+len(d_uniref_x_y['x'])]

	master_x_y['x'] += d_uniref_go_x_y['x']
	master_x_y['y'] += d_uniref_go_x_y['y']
	inds += [len(inds), len(inds)+len(d_uniref_go_x_y['x'])]
	
	[gp, mp] = calculate_priority(master_x_y)


	pyplot.figure()
	pyplot.xlabel(labels['xlabel'])
	pyplot.ylabel(labels['ylabel'])
	pyplot.title(labels['title'])

	pyplot.scatter(numpy.log(gp[inds[0]:inds[1]], dtype='float64'), \
				   numpy.log(mp[inds[0]:inds[1]], dtype='float64'), \
				   c='grey', \
				   alpha=0.1, \
				   linewidths=0.0, \
				   zorder=zord[0], \
				   marker='o',\
				   label='Low Priority')
	pyplot.hold(True)

	pyplot.scatter(numpy.log(gp[inds[2]:inds[3]], dtype='float64'), \
				   numpy.log(mp[inds[2]:inds[3]], dtype='float64'), \
				   c='grey', \
				   alpha=0.1,\
				   linewidths =0.0, \
				   zorder=zord[0], \
				   marker='o',
				   label='Low Priority')
	##WithoutGO/UniRef90
	pyplot.scatter(numpy.log(gp[inds[4]:inds[5]], dtype='float64'), \
				   numpy.log(mp[inds[4]:inds[5]], dtype='float64'), \
				   c='yellow', \
				   alpha=0.5,\
				   linewidths =0.0, \
				   zorder=zord[3], \
				   marker= 'o', \
				   edgecolors='none',\
				   label='Prioritized unannotated')
	##WithUniRef90/WithoutGO
	pyplot.scatter(numpy.log(gp[inds[6]:inds[7]], dtype='float64'), \
				   numpy.log(mp[inds[6]:inds[7]], dtype='float64'), \
				   c='Red', \
				   alpha=0.5,\
				   linewidths =0.0, \
				   zorder=zord[1], \
				   marker='o', \
				   edgecolors='none',\
				   label='Prioritized UniRef90 annotated')
	##WITH UniRef AND GO
	pyplot.scatter(numpy.log(gp[inds[8]:inds[9]], dtype='float64'), \
	 			   numpy.log(mp[inds[8]:inds[9]], dtype='float64'), \
	 			   c='Blue', \
	 			   alpha=0.5,\
	 			   linewidths =0.0, \
	 			   zorder=zord[2],
	 			   edgecolors='none',\
	 			   label='Prioritized UniRef90 and GO annotated')
	
	pyplot.legend( loc=4, \
				   fontsize='x-small', \
				   framealpha=0.4, )	

	pyplot.savefig(labels['filename'])
	pyplot.savefig(labels['filename']+'.png')

if __name__ == '__main__':
	##UniRef=red; Unannotated; blue
	## argv[1]= original map; argv[2]=selected
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--input_table', help='Gene abundance table with metadata', required=True)
	parser.add_argument('--original_table', default=False, help='Gene abundance table with metadata')
	parser.add_argument('--bypass_cloud', action='store_true', default=False, help='To draw Abundance Prevalence Cloud')
	parser.add_argument('--prev', default=False, help='Graph will be prevalence across centroids')
	parser.add_argument('--abund',default=False, help='Graph will be mean abundance across centroids')
	parser.add_argument('-m','--mapper', help='GO to UniRef mapper')
	parser.add_argument('--write_mapper', action='store_true', default=False, help='Gene to GO table written')
	parser.add_argument('--zorder', default='1,2,3,4', help='Zorder [1,2,3,4] [Old, UniRef, UniRef/GO, NA]')
	parser.add_argument('--hexplot', default=False, action='store_true', help='Plot HEXBIN')
	parser.add_argument('--bypass_priority', default=False, action='store_true', help='Generates Metagenome vs. Genome Priority plots')

	args = parser.parse_args()
	zorder =  [int(i) for i in args.zorder.split(',')]
	[genes_table_i, uniref_table_i, inds_i, keys_i] = parse_table(args.input_table)
	alpha_is = inds_i['alpha']
	abund_i = inds_i['abund']
	
	map_go_fname = args.mapper
	[uniref_i_x_y, uniref_go_x_y] = get_go_mapping(uniref_table_i, [alpha_is, abund_i], map_go_fname, args.input_table+'_GO_map.txt_tmp')
		
	if args.write_mapper:
		r = open(args.input_table+'_GO_map.txt_tmp')
		with open(args.input_table+'_GO_map.txt','w') as foo:
			for line in r:
				foo.writelines(line)
			for gene in genes_table_i:
				foo.writelines(gene+'\tNA\n')
		os.system('rm '+args.input_table+'_GO_map.txt_tmp')
	abund=  float(args.abund)
	prev = float(args.prev)+0.1
	
	if args.original_table:
		[genes_table, uniref_table, inds, keys] = parse_table(args.original_table)
		
		genes_table_o, uniref_table_o = {}, {}
		for gene in genes_table:
			if gene not in genes_table_i:
				genes_table_o[gene] = genes_table[gene]

		for gene in uniref_table:
			if gene not in uniref_table_i:
				uniref_table_o[gene] = uniref_table[gene]
		
		cloud_points = split_x_y([genes_table_o, uniref_table_o], \
								 [alpha_is, abund_i])
		data_points = split_x_y([genes_table_i, uniref_table_i], \
								[alpha_is, abund_i])
		labels = {'xlabel': keys_i[alpha_is], \
				  'ylabel': keys_i[abund_i], \
				  'title': 'Mean Abundance vs. '+keys_i[alpha_is], \
				  'filename': args.input_table+'_cloud.pdf'}
		data_points[1] = uniref_i_x_y
		data_points += [uniref_go_x_y]	

		if not args.bypass_cloud:
			draw_cloud(cloud_points, data_points, labels, [abund, prev], zorder)
		if not args.bypass_priority:
			labels = {'xlabel': 'Genomic priority', \
				  'ylabel': 'Metagenomic priority', \
				  'title': 'Metagenomic vs. Genomic Priority', \
				  'filename': args.input_table+'_mpgp.pdf'}
			plot_priority(cloud_points, data_points, zorder, labels)

	if args.prev:
		name = args.input_table
		title=name.split('/')[-1].split('.')[0]
		labels = {'xlabel': 'Gene Centroids', \
				  'ylabel': 'Alpha Prevalence', \
				  'title': title, \
				  'filename': name+'_prev.pdf'}
		draw_prev_plot(genes_table_i, uniref_table_i, alpha_is, labels, prev)

	if args.abund:
		name = args.input_table
		title=name.split('/')[-1].split('.')[0]
		labels = {'xlabel': 'Gene Centroids', \
				  'ylabel': 'Mean Abundance', \
				  'title': title, \
				  'filename': name+'_abund.pdf'}
		draw_prev_plot(genes_table_i, uniref_table_i, abund_i, labels, abund)

	if args.hexplot:
		name = args.input_table
		title=name.split('/')[-1].split('.')[0]
		labels = {'xlabel': 'Gene Prevalence', \
				  'ylabel': 'Mean Abundance', \
				  'title': title, \
				  'filename': name+'_hexbin.pdf'}
		draw_hexbin(genes_table_i, uniref_i_x_y, uniref_go_x_y, labels, inds_i)


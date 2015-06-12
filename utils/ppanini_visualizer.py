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

numpy.seterr(divide='ignore', invalid='ignore')

def parse_table(table_name):
	table_obj = open(table_name)
	table = table_obj.readlines()
	genes_table, uniref_table = {}, {}

	keys = [re.sub('[\r\t\n]','',i) for i in table[0].split('\t')[1:]]
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
	x = inds[0]
	y = inds[1]

	x_i, y_i = [], []
	for gene in genes_table:
		x_i += [genes_table[gene][x]]
		y_i += [genes_table[gene][y]]
	return {'x': x_i, 'y': y_i}

def split_x_y(table, inds):
	[genes_table, uniref_table] = table

	genes_x_y = split_table(genes_table, inds)
	uniref_x_y = split_table(uniref_table, inds)
	return [genes_x_y, uniref_x_y]

def draw_cloud(cloud_points, data_points, labels, margins, zord):
	pyplot.figure()
	pyplot.xlabel(labels['xlabel'])#'Alpha Prevalence_'+keys[i])
	pyplot.ylabel(labels['ylabel'])
	pyplot.title(labels['title']) #'Mean abundance')

	[c_genes_x_y, c_uniref_x_y] = cloud_points
	[d_genes_x_y, d_uniref_x_y, d_uniref_go_x_y] = data_points
	z= numpy.log(c_genes_x_y['y']+c_uniref_x_y['y'], dtype='float64')
	z1= numpy.percentile(z, 90)
	z2= numpy.percentile(z, 75)
	z3= numpy.percentile(z, 50)
	z4= numpy.percentile(z, 30)
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
	#pdb.set_trace()
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
	uniref, uniref_go = {}, {}
	map = open(map_go_fname,'r')
	mapper = {}
	t= time.time()
	for line in map:
		split_line = line.split('\t')
		for i in split_line[1:]:
			mapper[re.sub('[\r\t\n]','',i)] = split_line[0]
	t2=time.time()
	print (t2-t)/60
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
	#pdb.set_trace()
	return [uniref_x_y, uniref_go_x_y]

def draw_prev_plot(genes_table, uniref_table, inds, labels, check):
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
	# print per_90
	arr_prev.sort()
	prev = arr_prev
	prev.sort()
	#pdb.set_trace()
	pyplot.figure()
	x= range(1, len(prev)+1)
	#pyplot.ylim([min(prev), max(prev)])
	pyplot.scatter(x, numpy.log(prev), marker='.')
	pyplot.hold(True)
	pyplot.axhline(numpy.log(float(check)), alpha=0.5, color='red', label='>0.1+2SE')
	pyplot.legend(loc='upper left', fontsize='small')
	pyplot.xlabel(labels['xlabel'])
	pyplot.ylabel(labels['ylabel'])
	pyplot.title(labels['title'])
	pyplot.savefig(labels['filename'])
	pyplot.savefig(labels['filename']+'.png')


if __name__ == '__main__':
	##UniRef=red; Unannotated; blue
	## argv[1]= original map; argv[2]=selected
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--input_table', help='Gene abundance table with metadata', required=True)
	parser.add_argument('--original_table', help='Gene abundance table with metadata')
	parser.add_argument('--abund_prev', action='store_true', default=False, help='To draw Abundance Prevalence Cloud')
	parser.add_argument('--prev', default=False, help='Graph will be prevalence across centroids')
	parser.add_argument('--abund',default=False, help='Graph will be mean abundance across centroids')
	parser.add_argument('-m','--mapper', help='GO to UniRef mapper')
	parser.add_argument('--write_mapper', action='store_true', default=False, help='Gene to GO table written')
	parser.add_argument('--zorder', default='1,2,3,4', help='Zorder [1,2,3,4] [Old, UniRef, UniRef/GO, NA]')

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
	abund=  float(args.abund) #+200.0
	prev = float(args.prev)+0.1
	
	if args.abund_prev:	
		[genes_table, uniref_table, inds, keys] = parse_table(args.original_table)
		
		genes_table_o, uniref_table_o = {}, {}
		for gene in genes_table:
			if gene not in genes_table_i:
				genes_table_o[gene] = genes_table[gene]

		for gene in uniref_table:
			if gene not in uniref_table_i:
				uniref_table_o[gene] = uniref_table[gene]
		
		## Abundance vs. Alpha Prevalence ##
		
		cloud_points = split_x_y([genes_table_o, uniref_table_o], \
								 [alpha_is, abund_i])
		data_points = split_x_y([genes_table_i, uniref_table_i], \
								[alpha_is, abund_i])
		print time.time()
		t=time.time()
		labels = {'xlabel': keys_i[alpha_is], \
				  'ylabel': keys_i[abund_i], \
				  'title': 'Mean Abundance vs. '+keys_i[alpha_is], \
				  'filename': args.input_table+'_cloud.pdf'}
		# [uniref_i_x_y, uniref_go_x_y] = get_go_mapping(uniref_table_i, [alpha_is, abund_i], map_go_fname, args.input_table+'_GO_map.txt')
		print (time.time()-t)/60
		data_points[1] = uniref_i_x_y
		data_points += [uniref_go_x_y]	

		draw_cloud(cloud_points, data_points, labels, [abund, prev], zorder)

	if args.prev:
		name = args.input_table
		title=name.split('/')[-1].split('.')[0]
		labels = {'xlabel': 'Gene Centroids', \
				  'ylabel': 'Alpha Prevalence', \
				  'title': title, \
				  'filename': name+'_prev.pdf'}
		draw_prev_plot(genes_table, uniref_table, inds_i['alpha'], labels, prev)

	if bool(args.abund):
		name = args.input_table
		title=name.split('/')[-1].split('.')[0]
		labels = {'xlabel': 'Gene Centroids', \
				  'ylabel': 'Mean Abundance', \
				  'title': title, \
				  'filename': name+'_abund.pdf'}
		draw_prev_plot(genes_table, uniref_table, inds_i['abund'], labels, abund)

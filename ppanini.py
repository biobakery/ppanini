import os
import sys
import re
import argparse
import numpy
import logging
import scipy.stats
sys.path.append('/n/hutlab12_nobackup/data/ppanini/ppanini')#'/Users/rah/Documents/Hutlab/ppanini')#
import ppanini
from src import utilities
from src import annotate_genes
from src import config

logger = logging.getLogger(__name__)


beta = 0.5

numpy.seterr(divide='ignore', invalid='ignore')

def read_gene_table(gene_table_fname):
	'''Returns the different elements from the gene table

	Input: gene_table_fname = Filename of the gene_table

	Output: metadata = [metadata strings]; Rows with # as first character in table
			uniref_dm = {UniRef_XYZ: numpy.array(abundance), ...}
			gi_dm = {GENE_ID_UNKNOWN_UNIREF: numpy.array(abundance), ...}'''
	
	logger.debug('read_gene_table')

	gene_table = open(gene_table_fname) 
	metadata = []
	uniref_dm, gis_dm = {}, {}
	count = 0
	flag =False
	for line in gene_table:
		if count == 350000:
			break
		count +=1
		if line.startswith('#'):
			metadata += [line]
		else:
			split_i = line.split('\t')
			annot = split_i[0].split('|') #geneID column split 
			#print annot
			u90_annot = [i for i in annot if 'UniRef90' in i][0]
			u50_annot = [i for i in annot if 'UniRef50' in i][0]
			#print 'u50_annot: ',u50_annot
			data_row = numpy.array([float(i) for i in split_i[1:]])
			#print 'data_row: ', data_row
			if flag or u90_annot in config. essantial_genes_uniref90_id:
				flag = True
				print 'Essential!'
				#count +=1
				if 'UniRef90_unknown' == u90_annot:
					#if 'UniRef50_unknown' == u50_annot:
					try: #same name
						gis_dm[annot[0]] += data_row
					except:
						gis_dm[annot[0]] = data_row
					'''
					else: #same uniref90 id
						try:
							uniref_dm[u50_annot] += data_row
						except KeyError:
							uniref_dm[u50_annot] = data_row
					'''
				else: #same uniref90 id
					try:
						uniref_dm[u90_annot] += data_row
					except KeyError:
						#print "Data row",data_row
						uniref_dm[u90_annot] = data_row	
				
	if config.verbose == 'DEBUG':
		print ("Gene Table contains %s genes." % count)
	return [uniref_dm, gis_dm, metadata]

def get_centroids_fromUCLUST(genes):
	'''Returns the clusters dictionary

	Input: gene_centroid_clusters_file_path = Filename of the centroids UC file
		   genes = [List of genes that are unannotated with UniRef]

	Output: cluster_dict = {CENTROIDS: [list of genes], ...}'''
	gene_centroid_clusters_file_path = config.uclust_file
	logger.debug('get_centroids_fromUCLUST: '+gene_centroid_clusters_file_path)

	cluster_dict = {}
	genes_clustered = []
	foo = open(gene_centroid_clusters_file_path)

	for line in foo:
		if line.startswith('H'):
			split_line = [re.sub('[\r\t\n]','', i) for i in line.split('\t')[-2:]]
			if sum([1 for i in split_line if i in genes]):
				try:
					cluster_dict[split_line[1]] += [split_line[0]]
					genes_clustered += [split_line[0]]
				except:
					cluster_dict[split_line[1]] = [split_line[0], split_line[1]]
					genes_clustered += [split_line[0], split_line[1]]

	for i in genes:
		if not i in genes_clustered:
			cluster_dict[i] = [i]

	return cluster_dict

def get_centroids(uniref_dm, gi_dm):
	'''Returns the dict of all centroids containing clusters of gene IDs

	Input:	uniref_dm = {UniRef_XYZ: numpy.array(abundance), ...}
			gi_dm = {GENE_ID_UNKNOWN_UNIREF: numpy.array(abundance), ...}
			usearch_folder = Location of the USEARCH program or in path if not provided
			uclust_file = path to USEARCH_UCLUST FILE from PREPPANINI
			gene_catalog = path to gene catalog containing all genes
			nprocesses = Number of threads to run for clustering

	Output: all_centroids = {gene_centroid : numpy.array(abundance), ...}'''
	
	logger.debug('get_centroids')

	centroids_fasta = {}

	if config.uclust_file == '':
		centroid_gis = get_clusters() #all UniRef90_unknowns are clustered across samples
	else:
		centroid_gis = get_centroids_fromUCLUST(gi_dm.keys())

	gc_dm = {}
	for centroid in centroid_gis:
		for gene in centroid_gis[centroid]:
			if gene in gi_dm:
				try:
					gc_dm[centroid] += gi_dm[gene]
				except:
					gc_dm[centroid] = gi_dm[gene]
	
	for centroid in uniref_dm:
		gc_dm[centroid] = uniref_dm[centroid]
	
	return gc_dm


def get_clusters(): #ONLY FOR THE UNIREF UNANNOTATED
	'''Returns the dict of unannotated gene centroids containing clusters of genes at 90% similarity

	Input:	gene_catalog = path to all genes catalog}
			args = args object containing information about the clust program
			nprocesses = Number of threads to run for clustering

	Output: centroid_gis = {gene_centroid: [List of genes in the cluster]}'''

	logger.debug('get_clusters')

	allgenes_file_path = config.gene_catalog
	gene_centroids_file_path = config.temp_folder+'/'+config.basename+'_centroids.fasta'
	gene_centroid_clusters_file_path = config.temp_folder+'/'+config.basename+'_clusters.uc'
	
	clust_method = 'vsearch'

	if config.usearch:
		clust_method = config.usearch
		annotate_genes.run_uclust(clust_method, allgenes_file_path, gene_centroids_file_path, gene_centroid_clusters_file_path, 0.9)
	else:
		if config.vsearch:
			clust_method = config.usearch
		annotate_genes.run_vclust(clust_method, allgenes_file_path, gene_centroids_file_path, gene_centroid_clusters_file_path, 0.9)
		
	
	centroid_gis = annotate_genes.get_clusters_dict(gene_centroid_clusters_file_path)

	return centroid_gis


def get_centroids_table(all_centroids, metadata):
	'''Returns data matrix containing gene centroids and abundance per sample

	Input:	metadata = [metadata strings]; Rows with # as first character in table
			all_centroids = all_centroids = {gene_centroid : numpy.array(abundance), ...}

	Output: norm_data_matrix = numpy.array([0,0,0.0,...],[0,0,0.0,...],...)
			centroids_list = [List of centroids]'''
	
	logger.debug('get_centroids_table')

	centroids_data_matrix = []
	centroids_list = []
	
	for centroid in all_centroids:
		centroids_list +=[centroid]
		centroids_data_matrix +=[all_centroids[centroid]]

	#NORMALIZATION PER SAMPLE
	centroids_data_matrix = numpy.array(centroids_data_matrix)
	norm_data_matrix = centroids_data_matrix/sum(centroids_data_matrix)
	norm_data_matrix = norm_data_matrix*1e6

	gene_centroids_table_file_path = config.temp_folder+'/'+config.basename+'_gene_centroids_table.txt'
	
	with open(gene_centroids_table_file_path,'w') as foo:
		foo.writelines(metadata)
		for i, val in enumerate(centroids_list):
			foo.writelines([str.join('\t', [val] + [str(j) for j in norm_data_matrix[i]]) + '\n'])
			

	return [norm_data_matrix, centroids_list]


def get_prevalence_abundance(centroids_data_matrix, centroids_list, metadata, beta):
	'''Returns the dict of centroids with their prevalence and abundance

	Input:	centroids_data_matrix = {gene_centroid: [Gene centroid abundance across samples]}
			metadata = [metadata strings]; Rows with # as first character in table

	Output: centroid_prev_abund = {centroid: {'mean_abundance': mean abundance, 'prevalence': prevalence}}
			all_prevalence = [List of all observed gene centroid prevalence values (>0) across samples]
			all_abund = [List of all calculated gene centroid abundance across samples]
			flag = True (if NICHE PRESENT) or False(if NICHE ABSENT)'''
	
	logger.debug('get_prevalence_abundance')

	centroid_prev_abund_file_path = config.temp_folder+'/'+config.basename+'_centroid_prev_abund.txt'
	
	[niche_line, ind] = utilities.is_present(metadata, '#NICHE')
	#print niche_line
	if niche_line:
		niche_flag = True
		centroid_prev_abund = get_niche_prevalence_abundance (centroids_data_matrix, centroids_list, niche_line, beta)
		
		#return centroid_prev_abund#[centroid_prev_abund, all_alpha_prev, all_mean_abund, niche_flag]
	else:
		niche_flag = False
		centroid_prev_abund = {}
		all_prevalence = [] 
		all_abund = []
		#print centroids_data_matrix
		for iter, centroid in enumerate(centroids_data_matrix):
			'''
			befor was
			for centroid in centroids_data_matrix:
			'''
			#abund only where the gene is present in sample
			#print 'cenetroid_id', centroids_list[iter], 'cenroid', centroid
			abund_i = [i for i in  centroid if i > 0]
			'''
			before was  
			abund_i = [i for i in  centroids_data_matrix[centroid] if i > 0]
			'''
			abund_score = numpy.mean(abund_i)
			prev_score = float(sum(centroid > 0)/\
						 float(len(centroid)))
			''' before was 
			prev_score = float(sum(numpy.array(centroids_data_matrix[centroid]) > 0)/\
						 float(len(centroids_data_matrix[centroid])))
						 '''

			centroid_prev_abund[centroids_list[iter]] = {'mean_abundance': abund_score, \
											 'prevalence': prev_score}
			
			all_prevalence += [centroid_prev_abund[centroids_list[iter]]['prevalence']]
			all_abund += abund_i
			#config.all_prevalence = all_prevalence
			#config.all_mean_abund = all_abund
		
		all_prevalence = sorted(all_prevalence)
		all_abund = sorted(all_abund)
		
		for centroid in centroid_prev_abund:
			p_score = 1/((1/(beta*(scipy.stats.percentileofscore(all_prevalence, centroid_prev_abund[centroid]['prevalence']))))+\
				     (1/((1-beta)*(scipy.stats.percentileofscore(all_abund, centroid_prev_abund[centroid]['mean_abundance'])))))
			centroid_prev_abund[centroid]['ppanini_score'] = p_score
		
	write_prev_abund_matrix(centroid_prev_abund, centroid_prev_abund_file_path)
	
	config.niche_flag = niche_flag
	return centroid_prev_abund #[centroid_prev_abund, all_prevalence, all_abund, niche_flag]

def get_niche_prevalence_abundance(centroids_data_matrix, centroids_list, niche_line, beta):
	'''Returns the dict of centroids with their prevalence and abundance

	Input:	centroids_data_matrix = {gene_centroid: [Gene centroid abundance across samples]}
			metadata = [metadata strings]; Rows with # as first character in table
			line = line from mapper file containing NICHE information

	Output: centroid_prev_abund = {centroid: {'mean_abundance': float_mean abundance, 
											  'alpha_prevalence': {'NICHE**': prevalence of centroid across samples within niche, ...}},
											  'b_prev': float_median of alpha prevalence within niches for centroid}}
			all_alpha_prev = {niche: List of all observed gene centroid alpha prevalence values (>0) across samples within each niche}
			all_mean_abund = [List of all calculated mean gene centroid abundance across samples]'''

	logger.debug('get_niche_prevalence_abundance')
	
	centroid_prev_abund_file_path = config.temp_folder+'/'+config.basename+'_centroid_prev_abund.txt'
	
	niches = {}
	split_i = [re.sub('[\r\t\n]', '', i) for i in niche_line.split('\t')[1:]]
	for i, val in enumerate(split_i):
		try:
			niches[val] += [i]
		except KeyError:
			niches[val] = [i]
	
	niches_label = niches.keys()
	centroid_prev_abund = {}

	all_alpha_prev = {}
	all_alpha_abund = []

	for niche in niches_label:
		all_alpha_prev[niche] = []

	for i, centroid in enumerate(centroids_list):
		centroid_prev_abund[centroid] = {}
		a_prev = {}
		a_abund = {}
		for niche in niches_label:
			a_prev[niche] = 0
			a_abund[niche] = []
 			for ind in niches[niche]:
 				check_i = int(centroids_data_matrix[i][ind] > 0) #present in a sample
				a_prev[niche] += check_i
				if check_i:
					a_abund[niche] += [centroids_data_matrix[i][ind]]

			a_prev[niche] = float(a_prev[niche])/float(len(niches[niche]))
			all_alpha_prev[niche] += [a_prev[niche]]

			if not a_abund[niche]: #if gene completely empty throughout samples
				a_abund[niche] = [0]
			a_abund[niche] = numpy.mean(numpy.array(a_abund[niche]))
		
		max_mean_abund = max(a_abund.values())
		all_alpha_abund += [max_mean_abund]
		b_prev = numpy.median(a_prev.values())

		centroid_prev_abund[centroid]['alpha_prevalence'] = a_prev
		centroid_prev_abund[centroid]['mean_abundance'] = max_mean_abund #a_abund
		centroid_prev_abund[centroid]['beta_prevalence'] = b_prev
	
	#Percentile of score requires sorted vectors! Blekh!
	all_alpha_abund = sorted(all_alpha_abund)
	for niche in all_alpha_prev:
		all_alpha_prev[niche] = sorted(all_alpha_prev[niche])

	for centroid in centroid_prev_abund:
		p_score = {}
		for niche in centroid_prev_abund[centroid]['alpha_prevalence']:
			p_score[niche] = 1/((1/((beta)*scipy.stats.percentileofscore(all_alpha_prev[niche], centroid_prev_abund[centroid]['alpha_prevalence'][niche])))+\
															  (1/((1-beta)*scipy.stats.percentileofscore(all_alpha_abund, centroid_prev_abund[centroid]['mean_abundance']))))
			centroid_prev_abund[centroid]['ppanini_score'] = p_score

	dict_to_print = {}
	for centroid in centroid_prev_abund:
		dict_to_print[centroid] = {'beta_prevalence': centroid_prev_abund[centroid]['beta_prevalence'], \
								   'mean_abundance': centroid_prev_abund[centroid]['mean_abundance']}
		for niche in centroid_prev_abund[centroid]['alpha_prevalence']:
			dict_to_print[centroid]['alpha_prevalence_'+niche] = centroid_prev_abund[centroid]['alpha_prevalence'][niche]
			dict_to_print[centroid]['ppanini_score_'+niche] = centroid_prev_abund[centroid]['ppanini_score'][niche]
	#config.all_prevalence = all_prevalence
	#config.all_mean_abund = all_alpha_abund
	write_prev_abund_matrix(dict_to_print, centroid_prev_abund_file_path)
	
	return centroid_prev_abund

def get_important_niche_centroids():
	'''Returns the dict of important gene centroids [>= 10th percentile of alpha_prevalence and mean abundance]

	Input:	centroid_prev_abund = {centroid: {'mean_abundance': mean abundance, 'prevalence': prevalence}}
			all_alpha_prev = {niche: List of all observed gene centroid alpha prevalence values (>0) across samples within each niche}
			all_mean_abund = [List of all calculated mean gene centroid abundance across samples]
			output_folder = Location of the results folder

	Output: imp_centroids = {centroid: {'mean_abundance': mean abundance, 
										'beta_prevalence': median of alpha prevalences observed in each niche, 
										'alpha_prevalence_NICHEX': alpha_prevalence for niche X, ...}}'''

	logger.debug('get_important_niche_centroids')
	centroid_prev_abund = config.centroid_prev_abund
	beta = config.beta
	tshld = config.tshld
	output_folder = config.output_folder
	imp_centroid_prev_abund_file_path = config.basename+'_imp_centroid_prev_abund.txt'
	
	
	tshld_abund = config.tshld_abund
	tshld_prev = config.tshld_prev
	
	ppanini_score = 1/((1/(beta*tshld_prev)) + (1/((1-beta)*tshld_abund)))

	logger.debug('get_important_niche_centroids: tshld_prev:'+str(config.tshld_prev))
	logger.debug('get_important_niche_centroids: ppanini_score:'+str(ppanini_score))
	logger.debug('get_important_niche_centroids: tshld_abund:'+str(tshld_abund))
	
	imp_centroids = {}

	for centroid in centroid_prev_abund:
		check = sum([centroid_prev_abund[centroid]['ppanini_score'][niche]>=ppanini_score for niche in centroid_prev_abund[centroid]['ppanini_score']])
		if check:
			imp_centroids[centroid]={'mean_abundance': centroid_prev_abund[centroid]['mean_abundance'], \
									 'beta_prevalence': centroid_prev_abund[centroid]['beta_prevalence']}
			for niche in centroid_prev_abund[centroid]['alpha_prevalence']:
				imp_centroids[centroid]['alpha_prevalence_' + niche] = centroid_prev_abund[centroid]['alpha_prevalence'][niche]
				imp_centroids[centroid]['ppanini_score_'+niche] = centroid_prev_abund[centroid]['ppanini_score'][niche]
	
	write_prev_abund_matrix(imp_centroids, output_folder + '/' + imp_centroid_prev_abund_file_path)
	
	return imp_centroids


def get_important_centroids():
	'''Returns the dict of important gene centroids [value-2SE(prevalence and abundance) >0.1]

	Input:	centroid_prev_abund = {centroid: {'mean_abundance': mean abundance, 'prevalence': prevalence}}
			all_prevalence = [List of all observed gene centroid prevalence values (>0) across samples]
			all_mean_abund = [List of all calculated mean gene centroid abundance across samples]
			output_folder = Location of the results folder

	Output: imp_centroids = {centroid: {'mean_abundance': mean abundance, 'prevalence': prevalence}}'''
	
	logger.debug('get_important_centroids')
	beta = config.beta
	centroid_prev_abund = config.centroid_prev_abund
	imp_centroid_prev_abund_file_path = config.basename+'_imp_centroid_prev_abund.txt'
	
	tshld_prev = config.tshld_prev
	tshld_abund = config.tshld_abund
	ppanini_score =  1/((1/(beta*tshld_prev)) + (1/((1-beta)*tshld_abund)))
	
	imp_centroids = {}

	for centroid in centroid_prev_abund:
		check = centroid_prev_abund[centroid]['ppanini_score'] >= ppanini_score
		if check:
			imp_centroids[centroid]={'mean_abundance': centroid_prev_abund[centroid]['mean_abundance'], \
									 'prevalence': centroid_prev_abund[centroid]['prevalence'],\
									 'ppanini_score': centroid_prev_abund[centroid]['ppanini_score']}
	
	write_prev_abund_matrix(imp_centroids, config.output_folder + '/' + imp_centroid_prev_abund_file_path)
	
	return imp_centroids

def write_prev_abund_matrix(centroid_prev_abund, out_file):
	'''Writes the centroids prevalence and abundance information in text file

	Input: centroid_prev_abund = {centroids: {'mean_abundance': mean abundance, ...}}
		   out_file = output_filename

	Output: Writes the centroids dictionary to the output_filename'''

	logger.debug('write_prev_abund_matrix')

	keys = []
	for i in centroid_prev_abund:
		keys = centroid_prev_abund[i].keys()
		break
	
	with open(out_file,'w') as foo:
		foo.writelines(['#Centroids\t' + str.join('\t', keys) + '\n'])
		for centroid in centroid_prev_abund:
			foo.writelines([str.join('\t', [centroid] + [str(centroid_prev_abund[centroid][key]) for key in keys]) + '\n'])

def read_prevalence_abundance_table(input_table):
	foo = open(input_table)
	abund_i = 0
	beta_i = 0
	alphas_i = 0
	baseline = 0
	niche_flag = 0
	centroid_prev_abund = {}
	all_prevalence = []
	all_mean_abund = []
	for line in foo:
		if baseline:
			split_line = [re.sub('[\r\t\n]','',i) for i in line.split('\t')]
			centroid_prev_abund[split_line[0]] = {'mean_abundance': float(split_line[abund_i]),
											  	  'beta_prevalence': float(split_line[beta_i])}
			if niche_flag:
				for i in alphas_i:
					niche_i = i[1].split('_')[-1]
					try:
						centroid_prev_abund[split_line[0]]['alpha_prevalence'][niche_i] =float(split_line[i[0]])
					except:
						centroid_prev_abund[split_line[0]]['alpha_prevalence']= {niche_i:float(split_line[i[0]])}
					all_prevalence[niche_i] += [centroid_prev_abund[split_line[0]]['alpha_prevalence'][niche_i]]
			else:
				centroid_prev_abund[split_line[0]][alphas_i[1]] =float(split_line[alphas_i[0]])
				all_prevalence += [centroid_prev_abund[split_line[0]][alphas_i[1]]]
			all_mean_abund += [centroid_prev_abund[split_line[0]]['mean_abundance']]
		else:
			split_line = line.split('\t')			
			split_line = [re.sub('[\r\t\n]','',i) for i in line.split('\t')]

			abund_i = split_line.index('mean_abundance')
			beta_i = split_line.index('beta_prevalence')
			print split_line
		
			alphas_i = [(i,val) for i, val in enumerate(split_line) if 'alpha_prevalence' in val]
			if sum([len(i[1].split('_')) for i in alphas_i])>1:
				niche_flag =1
				all_prevalence = {}
				for i in alphas_i:
					niche_i = i[1].split('_')[-1]
					all_prevalence[niche_i] = []
			else:
				alphas_i = alphas_i[0]
			baseline +=1
	return [centroid_prev_abund, all_prevalence, all_mean_abund, niche_flag]
def read_parameters():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i','--input_table', help='REQUIRED: Gene abundance table with metadata', required=True)
	parser.add_argument('-o','--output-folder', dest = 'output_folder',  help='Folder containing results', default=False)
	parser.add_argument('--gene-catalog', dest = 'gene_catalog', default=config.gene_catalog, help='GENE CATALOG')
	parser.add_argument('--uc', default= config.uclust_file, help='UCLUST file containg centroids and clustered genes')
	parser.add_argument('--usearch', default = config.usearch, help='Path to USEARCH') #add to be in path?
	parser.add_argument('--vsearch', default = config.vsearch, help='Path to VSEARCH') #add to be in path?
	parser.add_argument('--basename', default = config.basename, help='BASENAME for all the output files')
	parser.add_argument('--log-level', dest = 'log_level',  default=config.log_level, help='Choices: [DEBUG, INFO, WARNING, ERROR, CRITICAL]')
	parser.add_argument('--threads', default= config.nprocesses, type=int,help='Number of threads')
	parser.add_argument('--tshld-abund', dest = 'tshld_abund', default=config.tshld_abund, type = float,help='[X] Percentile Cutoff for Abundance; Default=75th')
	parser.add_argument('--tshld-prev', dest = 'tshld_prev', default=config.tshld_prev, type =float, help='Percentile cutoff for Prevalence')
	parser.add_argument('--beta', default=config.beta, help='Beta parameter for weights on percentiles')
	# parser.add_argument('--bypass-prev-abund', dest = 'bypass_prev_abund', default=False, action='store_true', help='Bypass quantifying abundance and prevalence')

	args = parser.parse_args()
	config.nprocesses = args.threads
	config.basename = args.basename
	config.input_table = args.input_table
	config.uclust_file = args.uc
	config.gene_catalog = args.gene_catalog
	
	config.beta = args.beta

def run():
	if config.uclust_file == '' and config.gene_catalog == '':
		sys.exit("At least one of --uc or --gene-catalog should be provided!!!")
	# if not args.bypass_prev_abund:
	if config.basename =='':
		config.basename = config.input_table.split('.')[0].split('/')[-1]
	if config.output_folder == '':
		config.output_folder = config.basename
	config.temp_folder = config.output_folder+'/'+config.basename+'_temp'

	utilities.create_folders([config.temp_folder, config.output_folder])

	log_file = config.output_folder+'/'+config.basename+'.log'
	logging.basicConfig(filename=log_file, \
						format='%(asctime)s - %(name)s - %(levelname)s: %(message)s', \
						level=getattr(logging, config.log_level), \
						filemode='w', \
						datefmt='%m/%d/%Y %I:%M:%S %p')
	
	if config.verbose =='DEBUG':
		print "Reading the gene table..."
	[uniref_dm, gi_dm, metadata]= read_gene_table(config.input_table)
	if config.verbose =='DEBUG':
		print "Reading the gene table is done!"
	
	if config.verbose =='DEBUG':
		print "Getting centroids..."
	all_centroids = get_centroids(uniref_dm, gi_dm)
	if config.verbose =='DEBUG':
		print "Getting centroids is done!"
	
	if config.verbose =='DEBUG':
		print "Getting centroids table..."
	[centroids_data_matrix, centroids_list] = get_centroids_table(all_centroids, metadata)
	config.centroids_list = centroids_list
	if config.verbose =='DEBUG':
		print "Getting centroids table is done!"
	
	if config.verbose =='DEBUG':
		print "Getting prevelance abundnace..."
	centroid_prev_abund = get_prevalence_abundance(centroids_data_matrix, centroids_list = centroids_list, metadata = metadata, beta = config.beta)
	if config.verbose =='DEBUG':
		print "Getting prevelance abundnace is done!"
	# else:
	# 	[centroid_prev_abund, all_prevalence, all_mean_abund, niche_flag] = read_prevalence_abundance_table(input_table, config.beta)
	config.centroid_prev_abund = centroid_prev_abund
def  prioritize_centroids():
	if config.verbose =='DEBUG':
		print "Prioritize centroids..."
	if config.niche_flag:
		imp_centroids = get_important_niche_centroids()
	else:
		imp_centroids = get_important_centroids()
	if config.verbose =='DEBUG':
		print 'Prioritize centroids is done!'
	return imp_centroids

if __name__ == '__main__':
	read_parameters()
	run()
	prioritize_centroids()
	
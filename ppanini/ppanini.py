import os
import sys
import subprocess
import platform
import re
import pdb
import argparse
import numpy
import logging
import scipy.stats
import ppanini

from os.path import basename
from numpy import percentile
from collections import namedtuple
from operator import attrgetter, itemgetter
ppanini_table_row = namedtuple("ppanini_table_row", ["alpha_prevalence", "prevalence_percentile", "mean_abundance","abund_percentile", "beta_prevalence", "ppanini_score", "GO"], verbose=False, rename=False)
import pandas as pd
import pkg_resources

try:
    from . import utilities
    from . import config

except ImportError as e:
    sys.exit("CRITICAL ERROR: Unable to find the PPANINI python package." +
        " Please check your install."+str(e))

logger = logging.getLogger(__name__)

numpy.seterr(divide='ignore', invalid='ignore')


def read_gene_table(config=config):
    '''Returns the different elements from the gene table
    
    config.input_table = Filename of the gene_table
    
    Output: metadata = [metadata strings]; Rows with # as first character in table
    		uniref_dm = {UniRef_XYZ: numpy.array(abundance), ...}
    		gi_dm = {GENE_ID_UNKNOWN_UNIREF: numpy.array(abundance), ...}'''
    
    logger.debug('read_gene_table')
    #gene_table = pd.read_csv(config.input_table, sep='\t', index_col=0)

    gene_table = utilities.gzip_bzip2_biom_open_readlines(config.input_table)
    metadata = []
    uniref_dm, gis_dm = {}, {}
    count_metadata_lines = 0
    count_gene_lines = 0
    for line in gene_table:
    	if line.startswith('#') or line.startswith('\t') :
            metadata += [line]
            count_metadata_lines +=1
    	else:
            count_gene_lines +=1
            
            split_i = line.split('\t')
            annot = split_i[0].split('|') #geneID column split 
            try:
            	u90_annot = [i for i in annot if 'UniRef90' in i][0]
            except: #Incase Gene table is not annotated with UniRef90
            	u90_annot = 'UniRef90_unknown'
            try:
            	u50_annot = [i for i in annot if 'UniRef50' in i][0]
            except: #Incase Gene table is not annotated with UniRef50
            	u50_annot = 'UniRef50_unknown'
            
            data_row = numpy.array([float(i) for i in split_i[1:]])
    		
            if 'UniRef90_unknown' == u90_annot:
            	if 'UniRef50_unknown' == u50_annot:
            		try: #same name
            			gis_dm[annot[0]] += data_row
            		except:
            			gis_dm[annot[0]] = data_row
            	else: #same uniref90 id
            		try:
            			uniref_dm[u50_annot] += data_row
            		except KeyError:
            			uniref_dm[u50_annot] = data_row
            else: #same uniref90 id
            	try:
            		uniref_dm[u90_annot] += data_row
            	except KeyError:
            		uniref_dm[u90_annot] = data_row	
	
    if config.verbose == 'DEBUG':
    	print ("--- Gene Table contains %s metadata lines." % count_metadata_lines)
        print ("--- Gene Table contains %s gene families." % count_gene_lines)
    
    return [uniref_dm, gis_dm, metadata]

def summerize_centroids(uniref_dm, gi_dm, config=config):
    '''Returns the dict of all centroids containing clusters of gene IDs
    
    Input:	uniref_dm = {UniRef_XYZ: numpy.array(abundance), ...}
    		gi_dm = {GENE_ID_UNKNOWN_UNIREF: numpy.array(abundance), ...}
    		
    Output: gc_dm = {gene_centroid : numpy.array(abundance), ...}'''
    
    logger.debug('summerize_centroids')   
    centroid_gis = {}
    for gene in gi_dm:
        centroid_gis[gene] = [gene]
        
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
    print ("--- Number of gene families: %d"%(len(gc_dm)))
    return gc_dm


def normalize_centroids_table(all_centroids, metadata):
    '''Returns data matrix containing gene centroids and abundance per sample
    
    Input:	metadata = [metadata strings]; Rows with # as first character in table
    		all_centroids = all_centroids = {gene_centroid : numpy.array(abundance), ...}
    
    Output: norm_data_matrix = numpy.array([0,0,0.0,...],[0,0,0.0,...],...)
    		centroids_list = [List of centroids]'''
    
    logger.debug('get_centroids_table')
    #print all_centroids
    centroids_data_matrix = pd.DataFrame.from_dict(all_centroids).T
    samples = metadata[0].split('\t')
    centroids_data_matrix.columns = samples[1:]
    # remove centroids with all zero abundances  
    #centroids_data_matrix = centroids_data_matrix[(centroids_data_matrix.T != 0).any()]  
    
    # write abundance table
    abundance_table_path = config.temp_folder+'/'+config.basename+'_abundance_table.txt'
    centroids_data_matrix.to_csv(abundance_table_path, sep='\t', mode='w', header=True)
    
    # Normalize abundance per sample
    norm_data_matrix = centroids_data_matrix.apply(lambda x: x/x.sum(), axis=0)
    #centroids_data_matrix.div(centroids_data_matrix.sum(axis=1), axis=1)#centroids_data_matrix/sum(centroids_data_matrix)
    #norm_data_matrix = norm_data_matrix*1e6
    gene_centroids_table_file_path = config.temp_folder+'/'+config.basename+'_abundance_table_norm.txt'
    
    # write metadat
    #with open(gene_centroids_table_file_path,'w') as foo:
    #	foo.writelines(metadata)
        
    # write normalized abundance
    norm_data_matrix.to_csv(gene_centroids_table_file_path, sep='\t', mode='w', header=True)
    return norm_data_matrix

def get_prevalence_abundance(centroids_data_matrix, metadata):
    '''Returns the dict of centroids with their prevalence and abundance
    
    Input:	centroids_data_matrix = {gene_centroid: [Gene centroid abundance across samples]}
    		centroids_list = list of all the centroids
    		metadata = [metadata strings]; Rows with # as first character in table
    		beta = parameter value 
    
    Output: summary_table = {centroid: {'mean_abundance': mean abundance, 'prevalence': prevalence}}
    		all_prevalence = [List of all observed gene centroid prevalence values (>0) across samples]
    		all_abund = [List of all calculated gene centroid abundance across samples]
    		flag = True (if NICHE PRESENT) or False(if NICHE ABSENT)'''
    
    logger.debug('get_prevalence_abundance')
    
    set_niches = []
    [niche_line, ind] = utilities.is_present(metadata, '#NICHE')
    if len(niche_line) >0:
        set_niches = list(set(re.split(r'\t|\n|\r+', niche_line)))
        set_niches = list([item.title() for item in set_niches ])
    # use this if niche is specified as a row in the abundance table
    # and we have more than one niche (there is #Niche + '' as part of niche line + one niche 
        set_niches.remove('')
        set_niches.remove('#Niche')
        config.niches = set_niches
    if niche_line and len(set_niches) > 1  :
        print ("--- Niches have been provided in abundance table are:", set_niches)
        niche_flag = True
        summary_table = get_niche_prevalence_abundance (centroids_data_matrix, niche_line)
    
    # use this if niche is NOT specified as a row in the abundance table
    else:
    	niche_flag = False
        summary_table = get_no_niche_prevalence_abundance(centroids_data_matrix)
    config.niche_flag = niche_flag
    return summary_table 

def get_no_niche_prevalence_abundance(centroids_data_matrix):
    '''Returns the dataframe of centroids with their prevalence, abundance, percentiles, and PPANINI score
    
    Input:    centroids_data_matrix = {gene_centroid: [Gene centroid abundance across samples]}
            centroids_list = list of all the centroids
            metadata = [metadata strings]; Rows with # as first character in table
            beta = parameter value 
    
    Output: summary_table = {centroid: {'mean_abundance': mean abundance, 'prevalence': prevalence}}'''
            
    
    logger.debug('get_no_niche_prevalence_abundance')

    # cerate a data frame for one niche data
    summary_table = pd.DataFrame(index=centroids_data_matrix.index, columns=ppanini_table_row._fields)
       
    # Calculate mean abundance of non zero values
    df = centroids_data_matrix.replace(0, numpy.NaN)
    summary_table['mean_abundance'] = df.mean(axis = 1)

        # Alpha prevalence 
    summary_table['alpha_prevalence'] = df.count(axis = 1)/df.shape[1]
    
    # Calculate percentile for prevalence and abundance
    summary_table['prevalence_percentile'] = scipy.stats.rankdata(summary_table['alpha_prevalence'], method='average')/summary_table.shape[0]
    summary_table['abund_percentile'] = scipy.stats.rankdata(summary_table['mean_abundance'], method='average')/summary_table.shape[0]
    
    # Calculate PPANINI score
    summary_table['ppanini_score'] = 1/((config.beta/summary_table['prevalence_percentile'])+((1-config.beta)/summary_table['abund_percentile']))   
    
    return summary_table 

def get_niche_prevalence_abundance(centroids_data_matrix, niche_line):
    
    logger.debug('get_niche_prevalence_abundance')
    
    
    # cerate a data frame for niches 
    abundance_cols = []
    config.ppanini_niche_score_labels = []
    ppanini_columns = list(ppanini_table_row._fields)
    
    for niche in config.niches:
        ppanini_columns.append('alpha_prevalence_'+niche)
        ppanini_columns.append('prevalence_percentile_'+niche)
        abundance_cols.append('alpha_abundance_'+niche)
        ppanini_columns.append('alpha_abundance_'+niche)
        ppanini_columns.append('abundance_percentile_'+niche)
        
    for niche in config.niches:
        ppanini_columns.append('ppanini_score_'+niche)
        config.ppanini_niche_score_labels.append('ppanini_score_'+niche)
        
    summary_table = pd.DataFrame(index=centroids_data_matrix.index, columns=ppanini_columns)
    # Calculate mean abundance of non zero values
    #print centroids_data_matrix
    df = centroids_data_matrix.replace(0, numpy.NaN)
    summary_table['mean_abundance'] = df.mean(axis = 0)
   
    # Alpha prevalence 
    summary_table['alpha_prevalence'] = centroids_data_matrix.astype(bool).sum(axis=0)/centroids_data_matrix.shape[1]
    
    # get the index of columns for each niche
    niches = {}
    split_i = [re.sub('[\r\t\n]', '', i) for i in niche_line.split('\t')[1:]]
    for i, val in enumerate(split_i):
    	try:
    		niches[val] += [i]
    	except KeyError:
    		niches[val] = [i]
    
    niches_label = niches.keys()
        
    # Calculate alpha prevalence and abundance for each niche
    summary_table['mean_abundance'] = 0
    for niche in config.niches:
        #print niches[niche], niche, centroids_data_matrix
        summary_table['alpha_prevalence_'+niche] = centroids_data_matrix.loc[:,niches[niche]].astype(bool).sum(axis=1)/len(niches[niche])
        summary_table['prevalence_percentile_'+niche] = scipy.stats.rankdata(summary_table['alpha_prevalence_'+niche], method='average')/len(niches[niche]) 
        
        summary_table['alpha_abundance_'+niche] = df.loc[:, niches[niche]].mean(axis=1)
        summary_table['abundance_percentile_'+niche] = scipy.stats.rankdata(summary_table['alpha_prevalence_'+niche], method='average')/len(niches[niche]) 

    # update mean abundance
    summary_table.loc['mean_abundance'] = summary_table.loc[:,abundance_cols].max(axis=1) 
    
    summary_table['abund_percentile'] = scipy.stats.rankdata(summary_table['mean_abundance'], method='average')/summary_table.shape[0]

    
    # Calculate PPANINI score of reach niche
    
    for niche in config.niches:
        summary_table['ppanini_score_'+niche] = 1/((config.beta/summary_table['prevalence_percentile_'+niche])+((1-config.beta)/summary_table['abund_percentile']))   
        
    return summary_table 

def impotrance_measure(scores,q=.1, n = 10):
    max_rank = len(scores)
    random_prev = numpy.random.randint(max_rank, size=n)
    random_abund = numpy.random.randint(max_rank, size=n)
    print random_prev, random_abund
    random_score = [(.5 * prev + .5 * abund) for prev in random_prev for abund in random_abund]
    pvalue = [None]*max_rank
    for i in range(max_rank):
       pvalue[i] = (sum([1.0 if item <= scores[i] else 0.0 for item in random_score])+1)/n 
    return pvalue

def get_important_centroids(config=config):
    
    logger.debug('get_important_centroids')
    beta = config.beta
    summary_table = config.summary_table
    imp_summary_table_file_path = config.basename+'_ppanini_table.txt'
    
    tshld_prev = config.tshld_prev
    tshld_abund = config.tshld_abund
    ppanini_score = -1# 1/((1/(beta*tshld_prev)) + (1/((1-beta)*tshld_abund)))
    #summary_table["pvalue"] = impotrance_measure(summary_table['ppanini_score'])
    #Get important centroids based on their PPANINI score
    if config.niche_flag:
        imp_centroids = summary_table[summary_table.loc[:,config.ppanini_niche_score_labels].max(axis=1) >= ppanini_score]
        try:
            # for  pandas >= 0.17.0
            imp_centroids = imp_centroids.sort_values(by=config.ppanini_niche_score_labels, ascending=False)
        except:
            imp_centroids = imp_centroids.sort(ppanini_niche_score_labels, ascending=False)
    else:
        imp_centroids = summary_table[summary_table['ppanini_score'] >= ppanini_score]
        try:
            # for  pandas >= 0.17.0
            imp_centroids = imp_centroids.sort_values(by='ppanini_score', ascending=False)
        except:
            imp_centroids = imp_centroids.sort('ppanini_score', ascending=False)
        
    #imp_centroids.to_csv( config.temp_dir + '/' + imp_summary_table_file_path, sep='\t')
    return imp_centroids

def write_prev_abund_matrix(summary_table, out_file):
    '''Writes the centroids prevalence and abundance information in text file
    
    Input: summary_table = {centroids: {'mean_abundance': mean abundance, ...}}
    	   out_file = output_filename
    
    Output: Writes the centroids dictionary to the output_filename'''
    
    logger.debug('write_prev_abund_matrix')
    keys = []
    for i in summary_table:
    	keys = summary_table[i]._asdict().keys()
    	break
    with open(out_file,'w') as foo:
    	foo.writelines(['#gene_families\t' + str.join('\t', list(keys)) + '\n'])
    	for centroid in summary_table:
    		foo.writelines([str.join('\t', [centroid] + [str(summary_table[centroid][key]) for key in range(len(keys))]) + '\n'])

def read_prevalence_abundance_table(input_table):
	'''Need to redo this'''
	foo = open(input_table)
	abund_i = 0
	beta_i = 0
	alphas_i = 0
	baseline = 0
	niche_flag = 0
	summary_table = {}
	all_prevalence = []
	all_mean_abund = []
	for line in foo:
		if baseline:
			split_line = [re.sub('[\r\t\n]','',i) for i in line.split('\t')]
			summary_table[split_line[0]] = {'mean_abundance': float(split_line[abund_i]),
											  	  'beta_prevalence': float(split_line[beta_i])}
			if niche_flag:
				for i in alphas_i:
					niche_i = i[1].split('_')[-1]
					try:
						summary_table[split_line[0]]['alpha_prevalence'][niche_i] =float(split_line[i[0]])
					except:
						summary_table[split_line[0]]['alpha_prevalence']= {niche_i:float(split_line[i[0]])}
					all_prevalence[niche_i] += [summary_table[split_line[0]]['alpha_prevalence'][niche_i]]
			else:
				summary_table[split_line[0]][alphas_i[1]] =float(split_line[alphas_i[0]])
				all_prevalence += [summary_table[split_line[0]][alphas_i[1]]]
			all_mean_abund += [summary_table[split_line[0]]['mean_abundance']]
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
	return [summary_table, all_prevalence, all_mean_abund, niche_flag]

def read_parameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input_table', help='REQUIRED: Gene abundance table with metadata', required=True)
    parser.add_argument('-o','--output-folder', dest = 'output_folder',  help='Folder containing results', required=False, default=config.temp_dir)
    #parser.add_argument('--gene-catalog', dest = 'gene_catalog', default=config.gene_catalog, help='GENE CATALOG')
    #parser.add_argument('--uc', default= config.uclust_file, help='UCLUST file containg centroids and clustered genes')
    #parser.add_argument('--usearch', action="store_true", default = config.usearch, help='Path to USEARCH') #add to be in path?
    #parser.add_argument('--vsearch', action="store_true", default = config.vsearch, help='Path to VSEARCH') #add to be in path?
    parser.add_argument('--basename', default = config.basename, help='BASENAME for all the output files')
    parser.add_argument('--uniref2go', default = config.uniref2go, help='uniref to GO term mapping file')
    parser.add_argument('--log-level', dest = 'log_level',  default=config.log_level, help='Choices: [DEBUG, INFO, WARNING, ERROR, CRITICAL]')
    #parser.add_argument('--threads', default= config.nprocesses, type=int,help='Number of threads')
    parser.add_argument('--tshld-abund', dest = 'tshld_abund', default=config.tshld_abund, type = float,help='[X] Percentile Cutoff for Abundance; Default=75th')
    parser.add_argument('--tshld-prev', dest = 'tshld_prev', default=config.tshld_prev, type =float, help='Percentile cutoff for Prevalence')
    parser.add_argument('--beta', default=config.beta, help='Beta parameter for weights on percentiles')
    #parser.add_argument('--bypass-clustering', dest = 'bypass_clustering', default=config.bypass_clustering, action='store_true', help='Bypass clustering')
    parser.add_argument('--version', action="version", version="%(prog)s "+config.version, help='prints the version')
    # parser.add_argument('--bypass-prev-abund', dest = 'bypass_prev_abund', default=False, action='store_true', help='Bypass quantifying abundance and prevalence')
    
    args = parser.parse_args()
    config.beta = args.beta 
    config.temp_dir = args.output_folder
    config.basename = args.basename
    config.input_table = args.input_table
    config.bypass_clustering = True
    config.uniref2go = args.uniref2go
    
def run():
 
    # check for sofwares here
    
    #
    
    if config.basename=='':
        #config.basename = basename(config.input_table).split('.')[0]
        #print config.basename
        pass
    if config.temp_dir == '':
    	config.temp_dir = config.basename
    
    config.temp_folder = config.temp_dir+'/temp'

    utilities.create_folders([config.temp_dir, config.temp_folder])
    
    log_file = config.temp_dir+'/'+config.basename+'.log'
    logging.basicConfig(filename=log_file, \
    					format='%(asctime)s - %(name)s - %(levelname)s: %(message)s', \
    					level=getattr(logging, config.log_level), \
    					filemode='w', \
    					datefmt='%m/%d/%Y %I:%M:%S %p')
    
    if config.verbose =='DEBUG':
    	print "--- Reading the gene table..."
    [uniref_dm, gi_dm, metadata]= read_gene_table()
    #print uniref_dm, gi_dm, metadata
    
    if config.verbose =='DEBUG':
    	print "--- Summarize gene families table ..."
    all_centroids = summerize_centroids(uniref_dm, gi_dm)
    
    
    if config.verbose =='DEBUG':
    	print "--- Normalize gene families table ..."
    centroids_data_table = normalize_centroids_table(all_centroids, metadata)
    #config.centroids_list = centroids_list
    
   
    if config.verbose =='DEBUG':
    	print "--- Getting prevalence abundance ..."
    summary_table = get_prevalence_abundance(centroids_data_table, \
    												metadata = metadata)
        
    if config.genomic_score: # an option should be added for this
        metagenomic_table  = utilities.read_parsed("/Users/rah/Documents/PPANINI/ppanini_old_files/PARSED_BLAST_RESULTS/AN_mg.m8")
        metagenomic_table.to_csv(config.temp_dir + '/' +config.basename+'_metagenomic_table.txt', sep='\t')
    # else:
    # 	[summary_table, all_prevalence, all_mean_abund, niche_flag] = read_prevalence_abundance_table(input_table, config.beta)
    config.summary_table = summary_table
   
    # add Go terms to the table
    if not config.uniref2go == '':
        if config.verbose =='DEBUG':
            print("--- Mapping UniRef90 to GO terms!")
        utilities.uniref2go(config.summary_table, uniref_go_path = config.uniref2go)
    else:
        if config.verbose =='DEBUG':
            print("--- Mapping UniRef90 to GO terms!")
        resource_package = __name__  # Could be any module/package name
        resource_path = '/'.join(('data', 'map_uniref90_infogo1000.txt.gz'))
        template = pkg_resources.resource_filename(resource_package, resource_path)
        utilities.uniref2go(config.summary_table, uniref_go_path = template)   
    
def _main():
    read_parameters()
    run()
    
    if config.verbose =='DEBUG':
        print "--- Prioritize gene families ..."
    imp_centroids = get_important_centroids()
    
    # set the path for ppanini table out
    ppanini_output_file_path = config.temp_dir+'/'+config.basename+'_table.txt' 
    # Write PPANINI output table
    imp_centroids.to_csv(ppanini_output_file_path, sep='\t')
    if config.verbose =='DEBUG':
        print ("--- The PPANINI output is written in %s ..." % (config.output_folder))
        print "--- PPANINI process is successfully completed ..."

if __name__ == '__main__':
	_main()

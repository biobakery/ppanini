import os
import sys
import matplotlib
import re
import time
import numpy as np
import argparse
import pdb
import pandas as pd

import matplotlib
matplotlib.use( "Agg" )
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["font.family"] = "Arial"
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import matplotlib.cm as cm

from matplotlib import font_manager
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("error")
    try:
        font_file = font_manager.findfont(font_manager.FontProperties(family='Arial'))
        matplotlib.rcParams["font.family"] = "Arial"
    except UserWarning:
        pass 


from .. import utilities


def plot_scatter(table, txt_filename, no_uniq_genomes):
	'''Plots Scatter plot for genome hits per gene
	Input:
	txt_filename = filename of blast results
	table = {gene: [List of genomes]}
	no_uniq_genomes = Number of Unique Genomes in Metagenomic niche'''

	labels = {'xlabel': 'Prioritized Centroids',\
			  'ylabel':'No. of Genomes (Log10)', \
			  'title':'Metagenome vs. Genome Prioritization',\
			  'filename': txt_filename+'_prioritizationScatter.pdf'}	

	all_genes = []
	for gene in table:
		all_genes +=[len(table[gene])/float(no_uniq_genomes)]

	all_genes.sort()
	plt.figure()
	plt.xlabel(labels['xlabel'])
	plt.ylabel(labels['ylabel'])
	plt.title(labels['title']) 
	plt.scatter(np.arange(len(all_genes)), \
				   np.log(all_genes), \
				   c='grey', \
				   alpha=0.1, \
				   linewidths=0.1, \
				   zorder=1, \
				   marker='o',\
				   label='All Centroids')
	'''plt.legend( loc=4, \
				   fontsize='x-small', \
				   framealpha=0.4, )
				   '''	
	plt.savefig(labels['filename'])
	plt.savefig(labels['filename']+'.png')

def plot_hexbin(table, txt_filename):
	'''Plots HexBin plots for the genome hits per gene
	Input:
	txt_filename = filename of blast results
	table = {gene: [List of genomes]}'''

	labels = {'xlabel': 'Prioritized Centroids',\
			  'ylabel':'No. of Genomes', \
			  'title':'Metagenome vs. Genome Prioritization',\
			  'filename': txt_filename+'_prioritizationHEXBIN.pdf'}
	plt.figure()
	plt.xlabel(labels['xlabel'])
	plt.ylabel(labels['ylabel'])
	plt.title(labels['title']) 

	all_genes = []
	for gene in table:
		all_genes +=[len(table[gene])/4189.0]
	all_genes.sort()
	plt.hold(True)
	image = plt.hexbin(np.arange(len(all_genes)), \
				   all_genes,\
				   bins='log', gridsize=30, mincnt=1, cmap= 'darkgoldenrod', zorder=1) 
	cb = plt.colorbar(image, spacing='uniform', extend='max')
	
	plt.savefig(labels['filename'])

def plot_hist(table, txt_filename, no_uniq_genomes):
	'''Plots histogram for the genome hits per gene
	Input:
	txt_filename = filename of blast results
	table = {gene: [List of genomes]}
	no_uniq_genomes = Number of Unique Genomes in Metagenomic niche'''

	labels = {'ylabel': 'Centroids',\
			  'xlabel':'Genomes', \
			  'title':'Metagenome vs. Genome Prioritization',\
			  'filename': txt_filename+'_prioritizationHIST.pdf'}
	
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
	
	
	hist, bins = np.histogram(all_genes, bins=200)
	offset = bins[1:]-bins[:-1]
	axarr.plot(bins[:-1]+offset, np.cumsum(hist), color='darkgoldenrod')
	plt.savefig(labels['filename']+'_cumsum.pdf')
	plt.savefig(labels['filename']+'_cumsum.png')


def master_plot(path, size = 5):
	data_scale = 'log'
	fig, axarr = plt.subplots(nrows=2, ncols=2, dpi=300)#, sharex=True, sharey=True)
	fig.set_size_inches(size, size)
		
	# Read stool samples information
	# Stool
	metagenomic_table_stool =  utilities.gene2genomes(path +'stool.txt')
	ppanini_output_stool = utilities.pd.DataFrame.from_csv(path +'/stool_ppanini_table.txt',
											 sep='\t', index_col=0, header =0)
	no_uniq_genomes_stool  = utilities.number_of_unique_genomes(path +'/stool.txt')
	scatter_plot_metagenomic_priority(axarr[0, 0], ppanini_output_stool, metagenomic_table_stool, no_uniq_genomes_stool, title = 'Stool', scale = data_scale)
	
	# Buccal mucosa
	metagenomic_table_mucosa =  utilities.gene2genomes(path +'/mucosa.txt')
	ppanini_output_mucosa = utilities.pd.DataFrame.from_csv(path +'/mucosa_ppanini_table.txt',
											 sep='\t', index_col=0, header =0)
	no_uniq_genomes_mucosa  = utilities.number_of_unique_genomes(path +'mucosa.txt')
	scatter_plot_metagenomic_priority(axarr[0, 1], ppanini_output_mucosa, metagenomic_table_mucosa, no_uniq_genomes_mucosa, title = 'Buccal mucosa', scale = data_scale)
	
	# Posterior fornix 
	metagenomic_table_fornix =  utilities.gene2genomes(path +'fornix.txt')
	ppanini_output_fornix = utilities.pd.DataFrame.from_csv(path +'/fornix_ppanini_table.txt',
											 sep='\t', index_col=0, header =0)
	no_uniq_genomes_fornix  = utilities.number_of_unique_genomes(path +'fornix.txt')
	scatter_plot_metagenomic_priority(axarr[1, 0], ppanini_output_fornix, metagenomic_table_fornix, no_uniq_genomes_fornix, title = 'Posterior fornix', scale = data_scale)
	
	
	# Anterior nares
	metagenomic_table_nares =  utilities.gene2genomes(path +'nares.txt')
	ppanini_output_nares = utilities.pd.DataFrame.from_csv(path +'/nares_ppanini_table.txt',
											 sep='\t', index_col=0, header =0)
	no_uniq_genomes_nares  = utilities.number_of_unique_genomes(path +'/nares.txt')
	scatter_plot_metagenomic_priority(axarr[1, 1], ppanini_output_nares, metagenomic_table_nares, no_uniq_genomes_nares, title = 'Anterior nares', scale = data_scale)
	
	# Invisible y axis labels for second column
	plt.setp([a.set_ylabel('') for a in axarr[:, 1]])
	
	# Invisible x axis labels for first row
	plt.setp([a.set_xlabel('') for a in axarr[0, :]])
	
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
	#st = fig.suptitle('Metagenomic vs. Genomic Priority', fontsize=12, fontweight='bold', va='top')
	
	# shift subplots up:
	#st.set_y(0.95)
	#fig.subplots_adjust(top=.9)
	fig.subplots_adjust()
	
	
	plt.savefig(path+'/fig4b'+'.pdf', bbox_inches='tight', pad_inches = 0, dpi=300) 
	plt.savefig(path+'/fig4b'+'.png', bbox_inches='tight', pad_inches = 0, dpi=300) 
	plt.savefig(path+'/fig4b'+'.svgz', bbox_inches='tight', pad_inches = 0, dpi=300)
	#fig, axarr = plt.subplots(nrows=2, ncols=2, dpi=300)#, sharex=False, sharey=False)
	#fig.set_size_inches(10, 10)
	#hexbin_plot_metagenomic_priority(axarr[0, 0], ppanini_output1, metagenomic_table1, no_uniq_genomes1, title = 'Stool', scale = data_scale)
	
	#hexbin_plot_metagenomic_priority(axarr[1, 1], ppanini_output2, metagenomic_table2, no_uniq_genomes2, title = 'Anterior nares', scale = data_scale)
	
	#hexbin_plot_metagenomic_priority(axarr[1, 0], ppanini_output3, metagenomic_table3, no_uniq_genomes3, title = 'Buccal mucosa', scale = data_scale)
	
	#im= hexbin_plot_metagenomic_priority(axarr[0, 1], ppanini_output4, metagenomic_table4, no_uniq_genomes4, title = 'Posterior fornix', scale = data_scale)
	
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
	#st = fig.suptitle('Metagenomic vs. Genomic Priority', fontsize=12, fontweight='bold', va='top')
	
	# shift subplots up:
	'''st.set_y(0.95)
	fig.subplots_adjust(top=.9)
	fig.subplots_adjust()
	
	
	plt.savefig(path + '/fig4b_s'+'sclae_'+str(data_scale)+'.pdf', bbox_inches='tight', pad_inches = 0, dpi=300) 
'''
def scatter_plot_metagenomic_priority(axe, ppanini_table, title = None, score_type = 'universe', pan_genome_score = None, scale = None, output_path = None, size = 3, draw_contour = True):
    if axe == None:
    	fig, axe = plt.subplots(1, figsize=(size, size))
    #mp_gp = {}
    genes = ppanini_table.index
    #abund = ppanini_table['abundance']
    #prev = ppanini_table['prevalence']
    mp = ppanini_table['ppanini_score']
    #abund = np.array(abund)/max(abund)
    #prev = np.array(prev)/max(prev)
    draw_contour = draw_contour
    if score_type == 'universe':
        
        if  pan_genome_score:
            gp_pangenome = [float(pan_genome_score[gene]) if gene in pan_genome_score else float(0.0) for gene in genes]
            #gp = [x + y for x, y in zip(gp_niche, gp_pangenome)]
            gp = gp_pangenome
        else:
            sys.exit("Please provide pangenome score dictionary")
    elif score_type == 'niche':
        gp_niche = np.zeros(len(metagenomic_table))
        uniq_genomes = []
        gp_niche = [len(metagenomic_table[gene])/float(n_uniq_genomes) for gene in genes]
        gp_niche = np.array(gp_niche)
        gp = gp_niche
    temp_gp=[]
    temp_mp=[]
     
    for i in range(len(gp)):
    	if np.isnan(gp[i]) or np.isnan(mp[i]) or gp[i]==0 or mp[i]==0 :
    		continue
    	temp_gp += [gp[i]]
    	temp_mp += [mp[i]]
    
    gp = np.array(temp_gp) 
    mp = np.array(temp_mp)
    if scale == 'log':
    	gp = np.log(gp)
    	#mp = np.log(mp)
    elif scale == 'sqrt':
    	gp = np.sqrt(gp)
    	#mp = np.sqrt(mp) 		
    #fig, axe = plt.subplots()
    def ncolors( n, colormap="jet" ):
    	"""utility for defining N evenly spaced colors across a color map"""
    	cmap = plt.get_cmap( colormap )
    	cmap_max = cmap.N
    	return [cmap( int( k * cmap_max / (n - 1) ) ) for k in range( n )]
    my_color = ncolors(4)
    #x_dic= {'Gut': '#b87333', 'Skin':'#ffff00', 'Oral':'#009fff', 'Vaginal':'#ff4d00'}
    color_dic= {'Stool': my_color[0], 'Anterior nares':my_color[2], 'Buccal mucosa':my_color[1], 'Posterior fornix':my_color[3]}
    
    # Calculate the point density
    x = np.array(gp) 
    y = np.array(mp)
    if draw_contour:
        xy = np.vstack([x,y])
        z = gaussian_kde(xy)(xy)
        
        # Sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        #print idx, x, y
        x, y, z = x[idx], y[idx], z[idx]
    if not draw_contour:
        axe.hist2d(x, \
                   y, \
                   #extent=[np.min(xp_all), np.max(xp_all), np.min(yp_all), np.max(yp_all)],
                   cmap='binary',#color_dic[characterization_cat],#'jet'#,'YlOrBr'
                   #alpha=0.1,
                   bins=15
                    )
    if draw_contour:
        axe.scatter(x, \
        		   y, \
        		   c=z,\
        		   s=30,\
        		   #c= color_dic[title],\
        		   cmap='jet', \
        		   #edgecolor='',\
        		   #cmap= color_dic[title],\
        		   #'darkgoldenrod', \
        		   ##'slategray'
        		   alpha=0.20, \
        		   linewidths=0.2, \
                   #edgecolors= 'black',
        		   zorder=0, \
        		   marker='o',\
        		   label='')
    def density_estimation(m1, m2):
    	values = np.vstack([m1, m2])
    	kernel = gaussian_kde(values)                                                                 
    	sf = kernel.scotts_factor()
    	bwx = sf * np.std(m1)
    	bwy = sf * np.std(m2)
    	X, Y = np.mgrid[(np.min(m1)-3*bwx):(np.max(m1)+3*bwx):100j, (np.min(m2)-3*bwy):(np.max(m2)+3*bwy):100j]                                                     
    	positions = np.vstack([X.ravel(), Y.ravel()])                                                       
    	Z = np.reshape(kernel(positions).T, X.shape)
    	return X, Y, Z
    if draw_contour:
        X, Y, Z = density_estimation(gp, mp)
        axe.contour(X, Y, Z, linewidths = 1, alpha = .75)
    if scale == 'log':
    	axe.set_ylabel('PPANINI score (log)', fontsize=6)
    	axe.set_xlabel('Genomic priority (log)', fontsize=6)
    elif scale == 'sqrt':
    	axe.set_ylabel('PPANINI score (sqrt)', fontsize=6)
    	axe.set_xlabel('Genomic priority (sqrt)', fontsize=6)
    else:
    	axe.set_ylabel('PPANINI score', fontsize=6)
    	axe.set_xlabel('Genomic priority', fontsize=6)
    axe.set_ylabel('PPANINI score', fontsize=6)
    axe.get_xaxis().set_tick_params(which='both', labelsize=4,top='off',  direction='out')
    axe.get_yaxis().set_tick_params(which='both', labelsize=4, right='off', direction='out')
    axe.yaxis.set_label_position('left') 
    if title:
    	axe.set_title(title, loc='left', fontdict={'fontsize':'8','fontweight' :'bold'})
    try:
        axe.xlim([min(gp) , max(gp)])
        axe.ylim(min(mp), max(mp))
    except:
        plt.xlim([min(gp) , max(gp)])
        plt.ylim(min(mp), max(mp))
        
    axe.autoscale_view('tight')
    plt.tight_layout()
    
    if output_path:
    	plt.savefig(output_path+'.pdf', bbox_inches='tight', pad_inches = 0, dpi=300)
    	plt.savefig(output_path+'.png', bbox_inches='tight', pad_inches = 0, dpi=300)
    	plt.savefig(output_path+'.svgz', bbox_inches='tight', pad_inches = 0, dpi=300)
    return axe
def hexbin_plot_metagenomic_priority(axe, ppanini_table, table, no_uniq_genomes, title, scale = None):
	if axe is None:
		axe = plt.gca()
	#mp_gp = {}
	genes = ppanini_table.index
	#abund = ppanini_table['abundance']
	#prev = ppanini_table['prevalence']
	ppanini_score = ppanini_table['ppanini_score']
	#abund = np.array(abund)/max(abund)
	#prev = np.array(prev)/max(prev)
	
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
	gp = np.array(gp)/float(no_uniq_genomes)
	if scale == 'log':
		gp = np.log(gp)
		mp = np.log(mp)
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
		axe.set_ylabel('PPANINI score (log)',  fontsize=10)
		axe.set_xlabel('Genomic Priority (log)', fontsize=10)
	else:
		axe.set_ylabel('PPANINI score', fontsize=10)
		axe.set_xlabel('Genomic Priority', fontsize=10)
	axe.set_title(title, loc='left', fontdict={'fontsize':'10','fontweight' :'bold'})
	axe.get_xaxis().set_tick_params(which='both', labelsize=8,top='off',  direction='out')
	axe.get_yaxis().set_tick_params(which='both', labelsize=8, right='off', direction='out')
	axe.yaxis.set_label_position('left') 
	return im
def priority_scatter(args):
	diamond_output =  utilities.gene2genomes(args.diamond_output)
	ppanini_output = utilities.pd.DataFrame.from_csv(args.ppanini_output,
											 sep='\t', index_col=0, header =0)
	no_uniq_genomes = utilities.number_of_unique_genomes(args.diamond_output)
	scatter_plot_metagenomic_priority(None, ppanini_output, diamond_output, no_uniq_genomes, title = None,scale = None , output_path=args.path+'/'+args.outfile, size = args.size)


def scatter_plot_prev_abund(axe, ppanini_table, title, characterization_cat ='', essential_genes = '', xscale = 'log', yscale = 'log', output_path = None, size = 3, num_rand = None):
    if axe == None:
        fig, axe = plt.subplots(1, figsize=(size, size))
    #mp_gp = {}
    import random
    if not num_rand:
        num_rand = 5000#len(ppanini_table.index) #

    idxs = random.sample(range(len(ppanini_table.index)), min(num_rand, len(ppanini_table.index)))
    genes = ppanini_table.index[idxs]
    #abund = ppanini_table['abundance']
    ppanaini_scores = ppanini_table['ppanini_score'][idxs]
    #print ppanaini_scores
    prev = np.array(ppanini_table['alpha_prevalence'][idxs]) #alpha_prevalence
    abund = np.array(ppanini_table['mean_abundance'][idxs])
    if xscale == 'log':
        all_prev = np.log(np.array(ppanini_table['alpha_prevalence']))
    else: 
        all_prev = np.array(ppanini_table['alpha_prevalence'])
    if yscale == 'log':
        all_abund = np.log(np.array(ppanini_table['mean_abundance']))
    else:
        all_abund = np.array(ppanini_table['mean_abundance'])
    
    go_term = ppanini_table['GO'][idxs]
    
    #color_dic_cmap= {"Only Uniref": 'Blues', "Uniref and GO term":'Greens', 
    #            "Unannotated":'YlOrRd'}
    color_dic= {"Only Uniref": 'blue', "Uniref and GO term":'limegreen', 
                "Unannotated":'gold'}
    
    if xscale == 'log':
        prev = np.log(prev)
    if yscale == 'log':
        abund = np.log(abund)
    
    # Calculate the point density
    x = np.array(prev)#[val for (val,val2) in zip(prev, abund) if val != 'NaN' and val2 != 'NaN']
    y = np.array(abund) #[val2 for (val,val2) in zip(prev, abund) if val != 'NaN' and val2 != 'NaN']
    xy = np.vstack([x,y])
    #z = gaussian_kde(xy)(xy)   
    points_color = np.array([color_dic["Uniref and GO term"] if go_val == go_val else 
    color_dic["Only Uniref"] if gene.startswith('UniRef')  else color_dic["Unannotated"] for
    (go_val, gene) in zip(go_term, genes)])
    
    
    #display_essential = np.array([True if characterization_cat == "Uniref and GO term" and go_val == go_val  and gene in essential_genes else 
    #True if characterization_cat == "Only Uniref" and gene.startswith('UniRef') and go_val != go_val and gene in essential_genes else False for
    #(go_val, gene) in zip(go_term, genes)])
    display_essential = np.array([True if gene in essential_genes else  False for  gene in ppanini_table.index])

    
    points_marker = np.array(['^' if val >= 0.75 else 'v' if val < 0.75 else '+'
                               for val in ppanaini_scores])
    xp_all =  []
    yp_all =[]
    points_alpha = np.array([.35 if val >= 0.75 else .1 for val in ppanaini_scores])
    ''' 
    for xp, yp, mp, cp, ap in zip(x, y, points_marker, points_color, points_alpha):
        if xp == xp and yp == yp and xp and yp: #cp == color_dic[characterization_cat] and
            xp_all.append(xp)
            yp_all.append(yp)
           
    axe.hist2d(xp_all, \
                   yp_all, \
                   #extent=[np.min(xp_all), np.max(xp_all), np.min(yp_all), np.max(yp_all)],
                   cmap=color_dic[characterization_cat],#'jet'#,'YlOrBr'
                   #alpha=0.1,
                   bins=15
                   )
    
    
    axe.hexbin(xp_all, \
                   yp_all, \
                   #extent=[np.min(xp_all), np.max(xp_all), np.min(yp_all), np.max(yp_all)],
                   cmap=color_dic[characterization_cat],#'jet'#,'YlOrBr'
                   #alpha=0.1,
                   gridsize=15,
                   marginals=False)
    
    '''
    for xp, yp, mp, cp, ap in zip(x, y, points_marker, points_color, points_alpha):
        if cp != color_dic["Unannotated"]:
            continue
        axe.scatter(xp, \
                   yp, \
                   s=20,\
                   c= cp,\
                   cmap='jet', \
                   #edgecolor='',\
                   #cmap= color_dic[title],\
                   #'darkgoldenrod', \
                   ##'slategray'
                   alpha=ap, \
                   linewidths=0.2, \
                   #edgecolors= 'black',
                   zorder=0, \
                   marker=mp,\
                   label='')
    for xp, yp, mp, cp, ap in zip(x, y, points_marker, points_color, points_alpha):
        if cp != color_dic["Only Uniref"]:
            continue
        axe.scatter(xp, \
                   yp, \
                   s=20,\
                   c= cp,\
                   cmap='jet', \
                   #edgecolor='',\
                   #cmap= color_dic[title],\
                   #'darkgoldenrod', \
                   ##'slategray'
                   alpha=ap, \
                   linewidths=0.2, \
                   #edgecolors= 'black',
                   zorder=0, \
                   marker=mp,\
                   label='')
    for xp, yp, mp, cp, ap in zip(x, y, points_marker, points_color, points_alpha):
        if cp != color_dic["Uniref and GO term"]:
            continue
        axe.scatter(xp, \
                   yp, \
                   s=20,\
                   c= cp,\
                   cmap='jet', \
                   #edgecolor='',\
                   #cmap= color_dic[title],\
                   #'darkgoldenrod', \
                   ##'slategray'
                   alpha=ap, \
                   linewidths=0.2, \
                   #edgecolors= 'black',
                   zorder=0, \
                   marker=mp,\
                   label='')
    #display_essential = display_essential[random.sample(range(len(ppanini_table.index)), min(2000, len(display_essential)))]
    for xp, yp, de in zip(all_prev, all_abund, display_essential):
        if de:
            axe.scatter(xp, \
                       yp, \
                       s= 20,\
                       c= 'silver',\
                       cmap='jet', \
                       #edgecolor='',\
                       #cmap= color_dic[title],\
                       #'darkgoldenrod', \
                       ##'slategray'
                       alpha=0.25, \
                       linewidths=0.25, \
                       edgecolors= 'black',
                       #zorder=0, \
                       marker='*',\
                       label='')
    #print np.nanpercentile(x, 75), np.nanpercentile(y, 75)
    axe.text(0.55, 0.1, "number of gene families: %d" % (len(ppanini_table.index)), fontsize=5, va="center", ha="center", transform=axe.transAxes)
    axe.axvline(x= np.nanpercentile(all_prev[~np.isnan(all_prev)], 75), color='black', linestyle='--', linewidth = .75)
    axe.axhline(y= np.nanpercentile(all_abund[~np.isnan(all_abund)], 75), color='black', linestyle='--', linewidth = .75)

    def density_estimation(m1, m2):
        values = np.vstack([m1, m2])
        kernel = gaussian_kde(values)                                                                 
        sf = kernel.scotts_factor()
        bwx = sf * np.std(m1)
        bwy = sf * np.std(m2)
        X, Y = np.mgrid[(np.min(m1)-3*bwx):(np.max(m1)+3*bwx):100j, (np.min(m2)-3*bwy):(np.max(m2)+3*bwy):100j]                                                     
        positions = np.vstack([X.ravel(), Y.ravel()])                                                       
        Z = np.reshape(kernel(positions).T, X.shape)
        return X, Y, Z
    '''import math
    test = [math.isnan(a) or math.isnan(b) for a,b in zip(prev,abund)]
    new_X= [a for a,b in zip (prev,test) if ~b]
    new_Y= [a for a,b in zip (abund,test) if ~b]
    print new_X , new_Y'''
    #X, Y, Z = density_estimation(prev, abund)
    #axe.contour(X, Y, Z, linewidths = .5, alpha = .75)
    if yscale == 'log':
        axe.set_ylabel('Relative abundance (log)', fontsize=6)
    if xscale == 'log':
        axe.set_xlabel('prevalence (log)', fontsize=6)
    elif yscale == 'sqrt':
        axe.set_ylabel('Relative abundance (sqrt)', fontsize=6)
    elif xscale == 'sqrt':
        axe.set_xlabel('prevalence (sqrt)', fontsize=6)
    else:
        axe.set_ylabel('Relative abundance', fontsize=6)
        axe.set_xlabel('prevalence', fontsize=6)
    #axe.set_ylabel('Relative abundance', fontsize=6)
    axe.get_xaxis().set_tick_params(which='both', labelsize=4,top='off',  direction='out')
    axe.get_yaxis().set_tick_params(which='both', labelsize=4, right='off', direction='out')
    axe.yaxis.set_label_position('left') 
    if title:
        axe.set_title(title, loc='left', fontdict={'fontsize':'8','fontweight' :'bold'})
    try:
        axe.xlim([min(prev) , max(prev)])
        axe.ylim(min(abund), max(abund))
    except:
        plt.xlim([min(prev) , max(prev)])
        plt.ylim(min(abund), max(abund))
        
    #axe.autoscale_view('tight')
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path+'.pdf', bbox_inches='tight', pad_inches = 0, dpi=300)
        plt.savefig(output_path+'.png', bbox_inches='tight', pad_inches = 0, dpi=300)
        plt.savefig(output_path+'.svgz', bbox_inches='tight', pad_inches = 0, dpi=300)
    return axe	

def get_args( ):
    parser = argparse.ArgumentParser(
        description="PPANINI plotting tool",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument( "-i", "--ppanini-output",
                         metavar = "<input table>",
                         dest = "ppanini_output",
                         required = False, 
                         help="PPANINI output table", )
    parser.add_argument( "-m8", "--diamond-output ",
                         dest = "diamond_output",
                         metavar = "<feature id>",
                         help="a mapping file of gene-metagenom)", )
    parser.add_argument( "--master-plot", 
                         dest = "master",
                         help="plotting master figure of the paper",
                         action="store_true", )
    parser.add_argument( "--path",
                         dest = "path",
                         help="path for inputs and/or outputs", )
    parser.add_argument( "--outfile",
                         dest = "outfile",
                         help="output file", )
    parser.add_argument( "--size",
                         dest = "size",
                         default =3,
                         type =int,
                         help="size of the plot in inches", )
    
    
    
    return parser.parse_args()
def main():
	user_args = get_args()
	if user_args.master:
	    master_plot(user_args.path, size = user_args.size)
	else:
		priority_scatter(user_args)
		
if __name__ == '__main__':
	main()
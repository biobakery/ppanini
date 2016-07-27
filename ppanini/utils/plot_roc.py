
#print(__doc__)
 
import numpy as np
import matplotlib
from random import randint
from scipy import random
from sys import argv
from os.path import basename

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.cross_validation import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
import scipy.stats
import math
from scipy.stats import percentileofscore
import sys
import csv
from . import plot_metagenome_genome 
from .. import utilities
#sys.path.append('/Users/rah/Documents/Hutlab/ppanini')#/n/hutlab12_nobackup/data/ppanini/ppanini')

from .. import config
from . import plot_metagenome_genome
def fpr_tpr_genome(metagenomic_table, no_uniq_genomes, essential_genes):
    # calculate genomic score 
    gp = np.zeros(len(metagenomic_table))
    uniq_genomes = []
    genes = metagenomic_table.keys()#['gene']
    #print metagenomic_table.values()[1:100]
    for i in range(len(genes)):
        gene = genes[i]
        gp[i] += [len(metagenomic_table[gene])]
        
    gp = np.array(gp)/float(no_uniq_genomes)
    #print gp 
    ground_truth = [1 if (gene_id  in essential_genes) else 0 for gene_id in genes ]
    #print ground_truth
    fpr, tpr, _  = roc_curve( ground_truth, gp, pos_label = 1)
    return fpr, tpr
    
    
def main():
    plot_g = True
    plot_mg = True
    roc_info = [] 
    config.output_folder = '/Users/rah/Documents/Hutlab/ppanini/myOutput2'
    with open('/Users/rah/Documents/Hutlab/UniRef90_output_deg_p_gene.m8') as f:
        lines = f.read().splitlines()
    essential_genes_uniref90_id_deg = [line.split('\t')[1] for line in lines]
    with open('/Users/rah/Documents/Hutlab/UniRef90_output_299_gene.m8') as f:
        lines = f.read().splitlines()
    essential_genes_uniref90_id_299_eco = [line.split('\t')[1] for line in lines]
    essential_genes = essential_genes_uniref90_id_299_eco + essential_genes_uniref90_id_deg
    niches = ["stool", "buccal", "nares", "soil"]
    categories= ["metagenomic", "genomic"] 
    Betas = [.5]
    
    for niche in niches:
        for beta in Betas:
            for category in categoreis:
                if category == "metagenomic":
                    try:
                        
                        fpr, tpr = get_fpr_tpr(input_ppanini='/Users/rah/Dropbox (Huttenhower Lab)/PPANINI/ASM/ranks/nares.ranks',\
                                               essential_genes= essential_genes, beta =0.1)
                                               #input_ppanini='/Users/rah/Documents/Hutlab/ppanini//output_tables/nares_table.txt',\
                        roc_info.append(['Anterior nares - Metagenomic Priority beta = .1',fpr, tpr])
            
                    except ValueError:
                        print "ValueError for "+niche+ " "+ category + " " + beta+" calculation"  
                if category == "genomic":
                    try:
                        fpr, tpr = get_fpr_tpr(input_ppanini='/Users/rah/Documents/Hutlab/ppanini/PF_final_gene_centroids_table/PF_final_gene_centroids_table_imp_centroid_prev_abund.txt',\
                                                  essential_genes= essential_genes, beta =.5)     
                        roc_info.append(['Posterior fornix - Metagenomic Priority',fpr, tpr])
                    
                    except ValueError:
                        print "ValueError for Posterior fornix roc calculation"
# Stool
    if plot_mg:
        try:
            
            fpr, tpr = get_fpr_tpr(input_ppanini='/Users/rah/Dropbox (Huttenhower Lab)/PPANINI/ASM/ranks/stool.ranks',\
                                   essential_genes= essential_genes, beta =0.5)
                                   #input_ppanini='/Users/rah/Documents/Hutlab/ppanini//output_tables/nares_table.txt',\
            roc_info.append(['Stool - Metagenomic Priority',fpr, tpr])
            
            fpr, tpr = get_fpr_tpr(input_ppanini='/Users/rah/Dropbox (Huttenhower Lab)/PPANINI/ASM/ranks/buccal.ranks',\
                                   essential_genes= essential_genes, beta =.5)
                                   #input_ppanini='/Users/rah/Documents/Hutlab/ppanini//output_tables/nares _table.txt',\
            roc_info.append(['Buccal mucosa - Metagenomic Priority',fpr, tpr])
            
            fpr, tpr = get_fpr_tpr(input_ppanini='/Users/rah/Dropbox (Huttenhower Lab)/PPANINI/ASM/ranks/fornix.ranks',\
                                   essential_genes= essential_genes, beta =.5)
                                   #input_ppanini='/Users/rah/Documents/Hutlab/ppanini//output_tables/nares _table.txt',\
            roc_info.append(['Posterior fornix - Metagenomic Priority',fpr, tpr])
            fpr, tpr = get_fpr_tpr(input_ppanini='/Users/rah/Dropbox (Huttenhower Lab)/PPANINI/ASM/ranks/nares.ranks',\
                                   essential_genes= essential_genes, beta =.5)
                                   #input_ppanini='/Users/rah/Documents/Hutlab/ppanini//output_tables/nares _table.txt',\
            roc_info.append(['Anterior nares - Metagenomic Priority',fpr, tpr])
            fpr, tpr = get_fpr_tpr(input_ppanini='/Users/rah/Dropbox (Huttenhower Lab)/PPANINI/ASM/ranks/soil.ranks',\
                                   essential_genes= essential_genes, beta =.5)
                                   #input_ppanini='/Users/rah/Documents/Hutlab/ppanini//output_tables/stool _table.txt',\
            roc_info.append(['Prairie soil - Metagenomic Priority',fpr, tpr])

            
        except ValueError:
            print "ValueError for Buccal  roc calculation"  
        
    try:
        fpr, tpr = get_fpr_tpr(input_ppanini='/Users/rah/Documents/Hutlab/ppanini/output_tables/Stool_table.txt',\
                                  essential_genes= essential_genes, beta =.5)     
        roc_info.append(['Stool - Metagenomic Priority',fpr, tpr])      
    except ValueError:
        print "ValueError for Stool roc calculation"
    if plot_g:
        try:
            metagenomic_table1, ppanini_output1, no_uniq_genomes1  = utilities.read_data('/Users/rah/Documents/Hutlab/ppanini/PARSED_BLAST_RESULTS/stool_mg.m8', '/Users/rah/Documents/Hutlab/ppanini/output_tables/stool_table.txt')
            #print  metagenomic_table1
            fpr, tpr =fpr_tpr_genome(metagenomic_table1, no_uniq_genomes1, essential_genes)
            roc_info.append(['Stool - Genomic Priority',fpr, tpr])
        except ValueError:
            print "ValueError for Stool - Genomic Priority  roc calculation"
#Buccal mucosa    
    if plot_mg:
        try:
            fpr, tpr = get_fpr_tpr(input_ppanini='/Users/rah/Documents/Hutlab/ppanini/output_tables/BM_table.txt',\
                                      essential_genes= essential_genes, beta =.5)     
            roc_info.append(['Buccal mucosa - Metagenomic Priority',fpr, tpr])      
        except ValueError:
            print "ValueError for Buccal mucosa roc calculation"
    if plot_g:
        try:
            metagenomic_table3, ppanini_output3, no_uniq_genomes3 = utilities.read_data('/Users/rah/Documents/Hutlab/ppanini/PARSED_BLAST_RESULTS/BM_mg.m8', '/Users/rah/Documents/Hutlab/ppanini/output_tables/BM_table.txt')
            fpr, tpr = fpr_tpr_genome(metagenomic_table3, no_uniq_genomes3, essential_genes)
            roc_info.append(['Buccal mucosa - Genomic Priority',fpr, tpr])
        except ValueError:
            print "ValueError for Buccal mucosa  roc calculation"
#Anterior nares
    
    if plot_mg:
        try:
            fpr, tpr = get_fpr_tpr(input_ppanini='/Users/rah/Documents/Hutlab/ppanini/output_tables/AN_table.txt',\
                                      essential_genes= essential_genes, beta =.5)     
            roc_info.append(['Anterior nares - Metagenomic Priority',fpr, tpr])
        except ValueError:
            print "ValueError for Anterior nares roc calculation" 
        try:
            fpr, tpr = get_fpr_tpr(input_ppanini='/Users/rah/Documents/Hutlab/ppanini/output_tables/AN_table_beta25.txt',\
                                      essential_genes= essential_genes, beta =.25)     
            roc_info.append(['Anterior nares - Metagenomic Priority Beta = .25',fpr, tpr])
        except ValueError:
            print "ValueError for Anterior nares roc calculation" 
        try:
            fpr, tpr = get_fpr_tpr(input_ppanini='/Users/rah/Documents/Hutlab/ppanini/output_tables/AN_table_beta75.txt',\
                                      essential_genes= essential_genes, beta =.25)     
            roc_info.append(['Anterior nares - Metagenomic Priority Beta = .75',fpr, tpr])
        except ValueError:
            print "ValueError for Anterior nares roc calculation"
            
        try:
            fpr, tpr = get_fpr_tpr(input_ppanini='/Users/rah/Documents/Hutlab/ppanini/output_tables/AN_table_from_gene_catalog.txt',\
                                      essential_genes= essential_genes, beta =.25)     
            roc_info.append(['Anterior nares - Metagenomic Priority from gene catlog',fpr, tpr])
        except ValueError:
            print "ValueError for Anterior nares gene-catalog roc calculation"
    if plot_g:
        try:
            metagenomic_table2, ppanini_output2, no_uniq_genomes2 = utilities.read_data('/Users/rah/Documents/Hutlab/ppanini/PARSED_BLAST_RESULTS/AN_mg.m8', '/Users/rah/Documents/Hutlab/ppanini/output_tables/AN_table.txt')
            fpr, tpr = fpr_tpr_genome(metagenomic_table2, no_uniq_genomes2, essential_genes)
            roc_info.append(['Anterior nares - Genomic Priority',fpr, tpr])
        except ValueError:
            print "ValueError for Anterior nares  roc calculation"
   
    if plot_mg:
        try:
            metagenomic_table4, ppanini_output4, no_uniq_genomes4 = utilities.read_data('/Users/rah/Documents/Hutlab/ppanini/PARSED_BLAST_RESULTS/PF_mg.m8', '/Users/rah/Documents/Hutlab/ppanini/output_tables/PF_table.txt')
            fpr, tpr = fpr_tpr_genome(metagenomic_table4, no_uniq_genomes4,essential_genes)
            roc_info.append(['Posterior fornix Genomic Priority',fpr, tpr])
        
        except ValueError:
            print "ValueError for Posterior fornix  roc calculation"
    
    if plot_g:
        try:
            fpr, tpr = get_fpr_tpr(input_ppanini='/Users/rah/Documents/Hutlab/ppanini/PF_final_gene_centroids_table/PF_final_gene_centroids_table_imp_centroid_prev_abund.txt',\
                                      essential_genes= essential_genes, beta =.5)     
            roc_info.append(['Posterior fornix - Metagenomic Priority',fpr, tpr])
        
        except ValueError:
            print "ValueError for Posterior fornix roc calculation"
    try:
        print roc_info
        roc_plot(roc_info, figure_name=config.output_folder+'/nares_roc_plot_ppanini') 
    except ValueError:
        print "ValueError for general roc plot"
    print "The evaluation is successfully done!" 
def get_fpr_tpr(input_ppanini, essential_genes, beta =.5):
    fpr = dict()
    tpr = dict()
    true = dict()
    score = dict()
    #config.input_table = '/Users/rah/Documents/Hutlab/stool_ppanini050715.txt'#'/n/hutlab12_nobackup/data/ppanini/DATA/PPANINI_INPUT/stool_ppanini.txt' 
    #config.uclust_file = '/Users/rah/Documents/Hutlab/stool_ppanini/stool_final_clusters.uc'
    
    #config.output_folder = '/Users/rah/Documents/Hutlab/ppanini/myOutput2'
    
    #uniref_id_list = []
    #ground_truth = [1 if (uniref_id in essential_genes_uniref_id) else 0 for uniref_id in uniref_id_list ] # this an example for each gene if it's important use 1 otherwise 0
    #print config.input_table
    if config.verbose =='DEBUG':
        print "PPANINI Score and ROC Plot!"
    #ppanini.run()
    '''try:
        input_ppanini = str(sys.argv[1])
    except ImportError:
       input_ppanini  = '/Users/rah/Documents/Hutlab/ppanini/myOutput2/AN_final_centroid_prev_abund.txt '
       #stool_ppanini050715_imp_centroid_prev_abund.txt'
    '''
    '''with open(input_ppanini) as f:
        lines2 = f.read().splitlines()
    config.centroids_list = [line.split('\t')[0] for line in lines2]
    prev = [line.split('\t')[1] for line in lines2]
    abun = [line.split('\t')[4] for line in lines2]
    ppanini_score = [line.split('\t')[2] for line in lines2]
    #ppanini_score = [0 for line in lines2]'''
    ppanini_table = utilities.read_ppanini_imp_genes_table(input_ppanini)
    #print ppanini_table['ppanini_score'][1:10]
    
    n = len(ppanini_table['genes'])-1
    config.centroids_list = ppanini_table['genes'][0:n]#config.centroids_list[1:n]
    
    #print config.centroids_list[0:10]
    prev = ppanini_table['prevalence_rank'][0:n]#prev[1:n]
    prev = [float(val) for val in prev]
    
    sorted_prev = sorted(prev)
    abun = ppanini_table['abundance_rank'][0:n]#abun[1:n]
    ppanini_table['ppanini_score'] = [1/((1/(beta*prev[i])+(1/((1-beta)* abun[i])))) for i in range(len(prev))] #
    ppanini_score = ppanini_table['ppanini_score'][0:n]#ppanini_score[1:n]
    abun = [float(val) for val in abun]
    sorted_abun = sorted(abun)
    #print prev[0:101]
    #print abun[0:]
    ground_truth = [1 if (gene_id  in essential_genes) else 0 for gene_id in config.centroids_list ]
    #print ground_truth
    eval_file = open(config.output_folder+"/"+basename(input_ppanini)+'.txt', 'w') 
    csvw = csv.writer(eval_file, csv.excel_tab, delimiter='\t')
    csvw.writerow(["centroid", "prevalence", "abundance",  "ppanini score","is_essential"])
    for i in range(len(config.centroids_list)):
        csvw.writerow([config.centroids_list[i], prev[i], abun[i], ppanini_score[i] ,ground_truth[i]])
    eval_file.close()
        
    
    #for b in ["DEG", "ECO"]:#range(5, 6, 1):
    #beta = beta# float(b/10.0)         
    #config.beta = beta
    #if b == "DEG":
    #    essential_genes = essential_genes_uniref90_id_deg
    #else:
    #    essential_genes = essential_genes_uniref90_id_299_eco
    ground_truth = [1 if (gene_id  in essential_genes) else 0 for gene_id in config.centroids_list ]
    #scipy.stats.rankdata()
    score[beta] = ppanini_score#[1/((1/((beta)*scipy.stats.percentileofscore(sorted_prev, prev[i])))+\
                  #(1/((1-beta)*scipy.stats.percentileofscore(sorted_abun, abun[i])))) for i in range(len(prev))]
    #prioritize_results = ppanini.prioritize_centroids()
    #print "prioritize_results", prioritize_results
    #print config.centroids_list
    #print prioritize_gene_results
    #uniref90_id_list = [id.split('|')[1] for id in config.centroid_prev_abund]
    #uniref50_id_list = [id.split('|')[2] for id in config.centroid_prev_abund] 
    #prioritized_ids = config.centroids_list 
    #print "Essential genes: ",essential_genes_uniref_id
    #print "**********************************************************************"
    #print "Centroids list: ", config.centroids_list
    
    '''ground_truth = []
    for gene_id in config.centroids_list:
        if gene_id in essential_genes_uniref_id:
            print gene_id
            ground_truth.append(0)
        else:
            ground_truth.append(1)
    '''
    
    true[beta] = ground_truth
    #print "ground_truth", scipy.stats.rankdata(score[beta])
    '''if config.niche_flag:
        score[beta] =[config.centroid_prev_abund[gene_id]['ppanini_score'][config.centroid_prev_abund[gene_id]['ppanini_score'].keys()[0]] for gene_id in config.centroids_list ]
    else:
        score[beta] =[config.centroid_prev_abund[gene_id]['ppanini_score'] for gene_id in config.centroids_list ] # 
    '''
    #legend_tag = b
    #print score[0:10]
    assert(len(true[beta])==len(score[beta])) 
    fpr[beta], tpr[beta], _  = roc_curve( true[beta], score[beta], pos_label = 1)
    return fpr[beta], tpr[beta]
    
                                     
def roc_plot(roc_info=None, figure_name='roc_plot_ppanini'):
    """
    =======================================
    Receiver Operating Characteristic (ROC)
    =======================================

    Parameters
    ------------
        roc_info : List of lists (method name, fpr, tpr) 

    Returns 
    -------------
        Plots ROC curves
        save as pdf and show  
    """
    
    # Compute ROC curve and ROC area for each class
    fpr = dict()
    tpr = dict()
    roc_name = ''
    roc_auc = dict()
    for i in range(len(roc_info)):
        # print ('Hi', (roc_info[i][1]))
        fpr[roc_info[i][0]] = roc_info[i][1]
        # print ((roc_info[i][1])[0])
        tpr[roc_info[i][0]] = roc_info[i][2]
        roc_auc[roc_info[i][0]] = auc(fpr[roc_info[i][0]], tpr[roc_info[i][0]] )
        roc_name += '_' + roc_info[i][0] 
        
    # Plot ROC curve
    fig, axe = plt.subplots(figsize=(5, 5 ), dpi=300)#, sharex=False, sharey=False)
    #fig.set_size_inches(1, 10)
    #plt.figure(dpi= 300, figsize=(4, 4))
    for i in range(len(roc_info)):
        params = {'legend.fontsize': 6
        }
        #mpl.rcParams['lines.linewidth'] = 2
        plt.rcParams.update(params)
        axe.plot(fpr[roc_info[i][0]], tpr[roc_info[i][0]],  label='{0} (area = {1:0.2f})'
                                       ''.format(str(roc_info[i][0]), roc_auc[roc_info[i][0]]))   
    axe.plot([0, 1], [0, 1], 'k--')
    axe.set_xlim([0.0, 1.0])
    axe.set_ylim([0.0, 1.05])
    axe.legend(loc="lower right")
    axe.set_ylabel('True Positive Rate', fontsize = 10)
    axe.set_xlabel('False Positive Rate', fontsize = 10)
    axe.get_xaxis().set_tick_params(which='both', labelsize=8,top='off',  direction='out')
    axe.get_yaxis().set_tick_params(which='both', labelsize=8, right='off', direction='out')
    axe.yaxis.set_label_position('left') 
    axe.set_title('Receiver operating characteristic', fontsize=12, fontweight='bold')
    # plt.savefig('./test/'+roc_name+'foo.pdf')
    plt.tight_layout()
    plt.savefig(figure_name + '.pdf')
if __name__ == '__main__':
    main()
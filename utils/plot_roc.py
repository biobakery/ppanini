
#print(__doc__)
 
import numpy as np
import matplotlib
from random import randint
from scipy import random
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

sys.path.append('/Users/rah/Documents/Hutlab/ppanini')#/n/hutlab12_nobackup/data/ppanini/ppanini')
import ppanini
from src import config
 
def evaluation_multi_roc():
    fpr = dict()
    tpr = dict()
    true = dict()
    score = dict()
    roc_info = []
    # list of truth about data or associations, here, is this an important gene?
    '''
    config.input_table = '/Users/rah/Documents/Hutlab/output_ppanini.txt'#'/n/hutlab12_nobackup/data/ppanini/DATA/PPANINI_INPUT/stool_ppanini.txt' 
    config.uclust_file = '/Users/rah/Documents/Hutlab/output.uc'
    '''
    config.input_table = '/Users/rah/Documents/Hutlab/stool_ppanini050715.txt'#'/n/hutlab12_nobackup/data/ppanini/DATA/PPANINI_INPUT/stool_ppanini.txt' 
    config.uclust_file = '/Users/rah/Documents/Hutlab/stool_ppanini/stool_final_clusters.uc'
    
    config.output_folder = 'myOutput2'
    #with open('/Users/rah/Documents/UniRef50_299_genes.txt') as f:
    #    essantial_genes_uniref50_id = f.read().splitlines()
    with open('/Users/rah/Documents/Hutlab/UniRef90_output_deg_p_gene.m8') as f:
        lines = f.read().splitlines()
    essantial_genes_uniref90_id_deg = [line.split('\t')[1] for line in lines]
    with open('/Users/rah/Documents/Hutlab/UniRef90_output_299_gene.m8') as f:
        lines = f.read().splitlines()
    essantial_genes_uniref90_id_299_eco = [line.split('\t')[1] for line in lines]
    #config.essantial_genes_uniref90_id = essantial_genes_uniref90_id
    #essantial_genes_uniref_id = essantial_genes_uniref90_id +essantial_genes_uniref50_id
    #print len(essantial_genes_uniref90_id)
    uniref_id_list = []
    #ground_truth = [1 if (uniref_id in essantial_genes_uniref_id) else 0 for uniref_id in uniref_id_list ] # this an example for each gene if it's important use 1 otherwise 0
    #print config.input_table
    if config.verbose =='DEBUG':
        print "Start Evaluating PPANINI Score !!!"
    #ppanini.run()
    with open('/Users/rah/Documents/Hutlab/ppanini/myOutput2/stool_ppanini050715_imp_centroid_prev_abund.txt') as f:
        lines2 = f.read().splitlines()
    config.centroids_list = [line.split('\t')[0] for line in lines2]
    prev = [line.split('\t')[2] for line in lines2]
    abun = [line.split('\t')[3] for line in lines2]
    ppanini_score = [line.split('\t')[1] for line in lines2]
    n =1000
    config.centroids_list = config.centroids_list[1:n]
    #print config.centroids_list[0:200]
    prev = prev[1:n]
    prev = [float(val) for val in prev]
    
    sorted_prev = sorted(prev)
    abun = abun[1:n]
    abun = [float(val) for val in abun]
    sorted_abun = sorted(abun)
    #print prev[0:101]
    #print abun[0:]
    ground_truth = [1 if (gene_id  in essantial_genes_uniref90_id_299_eco and\
                           gene_id in essantial_genes_uniref90_id_deg) else 0 for gene_id in config.centroids_list ]
    print ground_truth
    eval_file = open('/Users/rah/Documents/Hutlab/ppanini/eval_result.txt', 'w') 
    csvw = csv.writer(eval_file, csv.excel_tab, delimiter='\t')
    csvw.writerow(["centroid", "ppanini score","prevalence", "abundance",  "ground true"])
    for i in range(len(config.centroids_list)):
        csvw.writerow([config.centroids_list[i], prev[i], abun[i], ppanini_score[i] ,ground_truth[i]])
    eval_file.close()
        
    
    for b in range(1, 2, 1):
        beta = float(b/10.0)         
        config.beta = beta
        #scipy.stats.rankdata()
        score[beta] = [1/((1/((beta)*scipy.stats.percentileofscore(sorted_prev, prev[i])))+\
                      (1/((1-beta)*scipy.stats.percentileofscore(sorted_abun, abun[i])))) for i in range(len(prev))]
        #prioritize_results = ppanini.prioritize_centroids()
        #print "prioritize_results", prioritize_results
        #print config.centroids_list
        #print prioritize_gene_results
        #uniref90_id_list = [id.split('|')[1] for id in config.centroid_prev_abund]
        #uniref50_id_list = [id.split('|')[2] for id in config.centroid_prev_abund] 
        #prioritized_ids = config.centroids_list 
        #print "Essential genes: ",essantial_genes_uniref_id
        #print "**********************************************************************"
        #print "Centroids list: ", config.centroids_list
        
        '''ground_truth = []
        for gene_id in config.centroids_list:
            if gene_id in essantial_genes_uniref_id:
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
        # score[beta] =[beta*2/(i+1) for i in range(len(true[beta]))]
        print "score", score[beta]
        assert(len(true[beta])==len(score[beta])) 
        fpr[beta], tpr[beta], _  = roc_curve( true[beta], score[beta], pos_label = 1)
        roc_info.append([str(beta),fpr[beta], tpr[beta]])
        #print roc_info
        #break;
    try:
        roc_plot(roc_info) 
    except ValueError:
        print "ValueError for roc plot"
        
                                     
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
        fpr[roc_info[i][0]] = roc_info[i][1]
        tpr[roc_info[i][0]] = roc_info[i][2]
        roc_auc[roc_info[i][0]] = auc(fpr[roc_info[i][0]], tpr[roc_info[i][0]] )
        roc_name += '_' + roc_info[i][0] 
        
    # Plot ROC curve
    plt.figure(dpi= 300, figsize=(4, 4))
    for i in range(len(roc_info)):
        params = {'legend.fontsize': 6,
        'legend.linewidth': 2}
        plt.rcParams.update(params)
        plt.plot(fpr[roc_info[i][0]], tpr[roc_info[i][0]],  label='{0} (area = {1:0.2f})'
                                       ''.format(str(roc_info[i][0]), roc_auc[roc_info[i][0]]))   
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic')
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(figure_name + '.pdf')
    #plt.show()
    # return plt
if __name__ == '__main__':
    evaluation_multi_roc()
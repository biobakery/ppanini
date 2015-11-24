
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

import sys
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
    config.input_table = '/Users/rah/Documents/Hutlab/stool_ppanini050715.txt'#'/n/hutlab12_nobackup/data/ppanini/DATA/PPANINI_INPUT/stool_ppanini.txt' 
    config.uclust_file = '/Users/rah/Documents/Hutlab/stool_ppanini/stool_final_clusters.uc'
    config.output_folder = 'myOutpu3'
    with open('/Users/rah/Documents/UniRef50_299_genes.txt') as f:
        essantial_genes_uniref50_id = f.read().splitlines()
    with open('/Users/rah/Documents/UniRef90_299_genes.txt') as f:
        essantial_genes_uniref90_id = f.read().splitlines()
    essantial_genes_uniref_id = essantial_genes_uniref90_id +essantial_genes_uniref50_id
    #print essantial_genes_uniref_id
    uniref_id_list = []
    #ground_truth = [1 if (uniref_id in essantial_genes_uniref_id) else 0 for uniref_id in uniref_id_list ] # this an example for each gene if it's important use 1 otherwise 0
    #print config.input_table
    print "Start running!"
    ppanini.run()
    for b in range(2, 3, 1):
        beta = float(b/10.0)          
        config.beta = beta
        prioritize_results = ppanini.prioritize_centroids()
        #print config.centroids_list
        #print prioritize_gene_results
        #uniref90_id_list = [id.split('|')[1] for id in config.centroid_prev_abund]
        #uniref50_id_list = [id.split('|')[2] for id in config.centroid_prev_abund] 
        #prioritized_ids = config.centroids_list 
        
        ground_truth = [1 if (gene_id in essantial_genes_uniref_id) else 0 for gene_id in config.centroids_list ]
        true[beta] = ground_truth
        
        #true[beta] = [ 0  ,  1,  0 ,  1,  1 ,  1,  0 ,  1,  1 ,1,  1,  1, 0 , 1, 0 ,  1, 1 , 1, 1 ,  0,  1  ]
        score[beta] =[config.centroid_prev_abund[gene_id]['ppanini_score'] for gene_id in config.centroids_list ] # 
        # score[beta] =[beta*2/(i+1) for i in range(len(true[beta]))]
        fpr[beta], tpr[beta], _  = roc_curve( true[beta], score[beta], pos_label = 1)
        roc_info.append([str(beta),fpr[beta], tpr[beta]])
        #print roc_info
        #break;
    roc_plot(roc_info)                              
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
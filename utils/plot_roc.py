
#print(__doc__)
import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.cross_validation import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier

import sys
sys.path.append('/Users/rah/Documents/Hutlab/ppanini')
import ppanini
from src import config
 
def evaluation_multi_roc():
    fpr = dict()
    tpr = dict()
    true = dict()
    score = dict() 
    
    # list of truth about data or associations, here, is this an important gene?
    #config.input_table = '/Users/rah/Documents/Hutlab/stool_ppanini.txt'#'/n/hutlab12_nobackup/data/ppanini/DATA/PPANINI_INPUT/stool_ppanini.txt' 
    truth = [1, 0, 1, 0, 0] # this an example for each gene if it's important use 1 otherwise 0
    print config.input_table
    ppanini.run()
    for b in range(.2, 1.0, .1):
         
        true[b] = truth 
        
        # score of of each association, ppanini score here   
        
        # ppanini_score should be implemented an returns list of scores.
        # Each gene has a score in the list in the same order as truth
        # get_important_centroids could be used by returning 
        # a list of 1, for prioritized, and 0, for unprioritized genes  
        #score[b] = ppanini_score(b)
        ppanini.prioritize_centroids
        continue
        # an example for scores
        score[b] = [.6, .25, .55, .15, .18] # 
        
        fpr[new_method], tpr[new_method], _  = roc_curve( true[new_method], score[new_method], pos_label = 1)
        roc_info.append([str(b),fpr[b], tpr[b]])
    
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
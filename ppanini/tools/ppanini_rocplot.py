
#print(__doc__)
 
import numpy as np
import matplotlib
from random import randint
from scipy import random
from sys import argv
import argparse
import pandas as pd
from os.path import basename
try:
    import matplotlib
    matplotlib.use( "Agg" )
    import matplotlib.pyplot as plt
    from matplotlib import colors
    import matplotlib.pyplot as plt
except:
    sys.exit( "This script requires the Python scientific stack: matplotlib." )
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

from matplotlib import font_manager
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("error")
    try:
        font_file = font_manager.findfont(font_manager.FontProperties(family='Arial'))
        matplotlib.rcParams["font.family"] = "Arial"
    except UserWarning:
        pass 

#from . import plot_metagenome_genome 
from .. import utilities

from .. import config
#from . import plot_metagenome_genome
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
    
    
def load_multiple_roc(args):
    roc_info = [] 
    path = args.path
   #UniRef90_output_deg_p_gene.m8'
    with open(args.essentail_gene_file1) as f:
        lines = f.read().splitlines()
    essential_genes_uniref90_id_deg = [line.split('\t')[1] for line in lines]
    
    #UniRef90_output_299_gene.m8'
    with open(args.essentail_gene_file2) as f:
        lines = f.read().splitlines()
    essential_genes_uniref90_id_299_eco = [line.split('\t')[1] for line in lines]
    essential_genes = essential_genes_uniref90_id_299_eco + essential_genes_uniref90_id_deg
    args.plot_mg = True
    args.plot_g = False
    
    if args.plot_mg:
        #try:
            
        fpr, tpr = get_fpr_tpr(input_ppanini=path + '/stool_ppanini_table.txt',\
                               essential_genes= essential_genes, beta =args.beta)
        roc_info.append(['Stool - metagenomic priority',fpr, tpr])
        
        fpr, tpr = get_fpr_tpr(input_ppanini=path + '/mucosa_ppanini_table.txt',\
                               essential_genes= essential_genes, beta =args.beta)
        roc_info.append(['Buccal mucosa - metagenomic priority',fpr, tpr])
        
        fpr, tpr = get_fpr_tpr(input_ppanini=path + '/fornix_ppanini_table.txt',\
                               essential_genes= essential_genes, beta =args.beta)
        roc_info.append(['Posterior fornix - metagenomic priority',fpr, tpr])
        
        fpr, tpr = get_fpr_tpr(input_ppanini=path + '/nares_ppanini_table.txt',\
                               essential_genes= essential_genes, beta =args.beta)
        roc_info.append(['Anterior nares - metagenomic priority',fpr, tpr])
        
        fpr, tpr = get_fpr_tpr(input_ppanini=path + '/soil_ppanini_table.txt',\
                               essential_genes= essential_genes, beta = args.beta)
        roc_info.append(['Prairie soil - metagenomic priority',fpr, tpr])
        #except ValueError:
         #   print "ValueError for roc calculation"  

    if args.plot_g:
        try:
            # Stool
            metagenomic_table_stool =  utilities.gene2genomes(path +'stool.m8')
            no_uniq_genomes_stool  = utilities.number_of_unique_genomes(path +'/stool.m8')
            fpr, tpr =fpr_tpr_genome(metagenomic_table_stool, no_uniq_genomes_stool, essential_genes)
            roc_info.append(['Stool - genomics priority',fpr, tpr])             
            
            # Buccal mucosa
            metagenomic_table_mucosa =  utilities.gene2genomes(path +'/mucosa.m8')
            no_uniq_genomes_mucosa  = utilities.number_of_unique_genomes(path +'mucosa.m8')
            fpr, tpr =fpr_tpr_genome(metagenomic_table_mucosa, no_uniq_genomes_mucosa, essential_genes)
            roc_info.append(['Buccal mucosa - genomics priority',fpr, tpr])             
            
            # Posterior fornix 
            metagenomic_table_fornix =  utilities.gene2genomes(path +'fornix.m8')
            no_uniq_genomes_fornix  = utilities.number_of_unique_genomes(path +'fornix.m8')
            fpr, tpr =fpr_tpr_genome(metagenomic_table_fornix, no_uniq_genomes_fornix, essential_genes)
            roc_info.append(['Posterior fornix - genomics priority',fpr, tpr]) 
           # Anterior nares
            metagenomic_table_nares =  utilities.gene2genomes(path +'/nares.m8')
            no_uniq_genomes_nares  = utilities.number_of_unique_genomes(path +'/nares.m8')
            fpr, tpr =fpr_tpr_genome(metagenomic_table_nares, no_uniq_genomes_nares, essential_genes)
            roc_info.append(['Anterior nares - genomics priority',fpr, tpr]) 
            
            # Prairie soil
            metagenomic_table_soil =  utilities.gene2genomes(path +'/soil.m8')
            no_uniq_genomes_soil  = utilities.number_of_unique_genomes(path +'/soil.m8')
            fpr, tpr =fpr_tpr_genome(metagenomic_table_soil, no_uniq_genomes_soil, essential_genes)
            roc_info.append(['Prairie soil - genomics priority',fpr, tpr])   
            
        except ValueError:
            print "ValueError for genomics priority  ROC calculation"
            
    try:
        #print roc_info
        roc_plot(roc_info, figure_name=path+ '/'+ args.output+'_ppanini_rocplot', size = args.size) 
    except ValueError:
        print "ValueError for general roc plot"
    print "The ROC evaluation is successfully done!" 

def get_fpr_tpr(input_ppanini, essential_genes, beta =.5):
    fpr = dict()
    tpr = dict()
    true = dict()
    score = dict()
   
    if config.verbose =='DEBUG':
        print "PPANINI ROC plot evaluation!"

    ppanini_table = pd.DataFrame.from_csv(input_ppanini, sep='\t', index_col=0, header =0)
    #print ppanini_table['ppanini_score'][1:10]
    
    #ppanini_table = ppanini_table[ppanini_table['ppanini_score'] > 18.75]
    centroids_list = list(ppanini_table.index  )
    #n = len(centroids_list)-1
    #print config.centroids_list[0:10]
    prev = ppanini_table['prevalence_percentile']
    #prev = [float(val) for val in prev]
    
    #sorted_prev = sorted(prev)
    abun = ppanini_table['abund_percentile']
    #ppanini_score = [1/((1/(beta*prev[i])+(1/((1-beta)* abun[i])))) for i in range(len(prev))] #
    ppanini_score = ppanini_table['ppanini_score']
    abun = [float(val) for val in abun]
    sorted_abun = sorted(abun)
    #print prev[0:101]
    #print abun[0:]
    ground_truth = [1 if (gene_id  in essential_genes) else 0 for gene_id in centroids_list ]
   
    ground_truth = [1 if (gene_id  in essential_genes) else 0 for gene_id in centroids_list ]
    #scipy.stats.rankdata()
    score[beta] = ppanini_score

    
    true[beta] = ground_truth
   
    assert(len(true[beta])==len(score[beta])) 
    fpr[beta], tpr[beta], _  = roc_curve( true[beta], score[beta], pos_label = 1)
    return fpr[beta], tpr[beta]
    
                                     
def roc_plot(roc_info=None, figure_name='roc_plot_ppanini', size= 5, title = ''):
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
    fig, axe = plt.subplots(figsize=(size, size ), dpi=300)#, sharex=False, sharey=False)
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
    axe.set_ylabel('True Positive Rate', fontsize = 8)
    axe.set_xlabel('False Positive Rate', fontsize = 8)
    axe.get_xaxis().set_tick_params(which='both', labelsize=6,top='off',  direction='out')
    axe.get_yaxis().set_tick_params(which='both', labelsize=6, right='off', direction='out')
    axe.yaxis.set_label_position('left') 
    axe.set_title(title, fontsize=9, fontweight='bold')
    # plt.savefig('./test/'+roc_name+'foo.pdf')
    plt.tight_layout()
    #plt.savefig(figure_name+'.pdf')
    plt.savefig(figure_name+'.png')
    #plt.savefig(figure_name+'.svg')


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
    parser.add_argument( "-e1", "--essential-genes1 ",
                         dest = "essentail_gene_file1",
                         metavar = "<feature id>",
                         help="a list of essential genes)", )
    parser.add_argument( "-e2", "--essential-genes2 ",
                         dest = "essentail_gene_file2",
                         metavar = "<feature id>",
                         help="a list of essential genes)", )
    parser.add_argument( "--master-plot", 
                         dest = "master",
                         help="plotting master figure of the paper",
                         action="store_true", )
    parser.add_argument( "--plot-metagenomics", 
                         dest = "plot_mg",
                         default = True,
                         help="plotting based of metagenomics priority score",
                         action="store_true", )
    parser.add_argument( "--plot-genomics", 
                         dest = "plot_g",
                         default = True,
                         help="plotting based of genomics priority score",
                         action="store_true", )
    parser.add_argument( "--path",
                         dest = "path",
                         help="path for inputs and/or outputs", )
    parser.add_argument( "--outfile",
                         dest = "outfile",
                         help="output file", )
    parser.add_argument( "--beta",
                         dest = "beta",
                         default =.5,
                         type =float,
                         help="beta is a weight of contribution, B * prevelence and (1-B) * abundance", )
    parser.add_argument( "--size",
                         dest = "size",
                         default =3,
                         type =int,
                         help="size of the plot in inches", )
    parser.add_argument( "--output",
                         dest = "output",
                         help="a name or output file", )    
    
    return parser.parse_args()
def main():
    user_args = get_args()
    #if user_args.master:
    load_multiple_roc(user_args)
    #else:
    #    priority_scatter(user_args)
        
if __name__ == '__main__':
    main()

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
try:
    from sklearn.model_selection import train_test_split
except:
    #for sklearn older version
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
from ..utilities import load_polymap_dic,load_polymap

from .. import config
#from . import plot_metagenome_genome
def fpr_tpr_genome(input_ppanini, essential_genes, metagenomic_table = None , n_uniq_genomes = None, score_type = 'universe', pan_genome_score = None):
    # calculate genomic score 
    
    ppanini_table = pd.DataFrame.from_csv(input_ppanini, sep='\t', index_col=0, header =0)
    genes = list(ppanini_table.index )
    #['gene']
    #print metagenomic_table.values()[1:100]
    #/float(n_uniq_genomes) 
    #genes = metagenomic_table.keys()
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
    ground_truth = [1 if gene_id  in essential_genes else 0 for gene_id in genes ]
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
    
    '''cog_uniref_map = load_polymap ('/Users/rah/Documents/Hutlab/ppanini_paper/input/essential_gene/keio_essential_uniref90.tsv')
    i= 0
    essential_genes =[]
    for cog in cog_uniref_map:
        i +=1
        if i < 6:
            essential_genes += cog_uniref_map[cog]
        else:
            break'''
    fpr, tpr = get_fpr_tpr(input_ppanini=args.ppanini_output,\
                           essential_genes= essential_genes, beta =args.beta)
    roc_info.append([args.output+' metagenomic priority',fpr, tpr])
     
    pan_genome_score = load_polymap_dic ( '/Users/rah/Documents/Hutlab/ppanini_paper/input/essential_gene/uniphlan90.tsv.gz',
                                      key_ind=0, value_ind = 2)
    metagenomic_table_soil =  utilities.gene2genomes(path + '/'+args.niche+'.txt')
    n_uniq_genomes_soil  = utilities.number_of_unique_genomes(path + '/'+args.niche+'.txt')
    fpr, tpr =fpr_tpr_genome(metagenomic_table_soil, n_uniq_genomes_soil, essential_genes, pan_genome_score = pan_genome_score)
    roc_info.append(['Prairie soil - genomics priority',fpr, tpr])         
    try:
        roc_plot(roc_info, figure_name=args.outfile, size = args.size) 
    except ValueError:
        print "ValueError for general roc plot"
    print "The ROC evaluation is successfully done!" 

def get_fpr_tpr(input_ppanini, essential_genes, beta =.5):
    fpr = dict()
    tpr = dict()
    true = dict()
    score = dict()
   
    #if config.verbose =='DEBUG':
    #    print "PPANINI ROC plot evaluation!"

    ppanini_table = pd.DataFrame.from_csv(input_ppanini, sep='\t', index_col=0, header =0)
    #print ppanini_table['ppanini_score'][1:10]
    
    #ppanini_table = ppanini_table[ppanini_table['ppanini_score'] > 18.75]
    gene_families = list(ppanini_table.index  )
    #n = len(gene_families)-1
    #print config.gene_families[0:10]
    prev = ppanini_table['prevalence_percentile']
    #prev = [float(val) for val in prev]
    
    #sorted_prev = sorted(prev)
    abun = ppanini_table['abund_percentile']
       
    if beta == 0.5:
        ppanini_score = ppanini_table['ppanini_score']
    else:
        #ppanini_score = 1.0/((beta/ppanini_table['prevalence_percentile'])+((1.0-beta)/ppanini_table['abund_percentile']))
        ppanini_score = ppanini_table[['prevalence_percentile','abund_percentile']].max(axis=1)
    abun = [float(val) for val in abun]
    sorted_abun = sorted(abun)
    ground_truth = [1 if gene_id  in essential_genes else 0 for gene_id in gene_families ]
    score[beta] = ppanini_score

    
    true[beta] = ground_truth
   
    assert(len(true[beta])==len(score[beta])) 
    fpr[beta], tpr[beta], _  = roc_curve( true[beta], score[beta], pos_label = 1)
    return fpr[beta], tpr[beta]
    
                                     
def roc_plot(roc_info=None, figure_name='roc_plot_ppanini', size= 5, title = '', axe = None):
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
    save_flag = False
    if axe == None:
        save_flag = True
        fig, axe = plt.subplots(figsize=(size, size ), dpi=300)#, sharex=False, sharey=False)
    #fig.set_size_inches(1, 10)
    #plt.figure(dpi= 300, figsize=(4, 4))
    for i in range(len(roc_info)):
        params = {'legend.fontsize': 3.5
        }
        #mpl.rcParams['lines.linewidth'] = 2
        plt.rcParams.update(params)
        axe.plot(fpr[roc_info[i][0]], tpr[roc_info[i][0]], linewidth=1.0, label='{0} (area = {1:0.2f})'
                                       ''.format(str(roc_info[i][0]), roc_auc[roc_info[i][0]]))   
    axe.plot([0, 1], [0, 1], 'k--', linewidth=0.5)
    axe.set_xlim([-0.01, 1.0])
    axe.set_ylim([-0.01, 1.01])
    axe.legend(loc="lower right")
    axe.set_ylabel('True Positive Rate', fontsize = 6)
    axe.set_xlabel('False Positive Rate', fontsize = 6)
    axe.get_xaxis().set_tick_params(which='both', labelsize=4,top='off',  direction='out')
    axe.get_yaxis().set_tick_params(which='both', labelsize=4, right='off', direction='out')
    axe.yaxis.set_label_position('left') 
    #axe.set_title(title, fontsize=8, fontweight='bold', pos = 'left')
    if title:
        axe.set_title(title, loc='left', fontdict={'fontsize':'8','fontweight' :'bold'})
    # plt.savefig('./test/'+roc_name+'foo.pdf')
    plt.tight_layout()
    #plt.savefig(figure_name+'.pdf')
    #if save_flag:
    plt.savefig(figure_name+'.png')
    #plt.savefig(figure_name+'.svg')
    axe.autoscale_view('tight')
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
    parser.add_argument( "--niche",
                         dest = "niche",
                         help="niche name", )  
    
    return parser.parse_args()
def main():
    user_args = get_args()
    #if user_args.master:
    load_multiple_roc(user_args)
    #else:
    #    priority_scatter(user_args)
        
if __name__ == '__main__':
    main()
import os
import sys
import re
import numpy
import time
import argparse
import pdb
import pandas as pd

import matplotlib
matplotlib.use( "Agg" )
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.patches as patches
import numpy as np
import scipy.cluster.hierarchy as sch

from matplotlib import font_manager
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("error")
    try:
        font_file = font_manager.findfont(font_manager.FontProperties(family='Arial'))
        matplotlib.rcParams["font.family"] = "Arial"
    except UserWarning:
        pass 

def get_args( ):
    parser = argparse.ArgumentParser(
        description="PPANINI plotting tool",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument( "-i1", "--ppanini-input",
                         metavar = "<input table>",
                         dest = "ppanini_input",
                         required = False, 
                         help="Gene abundance table", )
    parser.add_argument( "-i2", "--ppanini-output",
                         metavar = "<input table>",
                         dest = "ppanini_output",
                         required = False, 
                         help="PPANINI output table", )
    parser.add_argument( "--summary-table",
                         dest = "summary_table",
                         #action="store_true",
                         help="Summary table", )
    parser.add_argument( "--scale",
                         dest = "scale",
                         #action="store_true",
                         choices= ['abundance', 'counts'],
                         default = 'abundance',
                         help="Scale: abundance or counts [default: abundance]", )
    parser.add_argument( "-o", "--output",
                         dest = "plot_output",
                         metavar = "<feature id>",
                         help="output plot)", )
    
    return parser.parse_args()
def load_summary_table(summary_table_file_path):
    try: 
        df_in = pd.read_csv(str(summary_table_file_path), sep='\t', header=0, index_col =0)
    except ImportError:
        sys.exit("Input Error First File!") 
    return df_in

def summerize_gene_table(ppanini_input, ppanini_output, scale = 'abundance', output_path = None ):
    try: 
        df_in = pd.read_csv(str(ppanini_input), sep='\t', header=0, index_col =0)
        df_out = pd.read_csv(str(ppanini_output), sep='\t', header=0, index_col =0)
    except ImportError:
        sys.exit("Input Error First File!") 
    
    # add GO to the names of gene families 
    # that have a GO term in the PPANINI output file 
    mapper = []
    # check if the gene family has mapped to any GO term using PPANINI output file
    # add category names in the beginning of gene families
    for gene_family in df_in.index:
        go_term = df_out.loc[gene_family, 'GO']
        
        if go_term !=  go_term: # if mapper[gene_family] is nan then is not equal to itself :)
            if gene_family.startswith('UniRef') and not gene_family.endswith('unknown'):
                mapper.append('Protein_'+gene_family)
            elif gene_family == 'UNMAPPED': 
                mapper.append('Unannotated_'+gene_family)
            else:
                mapper.append('Unannotated_'+gene_family)
        else:
            mapper.append('Function_' + gene_family)
    # update gene family names 
    df_in.index = mapper        
    summary_table = pd.DataFrame(index = df_in.columns, columns=['Unannotated', 'UniRef', 'GO'], dtype=float)
    summary_table[:] = 0.0000
    sum1 =0.0
    for sample in list(df_in.columns):
        
        if scale == 'abundance':
            sum1 = sum(df_in[sample])
            # skip sample with all rows zero
            if sum1 == 0.0:
                continue
            summary_table.loc[sample, 'GO'] = df_in.loc[df_in.index.str.startswith('Function_'),sample].sum() * 1.00 /sum1 
            summary_table.loc[sample, 'UniRef'] = df_in.loc[df_in.index.str.startswith('Protein_'), sample].sum() * 1.00 /sum1 
            summary_table.loc[sample, 'Unannotated'] = df_in.loc[df_in.index.str.startswith('Unannotated_'), sample].sum() * 1.00/sum1
        elif scale == 'counts':
            sum1 = sum(df_in[sample]>0.00)
            # skip sample with all rows zero
            #print sum1
            if sum1 == 0.0:
                continue
            summary_table.loc[sample, 'GO'] = sum(df_in.loc[df_in.index.str.startswith('Function_'),sample]> 0.0) * 1.00 /sum1 
            summary_table.loc[sample, 'UniRef'] = sum(df_in.loc[df_in.index.str.startswith('Protein_'), sample]> 0.0) * 1.00 /sum1 
            summary_table.loc[sample, 'Unannotated'] = sum(df_in.loc[df_in.index.str.startswith('Unannotated_'), sample]>0.0) * 1.00/sum1
            
    summary_table.sort_values(by=['Unannotated', 'UniRef', 'GO' ], ascending=[0, 0 , 0], inplace= True) #
    
    if output_path == None:
        output_path = './ppanini_barplot'
    summary_table.to_csv(output_path+'.txt', sep='\t')
    return summary_table

def stack_barplot(df, output_path = None, axe = None, legend = True, legend_title = "Characterization:", title= None):
       
    # Create the general blog and the "subplots" i.e. the bars
    if axe == None:
        fig, axe = plt.subplots(1, figsize=(6,4))


    if title:
        axe.set_title(title, loc='left', fontdict={'fontsize':'10','fontweight' :'bold'})
    # Set the bar width
    bar_width = 1.0#0.75
    
    # positions of the left bar-boundaries
    bar_l = [i+1 for i in range(len(df.index))]
    
    # positions of the x-axis ticks (center of the bars as bar labels)
    tick_pos = [i+(bar_width/2) for i in bar_l]
    
    # Create a bar plot, in position bar_1
    axe.bar(bar_l,
            df['Unannotated'],
            # set the width
            width=bar_width,
            # with pre_score on the bottom
            #bottom=df['Unannotated'],
            # with the label mid score
            label='Novel protein', #hypothetical unannotated',
            # with alpha 0.5
            alpha=0.9,
            # with color
            color='gold',#'silver',#'#F1911E',
            linewidth=0)
    
    # Create a bar plot, in position bar_1
    axe.bar(bar_l,
            df['UniRef'],
            # set the width
            width=bar_width,
            # with the label pre score
            label='Homologous protein in UniRef',
            bottom= df['Unannotated'],
            # with alpha 0.5
            alpha=0.9,
            # with color
            color= 'blue',#'#29CF66',
            linewidth=0)
    
    # Create a bar plot, in position bar_1
    axe.bar(bar_l,
            # using the pre_score data
            df['GO'],
            # set the width
            width=bar_width,
            # with the label pre score
            label='Functionally known in GO',
            bottom=[i+j for i,j in zip(df['Unannotated'],df['UniRef'])],
            # with alpha 0.5
            alpha=0.9,
            # with color
            color='limegreen',#'#29CE66',
            linewidth=0)
   
    # set the x ticks with names
    plt.xticks(tick_pos, df.index)
    
    # Set the label and legends
    # Put a legend below current axis
    if legend:
        lgd = axe.legend(loc='upper left', bbox_to_anchor=(1, 1),
              fancybox=False, shadow=False,  frameon=False) #ncol=5
        if  legend_title != None:
            lgd.set_title(title = legend_title, prop={'weight':'bold'} )

    #lgd.get_title().set_ha('center')
    #lgd.get_title().set_position((20, 0))
    axe.set_ylabel("Fraction of abundance")
    axe.set_xlabel("Samples (N=%d)" % (len(df.index)))
    axe.get_xaxis().set_tick_params(which='both', labelsize=8,top='off', labelbottom='off', bottom= 'off', direction='out')
    axe.get_yaxis().set_tick_params(which='both', labelsize=8, right='off', direction='out')
    
    params = {'legend.fontsize': 8}#, 'legend.linewidth': 0}
    plt.rcParams.update(params)
    #plt.legend(loc='best', frameon=False )
    plt.xlim([min(tick_pos)-bar_width, max(tick_pos)])
    plt.ylim(0, 1)
    axe.autoscale_view('tight')
    plt.tight_layout()
    #fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
    #plt.show()
    # Set a buffer around the edge
    if output_path == None:
        output_path = './ppanini_barplot'
    plt.savefig(output_path+".pdf", dpi=300, format='pdf', bbox_inches='tight', pad_inches = 0) #dpi=300, format='png', bbox_extra_artists=(lgd,),
    plt.savefig(output_path+".png", dpi=300, format='png', bbox_inches='tight', pad_inches = 0)
    plt.savefig(output_path+".svgz", dpi=300, format='svgz', bbox_inches='tight', pad_inches = 0)
    return axe
def main():
    user_args = get_args()
    if user_args.summary_table:
        df = load_summary_table(user_args.summary_table)
    else:
        df = summerize_gene_table(ppanini_input = user_args.ppanini_input, ppanini_output = user_args.ppanini_output, scale = user_args.scale, output_path = user_args.plot_output)
    stack_barplot(df, user_args.plot_output )
if __name__ == '__main__':
    main()    

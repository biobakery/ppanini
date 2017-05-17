import os
import sys
import re
import numpy
import time
import argparse
import pdb
import pandas as pd
try:
    import matplotlib
    matplotlib.use( "Agg" )
    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["font.family"] = "Arial"
    import matplotlib.pyplot as plt
    from matplotlib import colors
    import matplotlib.patches as patches
    import numpy as np
    import scipy.cluster.hierarchy as sch
except:
    sys.exit( "This script requires the Python scientific stack: numpy, scipy, and matplotlib." )
  
#from . import utils

#from .. import utilities

def get_args( ):
    parser = argparse.ArgumentParser(
        description="HUMAnN2 plotting tool",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument( "-i1", "--ppanini-input",
                         metavar = "<input table>",
                         dest = "ppanini_input",
                         required = False, 
                         help="Gene abundance table", )
    parser.add_argument( "-i2", "--ppanini-output",
                         dest = "ppanini_output",
                         metavar = "<feature id>",
                         help="PPANINI output table)", )
    parser.add_argument( "-i", "--summary-table",
                         metavar = "<input table>",
                         dest = "summary_table",
                         required = False, 
                         help="Summary table", )
    parser.add_argument( "-o", "--output",
                         dest = "plot_output",
                         metavar = "<feature id>",
                         help="output plot)", )
    
    return parser.parse_args()
def load_summary_table(summary_table_file_path):
    try: 
        df_in = pd.read_csv(str(summary_table_file_path), sep='\t', header=0, index_col =0)#, nrows= 1, header=0,) #('/Users/rah/Documents/HMP/metadata/performance_modified.xlsx')
    except ImportError:
        sys.exit("Input Error First File!") 
    return df_in
def summerize_gene_table(ppanini_input_file_path, ppanini_output_file_path, output_path = None ):
     
    try: 
        df_in = pd.read_csv(str(ppanini_input_file_path), sep='\t', header=0, index_col =0)#, nrows= 1, header=0,) #('/Users/rah/Documents/HMP/metadata/performance_modified.xlsx')
        df_out = pd.read_csv(str(ppanini_output_file_path), sep='\t', header=0, index_col =0)
    except ImportError:
        sys.exit("Input Error First File!") 

    
    # summary table for annotations to be used for stacke bar plot (samples * characterization type
    summary_table = pd.DataFrame(index = df_in.columns, columns=['Unannotated', 'UniRef', 'GO'], dtype=float)
    summary_table[:] = 0.0000
    #print summary_table
    for sample in range(len(df_in.columns)):
        
        sum1 = sum(df_in[df_in.columns[sample]])
    
        # skip sample with all rows zero
        if sum1 == 0.00:
            continue
        
        for gene in range(len(df_in.index)):
            # check if nan
            if df_out.iat[gene,df_out.columns.get_loc('GO')] == df_out.iat[gene,df_out.columns.get_loc('GO')]:
                #print (df_out.loc[gene,'GO_term'])
                summary_table.iat[sample, 2] += df_in.iat[gene, sample]
            elif df_in.index[gene].find('UniRef')==0  :
                #print (df_out.loc[gene,'UniRef'])
                summary_table.iat[sample,1] += df_in.iat[gene, sample]
            else:
                summary_table.iat[sample,0] += df_in.iat[gene, sample]
        #print df
        summary_table.iat[sample, 0] /= sum1 #* 100
        summary_table.iat[sample, 1] /= sum1 #* 100
        summary_table.iat[sample, 2] /= sum1 #* 100         
    #print summary_table['Unannotated']
    
    #df = pd.DataFrame(summary_table, columns = summary_table.keys())#['Sample', 'EC', 'GO', 'UniRef', 'Unannotated', 'Unmaped'])
    summary_table.sort_values(by=['Unannotated', 'UniRef', 'GO' ], ascending=[0, 0 , 0], inplace= True) #
    #print summary_table
    if output_path == None:
        output_path = './ppanin_barplot'
    summary_table.to_csv(output_path+'.txt', sep='\t')
    return summary_table
def stache_barplot(df, output_path = None):
       
    # Create the general blog and the "subplots" i.e. the bars
    fig, ax1 = plt.subplots(1, figsize=(6,4))
    title =''# "Charactraztion of genes in NICHE"
    if title:
        ax1.set_title(title, loc='left', fontdict={'fontsize':'10','fontweight' :'bold'})
    # Set the bar width
    bar_width = 1.0#0.75
    
    # positions of the left bar-boundaries
    bar_l = [i+1 for i in range(len(df.index))]
    
    # positions of the x-axis ticks (center of the bars as bar labels)
    tick_pos = [i+(bar_width/2) for i in bar_l]
    
    # Create a bar plot, in position bar_1
    ax1.bar(bar_l,
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
    ax1.bar(bar_l,
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
    ax1.bar(bar_l,
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
    lgd = ax1.legend(loc='upper left', bbox_to_anchor=(1, 1),
          fancybox=False, shadow=False,  frameon=False) #ncol=5
    lgd.set_title(title = "Characterization:", prop={'weight':'bold'} )
    #lgd.get_title().set_ha('center')
    #lgd.get_title().set_position((20, 0))
    ax1.set_ylabel("Fraction of abundance")
    ax1.set_xlabel("Samples (N=%d)" % (len(df.index)))
    ax1.get_xaxis().set_tick_params(which='both', labelsize=8,top='off', labelbottom='off', bottom= 'off', direction='out')
    ax1.get_yaxis().set_tick_params(which='both', labelsize=8, right='off', direction='out')
    ax1.autoscale_view('tight')
    params = {'legend.fontsize': 8}#, 'legend.linewidth': 0}
    plt.rcParams.update(params)
    #plt.legend(loc='best', frameon=False )
    plt.xlim([min(tick_pos)-bar_width, max(tick_pos)])
    plt.ylim(0, 1)
    plt.tight_layout()
    #fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
    #plt.show()
    # Set a buffer around the edge
    if output_path == None:
        output_path = './ppanini_barplot'
    plt.savefig(output_path+".pdf", dpi=300, format='pdf', bbox_inches='tight', pad_inches = 0) #dpi=300, format='png', bbox_extra_artists=(lgd,),
    plt.savefig(output_path+".png", dpi=300, format='png', bbox_inches='tight', pad_inches = 0)
    plt.savefig(output_path+".svgz", dpi=300, format='svgz', bbox_inches='tight', pad_inches = 0)

    #fig.savefig("stacked_barplot.pdf")
def main():
    user_args = get_args()
    if user_args.summary_table:
        df = load_summary_table(user_args.summary_table)
    else:
        df = summerize_gene_table(user_args.ppanini_input, user_args.ppanini_output, user_args.plot_output)
    stache_barplot(df, user_args.plot_output )
if __name__ == '__main__':
    main()    

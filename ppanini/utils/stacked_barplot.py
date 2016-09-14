import os
import sys
import matplotlib
import re
import numpy
import time
import argparse
import pdb
import matplotlib.pyplot as plt
from matplotlib import colors
import pandas as pd
import numpy as np
#from . import utils

#from .. import utilities

def read_gene_table(gene_table_fname = None):
     
    try: 
        df = pd.read_csv(str(sys.argv[1]), sep='\t', header=0)# index_col =0)#, nrows= 1, header=0,) #('/Users/rah/Documents/HMP/metadata/performance_modified.xlsx')
    except ImportError:
        sys.exit("Input Error First File!") 
    #print df['Gene'][1]
    #df
    #print "len rows: ", df.index
    #print "len columns: ", df.columns
    df.rename(columns={'# Gene Family':'Gene'}, inplace=True)
    #print "dataset:\n", df
    #df = df.transpose()
    raw_data = {'Sample': [],
            'EC': [],
            'GO': [],
            'UniRef': [],
            'Unannotated': [],
            'Unmaped': []}
    #anotated = numpy.zeros(len(df.columns)-1)
    for i in range(1,len(df.columns)):
        raw_data['Sample'].append(i)
        raw_data['UniRef'].append(0)
        raw_data['EC'].append(0)
        raw_data['GO'].append(0)
        raw_data['Unannotated'].append(0)
        raw_data['Unmaped'].append(0)
        sum1 = sum(df[df.columns[i]][:])
        #print sum1
        for j in range(len(df.index)):
            if df['Gene'][j].find('EC')==0:
                raw_data['EC'][i-1] +=df[df.columns[i]][j]
            elif df['Gene'][j].find('GO')==0:
                raw_data['GO'][i-1] +=df[df.columns[i]][j]
            elif df['Gene'][j].find('UniRef')==0:
                raw_data['UniRef'][i-1] +=df[df.columns[i]][j]
            elif df['Gene'][j].find('UNMAPPED')==0:
                raw_data['Unmaped'][i-1] +=df[df.columns[i]][j]
            else:
                raw_data['Unannotated'][i-1] +=df[df.columns[i]][j]
        raw_data['UniRef'][i-1] = raw_data['UniRef'][i-1]/sum1
        raw_data['Unannotated'][i-1] = raw_data['Unannotated'][i-1]/sum1
        raw_data['Unmaped'][i-1] = raw_data['Unmaped'][i-1]/sum1
        raw_data['EC'][i-1] = raw_data['EC'][i-1]/sum1
        raw_data['GO'][i-1] = raw_data['GO'][i-1]/sum1         
    #print raw_data
    
    df = pd.DataFrame(raw_data, columns = ['Sample', 'EC', 'GO', 'UniRef', 'Unannotated', 'Unmaped'])
    df = df.sort('UniRef', ascending=0)
    return df
def stache_barplot(df):
       
    # Create the general blog and the "subplots" i.e. the bars
    f, ax1 = plt.subplots(1, figsize=(5,5))
    
    # Set the bar width
    bar_width = 0.75
    
    # positions of the left bar-boundaries
    bar_l = [i+1 for i in range(len(df['UniRef']))]
    
    # positions of the x-axis ticks (center of the bars as bar labels)
    tick_pos = [i+(bar_width/2) for i in bar_l]
    # Create a bar plot, in position bar_1
    ax1.bar(bar_l,
            # using the pre_score data
            df['EC'],
            # set the width
            #width=bar_width,
            # with the label pre score
            label='EC',
            # with alpha 0.5
            alpha=0.5,
            # with color
            color='#29CC66',
            linewidth=0)
    # Create a bar plot, in position bar_1
    ax1.bar(bar_l,
            # using the pre_score data
            df['GO'],
            # set the width
            #width=bar_width,
            # with the label pre score
            label='GO',
            # with alpha 0.5
            alpha=0.5,
            # with color
            color='#29CE66',
            linewidth=0)
    # Create a bar plot, in position bar_1
    ax1.bar(bar_l,
            # using the pre_score data
            df['UniRef'],
            # set the width
            #width=bar_width,
            # with the label pre score
            label='UniRef',
            # with alpha 0.5
            alpha=0.5,
            # with color
            color='#29CF66',
            linewidth=0)
    
    # Create a bar plot, in position bar_1
    ax1.bar(bar_l,
            # using the mid_score data
            df['Unannotated'],
            # set the width
            #width=bar_width,
            # with pre_score on the bottom
            bottom=df['Unannotated'],
            # with the label mid score
            label='Unannotated',
            # with alpha 0.5
            alpha=0.5,
            # with color
            color='#F1911E',
            linewidth=0)
    
    # Create a bar plot, in position bar_1
    ax1.bar(bar_l,
            # using the post_score data
            df['Unmaped'],
            # set the width
            #width=bar_width,
            # with pre_score and mid_score on the bottom
            bottom=[i+j for i,j in zip(df['UniRef'],df['Unannotated'])],
            # with the label post score
            label='Unmaped',
            # with alpha 0.5
            alpha=0.5,
            # with color
            color='#F1BD1A',
            linewidth=0)
    
    # set the x ticks with names
    #plt.xticks(tick_pos, df['Sample'])
    
    # Set the label and legends
    ax1.set_ylabel("Normalized Abundances")
    ax1.set_xlabel("Samples")
    ax1.get_xaxis().set_tick_params(which='both', labelsize=8,top='off', labelbottom='off', bottom= 'off', direction='out')
    ax1.get_yaxis().set_tick_params(which='both', labelsize=8, right='off', direction='out')
    params = {'legend.fontsize': 8,
        'legend.linewidth': 0}
    plt.rcParams.update(params)
    plt.legend(loc='best', frameon=False )
    plt.tight_layout()
    # Set a buffer around the edge
    #plt.xlim([min(tick_pos)-bar_width, max(tick_pos)+bar_width])
    plt.savefig("stacked_barplot.pdf")

if __name__ == '__main__':
    df = read_gene_table()
    stache_barplot(df)
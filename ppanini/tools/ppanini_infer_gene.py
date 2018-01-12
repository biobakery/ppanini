#!/usr/bin/env python

"""
Infer from abundances

This module will distinguish between genes with sufficient uniref90 map. 


"""
import os
import sys
import re
import argparse
import csv
import shutil
import json
import pandas as pd
from ..utilities import make_directory, iter_rows

def get_args():
    parser = argparse.ArgumentParser()
    #parser.add_argument( "orfs" )
    parser.add_argument( "hits" )
    parser.add_argument( "--output", required=True )
    parser.add_argument( "--min-percid", type=float, required=True )
    parser.add_argument( "--min-qcover", type=float, required=True )
    parser.add_argument( "--min-scover", type=float, required=True )
    parser.add_argument( "--all-valid-hits", action="store_true" )
    args = parser.parse_args()
    return args
def main():
    args = get_args()
    c_min_percid = args.min_percid
    c_min_qcover = args.min_qcover
    c_min_scover = args.min_scover
    
    """
    >NODE_1_length_628004_cov_22.174_1 # 2 # ...
    """
    
    
    
    """
    00 NODE_1_length_628004_cov_22.174_1
    01 UniRef90_R7JM20|1245
    02 94.3
    03 193
    04 11
    05 0
    06 1
    07 193
    08 141
    09 333
    10 4.9e-99
    11 367.1
    12 194
    13 415
    """
    abunds = {}
    # make the output directory
    make_directory(args.output)
    unirefs = {}
    no_unirefs = {}
    polymap_all = {}
    for row in iter_rows( args.hits ):
        gene = row[0]
        uniref = row[1].split( "|" )[0]
        percid = float( row[2] ) / 100.0
        qcover = (int( row[7] ) - int( row[6] ) + 1) / float( row[12] )
        scover = (int( row[9] ) - int( row[8] ) + 1) / float( row[13] )
        if percid >= c_min_percid \
                and qcover >= c_min_qcover \
                and scover >= c_min_scover:
            if gene not in unirefs or args.all_valid_hits:
                polymap_all.setdefault( uniref, {} )[gene] = 1 
                unirefs.setdefault( gene, set( ) ).add( uniref )
                
    # use gene-uniref map for genes that pass the thresholds
    #f1=open(args.output+'/map_uniref_gene.txt',"wt")
    #for uniref in unirefs:
    #    f1.write ("%s\t%s\n" % (uniref,str(';'.join(gene for gene in polymap_all.get(uniref)))))
    #
    # use gene-uniref map for genes that pass the thresholds
    with open(args.output+'/uniref_gene_map.txt', 'wt') as csv_file:
            writer = csv.writer(csv_file, delimiter='\t')
            for gene in unirefs:
               writer.writerow([list(unirefs[gene])[0]], gene)
    
    # make dictionary for genes that don't pass the threshold against uniref    
    for row in iter_rows( args.hits ):
        gene = row[0]
        if gene not in unirefs:
            no_unirefs.setdefault( gene, set( ) ).add( uniref )
    
    # use gene-uniref map for genes that pass the thresholds
    with open(args.output+'/hits.txt', 'wt') as csv_file:
            writer = csv.writer(csv_file, delimiter='\t')
            for gene in unirefs:
               writer.writerow([gene])
    # list of genes that doesn't pass the threshold for mapping to uniref  
    with open(args.output+'/no_hits.txt', 'wt') as csv_file:
            writer = csv.writer(csv_file, delimiter='\t')
            for gene in no_unirefs:
               writer.writerow([gene])  

if __name__ == '__main__':
    main()

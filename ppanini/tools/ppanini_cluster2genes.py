#!/usr/bin/env python
"""
This script make a map table for CD-Hit output
Cluster\tList of genes separated by ;\t Representative gene 
"""

import os
import sys
import re
import argparse
import csv
import shutil
import json
from ..utilities import make_directory

def  get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument( 
        '-i', '--cd-hit', dest ="cd_hit" ,
        help='CD-HIT clusters output', 
        )
    parser.add_argument( "-o", "--output", default = 'CD-Hit-Map', required=True )
    parser.add_argument( "--json", action="store_true" )
    
    args = parser.parse_args()
    return args

def main():
    
    args = get_args()
    
    # make the output directory
    make_directory(args.output)
    
    abunds = {}
    polymap_all = {}
    with open( args.cd_hit ) as fh: # it was args.orfs
        for line in fh:
            line = line.rstrip()
            if line[0] == ">":
               cluster = line.split(">")[1]
            else:
                gene = line.split(">")[1].split("...")[0]
                #print gene, line    
                polymap_all.setdefault( cluster, {} )[gene] = 1    
                if line.endswith('*'):
                    polymap_all[cluster]['representaive_gene']= gene
    if args.json:
        f=open(args.output+'/cd-hit_gene_map.txt.json',"wt")
        f.write(json.dumps(polymap_all))
    
    f1=open(args.output+'/cd-hit_gene_map.txt',"wt")
    # cluster_id  cluster_rep genes
    for cluster in polymap_all:
        f1.write ("%s\t%s\t%s\n" % (cluster, str(';'.join(gene for gene in polymap_all.get(cluster) )), polymap_all[cluster]['representaive_gene']))
if __name__ == '__main__':
    main() 

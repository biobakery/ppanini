#!/usr/bin/env python

import os
import sys
import re
#!/usr/bin/env python
"""
This script make a map table for CD-Hit output
Cluster\tList of genes separated by ;\t Representative gene 
"""

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
    parser.add_argument(
        '-f', '--fasta',
        help='fasta file for mapped genes to clusters (uniref90)', 
        )
    parser.add_argument( "--output", default = 'CD-Hit-Map', required=True )
    parser.add_argument( "--json", action="store_true" )
    
    args = parser.parse_args()
    return args

def main():
    
    args = get_args
    
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
                polymap_all.setdefault( cluster, {} )[gene] = 1    
                if line.endswith('*'):
                    polymap_all[cluster]['rep']= gene
    if args.json:
        f=open(args.output+'/json_map_uniref_gene.txt',"wt")
        f.write(json.dumps(polymap_all))
    
    f1=open(args.output+'/map_uniref_gene.txt',"wt")
    # cluster_id  cluster_rep genes
    for cluster in polymap_all:
        f1.write ("%s \t %s \t %s \n" % (cluster, str(';'.join(gene for gene in polymap_all.get(cluster))), polymap_all[cluster]['rep']))
    

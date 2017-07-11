
import os
import sys
import re
import argparse
import csv
import shutil
import json
import pandas as pd
from .. import utilities
from utilities import make_directory, iter_rows
parser = argparse.ArgumentParser()
#parser.add_argument( "orfs" )
parser.add_argument( "hits" )
parser.add_argument( "--output", required=True )
parser.add_argument( "--min-percid", type=float, required=True )
parser.add_argument( "--min-qcover", type=float, required=True )
parser.add_argument( "--min-scover", type=float, required=True )
parser.add_argument( "--all-valid-hits", action="store_true" )
parser.add_argument( "--json", action="store_true" )
args = parser.parse_args()

c_min_percid = args.min_percid
c_min_qcover = args.min_qcover
c_min_scover = args.min_scover

"""
>NODE_1_length_628004_cov_22.174_1 # 2 # ...
"""

abunds = {}
with open( args.hits ) as fh: # it was args.orfs
    for line in fh:
        if line[0] == ">":
           gene = line[1:].split( )[0]
           #abunds[gene] = float( gene.split( "_" )[5] )

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
            
if args.json:
    f=open(args.output+'/json_map_uniref_gene.txt',"wt")
    f.write(json.dumps(polymap_all))

f1=open(args.output+'/map_uniref_gene.txt',"wt")
for cluster in polymap_all:
    f1.write ("%s \t %s \n" % (cluster,str(';'.join(gene for gene in polymap_all.get(cluster)))))
        #f1.write([cluster, polymap_all.get(cluster)])

for row in iter_rows( args.hits ):
    gene = row[0]
    if gene not in unirefs:
        no_unirefs.setdefault( gene, set( ) ).add( uniref )

with open(args.output+'/hits.txt', 'wt') as csv_file:
        writer = csv.writer(csv_file, delimiter='\t')
        for gene in unirefs:
           writer.writerow([gene])  
with open(args.output+'/no_hits.txt', 'wt') as csv_file:
        writer = csv.writer(csv_file, delimiter='\t')
        for gene in no_unirefs:
           writer.writerow([gene])  
  
'''output = {}
for gene in abunds:
    if len( unirefs.get( gene, [] ) ) > 1:
        warn( unirefs[gene] )
    for u in unirefs.get( gene, ["UniRef90_unknown"] ):
        output[u] = output.get( u, 0 ) + abunds[gene]

tprint( "# GENE", "ABUND" )
for u in sorted( output, key=lambda x: -output[x] ):
    tprint( u, output[u] )

a = len( abunds )
u = len( unirefs )
warn( a, u, u / float( a ) )'''

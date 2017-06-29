#!/usr/bin/env python

import os
import sys
import re
import argparse
import csv
import shutil
#from zopy.utils import iter_rows, tprint, warn

#==== Helper function from zopy by Eric Franzosa ========
def warn ( *args ):
    script = "?"
    if sys.argv[0] != "":
        script = os.path.split( sys.argv[0] )[1].upper()
    args = ["WARNING ({}):".format( script )] + list( args )
    print >>sys.stderr, " ".join( map( str, args ) )

def iter_rows( path ):
    """ easy table loading """
    lens = []
    with try_open( path ) as fh:
        for row in reader( fh ):
            lens.append( len( row ) )
            yield row
    if len( set( lens ) ) != 1:
        warn( "rows didn't all have equal lengths:", set( lens ) )

def tprint( *args, **kwargs ):
    """ coerce list of items to strings then print with tabs between """
    print >>kwargs.get( "file", sys.stdout ), "\t".join( map( str, args ) )

def try_open( path, *args ):
    """ open an uncompressed or gzipped file; fail gracefully """
    fh = None
    try:
        if re.search( r".gz$", path ):
            print >>sys.stderr, "Treating", path, "as gzipped file"
            fh = gzip.GzipFile( path, *args )
        else:
            fh = open( path, *args )
    except:
        die( "Problem opening", path )
    return fh   
# ---------------------------------------------------------------
# text manipulation
# ---------------------------------------------------------------

def reader ( file_handle ):
    """ my favorite options for csv reader """
    for aItems in csv.reader( file_handle, delimiter="\t", quotechar="", quoting=csv.QUOTE_NONE ):
        yield aItems

def make_directory(output_dir):
    if not os.path.isdir(output_dir):
        try:
            print("Creating output directory: " + output_dir)
            os.mkdir(output_dir)
        except EnvironmentError:
            sys.exit("CRITICAL ERROR: Unable to create output directory.")
    else:
        try:
            print("Removing the old output directory: " + output_dir)
            shutil.rmtree(output_dir)
            print("Creating output directory: " + output_dir)
            os.mkdir(output_dir)
        except EnvironmentError:
            sys.exit("CRITICAL ERROR: Unable to create output directory.")
        
    
    if not os.access(output_dir, os.W_OK):
        sys.exit("CRITICAL ERROR: The output directory is not " + 
            "writeable. This software needs to write files to this directory.\n" +
            "Please select another directory.")
        
    print("Output files will be written to: " + output_dir) 
parser = argparse.ArgumentParser()
#parser.add_argument( "orfs" )
parser.add_argument( "hits" )
parser.add_argument( "--output", required=True )
parser.add_argument( "--min-percid", type=float, required=True )
parser.add_argument( "--min-qcover", type=float, required=True )
parser.add_argument( "--min-scover", type=float, required=True )
parser.add_argument( "--all-valid-hits", action="store_true" )
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
            unirefs.setdefault( gene, set( ) ).add( uniref )
    else:
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

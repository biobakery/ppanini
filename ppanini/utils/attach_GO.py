from __future__ import print_function # PYTHON 2.7+ REQUIRED
import os
import re
import pdb
import os
import sys
import re
import subprocess
import csv
import gzip
import bz2
import pandas as pd
'''Tmp file to parse results'''
def read_map(map_obj):
	csv_map = csv.reader(open(map_obj), csv.excel_tab)
	uniref_go = {}
	for line in csv_map:
		for uid in line[1:]:
			uniref_go[uid.strip()] = line[0]
	return uniref_go

# constants
# ---------------------------------------------------------------

c_strat_delim     = "|"
c_taxon_delim     = "."
c_name_delim      = ": "
c_multiname_delim = ";"
c_str_unknown     = "NO_NAME"
c_ungrouped       = "UNGROUPED"
c_unmapped        = "UNMAPPED"
c_unintegrated    = "UNINTEGRATED"
c_many_bytes      = 1e8
c_zip_multiplier  = 10
# the last line in the file with this indicator is the header
GENE_TABLE_COMMENT_LINE="#"

# the extension used for biom files
BIOM_FILE_EXTENSION=".biom"

c_topsort = {
    c_unmapped:0,
    c_ungrouped:1,
    c_unintegrated:2,
    "UniRef50_unknown":3,
    "UniRef90_unknown":4,
}

# ---------------------------------------------------------------
# helper functions
# ---------------------------------------------------------------
def size_warn( path ):
    m = 1 if ".gz" not in path else c_zip_multiplier
    if m * os.path.getsize( path ) > c_many_bytes:
		print( "  This is a large file, one moment please...", file=sys.stderr )

def try_zip_open( path, write=None ):
    """ 
    open an uncompressed or gzipped file; fail gracefully 
    """
    fh = None

    # set the open mode
    if write:
        open_mode = "w"
    elif path.endswith(".bz2"):
        open_mode = "r"
    else:
        open_mode = "rt"

    try:
        if path.endswith(".gz"):
            fh = gzip.open( path, open_mode )
        elif path.endswith(".bz2"):
            fh = bz2.BZ2File( path, open_mode )
        else:
            fh = open( path, open_mode )
    except EnvironmentError:
        sys.exit( "Problem opening file: " + path)
    return fh

def read_biom_table( path ):
    """
    return the lines in the biom file
    """

    try:
        import biom
    except ImportError:
        sys.exit("Could not find the biom software."+
            " This software is required since the input file is a biom file.")
        
    try:
        tsv_table = biom.load_table( path ).to_tsv().split("\n")
    except (EnvironmentError, TypeError):
        sys.exit("ERROR: Unable to read biom input file.")
        
    return tsv_table

def gzip_bzip2_biom_open_readlines( path ):
    """
    return the lines in the opened file for tab delimited text, gzip, bzip2 and biom files
    """

    # if the file is biom, convert to text and return lines
    if path.endswith(BIOM_FILE_EXTENSION):
        for line in read_biom_table(path):
            yield line
    else:
        with try_zip_open( path ) as file_handle:
            for line in file_handle:
                if path.endswith(".bz2"):
                    # convert the line to text from binary
                    yield line.decode('utf-8').rstrip()
                else:
                    yield line.rstrip()

def load_polymap ( path, start=0, skip=None, allowed_keys=None, allowed_values=None ):
    """
    Load a file like:
    A 1 2
    B 1
    B 3
    C 1 2 4
    To a nested dict structure:
    {A:{1:1, 2:1}, B:{1:1, 3:1}, C:{1:1, 2:2, 4:1}
    Inner values are not important (set to 1)
    """
    polymap = {}
    print( "Loading mapping file from:", path, file=sys.stderr )
    size_warn( path )
    for line in gzip_bzip2_biom_open_readlines( path ):
        row = line.split("\t")
        key = row[start]
        if allowed_keys is None or key in allowed_keys:
            for i, value in enumerate( row ):
                if i != start and (skip is None or i not in skip):
                    if allowed_values is None or value in allowed_values:
                        #polymap.setdefault( value, {} )[key] = 1 #polymap.setdefault( key, {} )[value] = 1
                        polymap[value] = str(key)
                        #if value in polymap.keys() :
                        #	polymap[value] += ("|"+str(key) )
                        #else:
                        #	polymap[value] = str(key)
    df_polymap = pd.DataFrame.from_dict(polymap,orient='index')#, columns =["GO_term"]) 
    df_polymap.to_csv("/Users/rah/Documents/Hutlab/ppanini/ppanini/data/map_uniref90_infogo1000.txt", sep='\t')                                       	
    #print("Mappping Uniref90 to infogo1000 is done")
    #print (polymap)
    return polymap
def uniref2go(ppanini_table, path = "/Users/rah/Documents/Hutlab/ppanini/ppanini/data/map_infogo1000_uniref90.txt.gz" ):
	go1000_uniref90_dic = load_polymap ( path )
	
	#print('Loading the mapping file is done!')
	#print (go1000_uniref90_dic.keys())
	#uniref_go_keys = go1000_uniref90_dic.keys()
	#go1000_uniref90_dic = pd.read_table("file.gz",compression='gzip',sep='\x01')
	for index, row in ppanini_table.iterrows():
	#	row["go_term"] = go1000_uniref90_dic.get(index, None)
		#if index == '':
		#	continue
		#print(go1000_uniref90_dic.get(index))
		ppanini_table.loc[index,'GO'] = go1000_uniref90_dic.get(index)# go1000_uniref90_dic.get(ppanini_table.index)
	#ppanini_table.loc[ppanini_table.index,'go_term'] = go1000_uniref90_dic[list(ppanini_table.index)][1]
	#print(ppanini_table['go_term']) 
	#return 	ppanini_table
	
if __name__ == '__main__':
	uniref_go = read_map(sys.argv[1])
	csv_table = csv.reader(open(sys.argv[2]), csv.excel_tab)
	check = 0
	new_foo = []
	for line in csv_table:
		new_line = ''
		if check:
			uid = line[0]
			if uid.strip() in uniref_go:
				head = [uid+'|'+uniref_go[uid]]
				new_line = head + line[1:]
				new_foo += [new_line]
			else:
				new_foo +=[line]
		else:
			new_foo += [line]
		check +=1
	with open(sys.argv[2]+'_GOadded.txt', 'w') as foo:
		for line in new_foo:
			foo.writelines(['\t'.join(line)+'\n'])



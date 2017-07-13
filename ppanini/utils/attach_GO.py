
import os
import re
import pdb
import os
import sys
import re
import subprocess
import csv
import bz2
import pandas as pd
from .. import utilities

def uniref2go(ppanini_table, uniref_go_path ):
    #go1000_uniref90_dic = load_polymap ( go_uniref_path )
    go1000_uniref90_dic = utilities.load_polymap_dic ( uniref_go_path )
    #print('Loading the mapping file is done!')
    #print (go1000_uniref90_dic.keys())
    #uniref_go_keys = go1000_uniref90_dic.keys()
    #go1000_uniref90_dic = pd.read_table("file.gz",compression='gzip',sep='\x01')
    for index, row in ppanini_table.iterrows():
    	ppanini_table.loc[index,'GO'] = go1000_uniref90_dic.get(index)
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
    
    

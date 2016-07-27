import os
import re
import sys

'''Takes a list of Uniref ids and Uniref_Go mapper file; 
   Extracts all the UniRefIDs with the GO id'''
   
def get_go_mapping(uniref_table, map_go_fname):
	map = open(map_go_fname,'r')
	mapper = {}
	for line in map:
		split_line = line.split('\t')
		for i in split_line[1:]:
			mapper[re.sub('[\r\t\n]','',i)] = split_line[0]
	for gene in uniref_table:
		if gene in mapper:
			print gene+'\t'+mapper[gene]
		else:
			print gene+'\tNA'
def main():
    cmd_h = ['-h', '--help']
    if sys.argv[1] in cmd_h:
        print 'usage: python write_mapper <uniref_ids> <map uniref_go_ids> > <uniref_ids_go_select>'
        sys.exit(0)
    input_table = sys.argv[1]
    uniref_table = [re.sub('[\r\t\n]','', line).strip() for line in open(input_table)]
    map_go_fname = sys.argv[2]
    get_go_mapping(uniref_table, map_go_fname)
    
if __name__ == '__main__':
	main()

	
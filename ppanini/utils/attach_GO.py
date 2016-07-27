import os
import re
import pdb
import sys
import csv
'''Tmp file to parse results'''
def read_map(map_obj):
	csv_map = csv.reader(open(map_obj), csv.excel_tab)
	uniref_go = {}
	for line in csv_map:
		for uid in line[1:]:
			uniref_go[uid.strip()] = line[0]
	return uniref_go

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


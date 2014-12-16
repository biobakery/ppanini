import os
import sys
import pdb
import re


def clean_make(file_path, tagger):
    '''tagger = faa or fna or txt '''
    files = os.listdir(file_path)
    files = [i for i in files if 'SRS' in i]
    for i in files:
        foo = open(file_path + '/' + i, 'r')
        foo_lines = foo.readlines()
        ind_tag = False
        for l in range(len(foo_lines)):
            if '>' == foo_lines[l][0]:
                ind_tag = l
                break
        if not tagger == 'txt':
            nuc_tag = False
            tags_list = [p for p in range(
                ind_tag, len(foo_lines)) if foo_lines[p][0] == '>'] + [len(foo_lines)]

            nucs = ['A', 'C', 'G', 'T', 'N']
            for m in range(0, len(tags_list) - 1):
                check = True
                sub_str = re.sub(
                    '[\r\n\t]', '', ''.join(foo_lines[tags_list[m] + 1:tags_list[m + 1]]))
                for n in sub_str:
                    if not n in nucs:
                        check = False
                        break
                if check:
                    nuc_tag = tags_list[m]
                    break
        print i + '\tcompleted'
	try:
        	os.mkdir(tagger+'_only')
	except:
		pass
        foo_new = open(tagger + '_only/' + i + '.' + tagger, 'w')
        if tagger == 'txt':
            foo_new.writelines(foo_lines[:ind_tag])
        elif tagger == 'faa':
            foo_new.writelines(foo_lines[ind_tag:nuc_tag])
        elif tagger == 'fna':
            foo_new.writelines(foo_lines[nuc_tag:])
        foo_new.close()
        foo.close()

if __name__ == '__main__':
    if len(sys.argv) == 3:
        clean_make(sys.argv[1], sys.argv[2])
    else:
        print 'Usage python %s <gff3_folder_path> <txt or faa or fna>'

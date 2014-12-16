import os
import re
import sys

if __name__ == '__main__':
	if not len(sys.argv) == 4:
		print 'Usage %s <abundance file> <uniref annotation file> <outputfile>' %sys.argv[0]
		raise IOError 
		
	foo1 = open(sys.argv[2]) #uniref file
	foo1 = foo1.readlines()

	foo2 = open(sys.argv[1]) #abundance file
	foo2 = [i.split('\t')[0:4] for i in foo2.readlines()]

	new_list = []
	for i in foo2:
		new_list += [[re.sub('[\r\n\t]','',j) for j in i]]

	uniref_dict = {}
	for i in foo1:
		split_i = i.split('\t')
		if i.split('\t')[0].split(' ')[0] not in uniref_dict:
			uniref_dict[split_i[0].split(' ')[0].strip()] = {}
			uniref_dict[split_i[0].split(' ')[0].strip()][split_i[1]] = [float(split_i[2])*float(split_i[3])/float(100)]
		else:
			uniref_dict[split_i[0].split(' ')[0].strip()][split_i[1]] = [float(split_i[2])*float(split_i[3])/float(100)]
	short_udict = {}
	for i in uniref_dict:
		if len(uniref_dict[i]) >1:
			uid = False
			for j in uniref_dict[i]:
				if not uid:
					uid = j
				elif uniref_dict[i][j]> uniref_dict[i][uid]:
					uid = j
			short_udict[i] = uid
		else:
			short_udict[i] = uniref_dict[i].keys()[0]
	if not sys.argv[4] = '':
		annot = {}
		for key in short_udict:
			if 'UniRef90' in short_udict[key]:
				fname = open('/n/huttenhower_lab_nobackup/downloads/HMP/HMGI/stool/'+key.split('.')[0]+'.with_fasta.gff3')
				for i in fname.readlines():
					if key in i and 'product' in i:
						if 'GO:' in i or 'EC:' in i:
							annot[key] = 'Characterized'
							break
						else:
							annot[key] = 'Uncharacterized'
							break
		elif 'UniRef50' in short_udict[key]:
			annot[key] = 'Uncertain'
	foo = open(sys.argv[3], 'w')
	for i in new_list:
		if i[0].strip()+'-T1-C' in short_udict:
			foo.writelines([str.join('\t', i)+'\t'+str(short_udict[i[0].strip()+'-T1-C'])+'\t'+annot[i[0].strip()+'-T1-C']+'\n'])
		else:
			foo.writelines([str.join('\t', i)+'\t \tNovel\n'])
	foo.close()



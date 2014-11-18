import os
import re
import sys
import pdb
import time
from Bio import SeqIO

if __name__ == '__main__':
	#python script.py results_folder output_folder
#	files = os.listdir(sys.argv[1])
	foo = sys.argv[1]
	existing_files = os.listdir(sys.argv[2])
	not_parsed = open(sys.argv[3])
	not_parse = [re.sub('[\r\n\t]','',i)+'.m8' for i in not_parsed.readlines()]
	t =time.time()
	fname= foo.rpartition('/')[-1]
	if not fname.strip() in existing_files and not fname.strip() in not_parse:
		print foo
		fname= foo.rpartition('/')[-1]
		foo1 = open(foo)#sys.argv[1]+'/'+foo)
		foo2 = open(sys.argv[2]+'/'+fname, 'w')
		foo1 = foo1.readlines()
		gis_u90 = []
		gis_u50 = []
                handle = open('/n/regal/huttenhower_lab/annot_pipe/repophlan_31122013_speciescentroids_ffn/'+fname[:-3], 'rU')
                record_dict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
                handle.close()
		for line in foo1:
			if not line[0]=='#':
				if ('UniRef90' in line.split('\t')[1] and not line.split('\t')[0] in gis_u90) or ('UniRef50' in line.split('\t')[1] and not line.split('\t')[0] in gis_u50):
					split_i = line.split('\t')
					id = float(split_i[2])/100.0
					aln = float(split_i[3])
#					q_len = split_i[0].rpartition('|')[-1]
#					try:
#						if not ',' in q_len:
#							q_len = q_len.split('-')
#							q_len = (abs(float(q_len[0].split('c')[-1].split(':')[-1])-float(q_len[1]))+1.0)/3.0
#						else:
#							q_len = q_len.split(',')
#							q_len1 = q_len[0].split('-')
#							q_len1 = abs(float(q_len1[0].split('c')[-1].split(':')[-1])-float(q_len1[1]))
#							q_len2 = q_len[1].split('-')
#							q_len = (abs(float(q_len2[0].split('c')[-1].split(':')[-1])-float(q_len2[1]))+q_len1+1.0)/3.0
#					except:
					q_len = float(len(record_dict[split_i[0].strip()].seq))/3.0
					num = id*aln*100.0/q_len
#					gis_uniref += [split_i[0]]
					if 'UniRef50' in split_i[1]:
						gis_u50 += [split_i[0]]
						check = num >= 50.0
					elif 'UniRef90' in split_i[1]:
						gis_u90 += [split_i[0]]
						check = num >=90.0
					if check:
						foo2.writelines([split_i[0]+'\t'+split_i[1]+'\t'+str(num)+'\t'+str(check)+'\n'])
		foo2.close()
		print str((time.time()-t)/60)+' mins elapsed'

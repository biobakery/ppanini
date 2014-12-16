import os
import re
import sys

def get_uniref_data(filename):	
	foo1 = open(filename)#sys.argv[2]) #uniref file
	foo1 = foo1.readlines()

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
	return short_udict

def get_abund_data(filename):
	foo2 = open(filename) #sys.argv[1]) #abundance file
	foo2 = [i.split('\t')[0:4] for i in foo2.readlines()]
	new_list = []
	for i in foo2:
		new_list += [[re.sub('[\r\n\t]','',j) for j in i]]
	return new_list
		
def write_annotated_table(new_list, short_udict, filename, annot):
	foo = open(filename, 'w') #sys.argv[3], 'w')
	if annot == []:
		for i in new_list:
			if i[0].strip()+'-T1-C' in short_udict:
					foo.writelines([str.join('\t', i)+'\t'+str(short_udict[i[0].strip()+'-T1-C'])+'\n'])
			elif i[0].strip() in short_udict: 
					foo.writelines([str.join('\t', i)+'\t'+str(short_udict[i[0].strip()])+'\n'])
			else:
				foo.writelines([str.join('\t', i)+'\n'])
	else:
	        for i in new_list:
			if not 'T1-C' in i[0]:
				if i[0].strip()+'-T1-C' in short_udict:
					foo.writelines([str.join('\t', i)+'\t'+str(short_udict[i[0].strip()+'-T1-C'])+'\t'+annot[i[0].strip()+'-T1-C']+'\n'])
	                	else:
        	                	foo.writelines([str.join('\t', i)+'\t \tNovel\n'])
			else:
				if i[0].strip() in short_udict:
        		                foo.writelines([str.join('\t', i)+'\t'+str(short_udict[i[0].strip()+'-T1-C'])+'\t'+annot[i[0].strip()]+'\n'])
	                	else:
        	                	foo.writelines([str.join('\t', i)+'\t \tNovel\n'])
	foo.close()

def get_annotations(short_udict, gi_annot):
	annot = {}
        for key in short_udict:
        	if 'UniRef90' in short_udict[key]:
                	val = gi_annot[key.split('.')[0]][key] 
			if 'GO:' in val or 'EC:' in val:
				annot[key] = 'Characterized'
                        else:   
                                annot[key] = 'Uncharacterized'
                elif 'UniRef50' in short_udict[key]:
                	annot[key] = 'Uncertain'
	return annot

def get_gi_annot(foldername):
	files = os.listdir(foldername)
	gi_annot = {}
	for i in files:
		fname = i.split('.')[0]
		foo = open(foldername+'/'+i)
		foo = foo.readlines()
		gi_annot[fname] = {}
		for j in foo:
			gi_annot[fname][j.split('ID=')[-1].split(';')[0]]= re.sub('[\r\t\n]','',j.split('product=')[-1])
	return gi_annot	

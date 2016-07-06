import sys
import re
import os
import argparse

#help if you select the wrong command line arguments 
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write("python selectGoTerms.py PPANINIOutput SwissProt OutputFile")
        sys.exit(2)
parser=MyParser()
parser.add_argument('foo', nargs='+')
args=parser.parse_args()

prioritized = open(sys.argv[1])
swissprot = open(sys.argv[2])

#some ugly code here for counting, on refactor to-do list. Finds # of PPANINI prioritized examples that are also in Uniref
pkeys = {}
j = 0
q = 0
for lines in prioritized:
	q += 1
	prekey = lines.rstrip().split('\t')[0]
	if prekey[0] == 'U':
		j += 1
		dval = prekey.split('_')[1]
		pkeys[dval] = 1
print str(q)+" in prioritized output."+'\n'
print str(j)+" also in Uniref."+'\n'

#Finds GO Terms in swissprot manually annotated data. SUPER ANNOYING as they don't have an API but a interactive downloader where columns are selected and downloaded to a \t deliminated file---in which they include \t's in some of their notes. REGEX use is to deal with this. 
prots = []
for lines in swissprot:
	m = re.findall(r"\GO:[0-9]*",lines.rstrip())
	prots.append((lines.rstrip().split('\t')[0],m))

#Using a dict to count frequence of GO terms occuring 
i = 0
godict = {}
for elems in prots:
	if pkeys.has_key(elems[0]):
		i += 1
		for terms in elems[1]:
			if godict.has_key(terms):
				godict[terms] += 1
			else:
				godict[terms] = 1


prioritized.close()
swissprot.close()
print str(i)+" also in SwissProt."+'\n'

#output the results to file
outputf =open(sys.argv[3],'w')
for elems in godict.keys():
	outputf.write(elems+'\t'+str(godict[elems])+'\n')

outputf.close()

#sort the outputted file 
os.system("sort -k2 -n "+sys.argv[3]+" >> "+sys.argv[3]+".sorted")
os.system("rm "+sys.argv[3])

print str(len(godict.keys()))+" GO Terms ranked for possible training"




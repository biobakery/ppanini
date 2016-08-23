import os
import sys
import re

'''Uses centroids from first file and 
produces mapped centroids to genes through mapper'''
foo = open(sys.argv[1])
centroids = [re.sub('[\t\r\n]','',i).strip() for i in foo]

mapper = open(sys.argv[2])
mapping = {}
for line in mapper:
	split_line = [re.sub('[\r\n\t]','', i) for i in line.split('\t')]
	try:
		mapping[split_line[0]] += [split_line[1]]
	except:
		mapping[split_line[0]] = [split_line[1]]

for i in centroids:
	if i in mapping:
		for j in mapping[i]:
			print '\t'.join([i,j])
	else:
		print i+'\tNA'
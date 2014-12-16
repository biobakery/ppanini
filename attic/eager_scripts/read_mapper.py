import os
import sys

indexer_files = os.listdir(sys.argv[1])
read_files = os.listdir(sys.argv[2])

#indexer_files = os.listdir('/n/huttenhower_lab_nobackup/downloads/HMP/HMASM/assemblies/PGAs/anterior_nares')

head = ['#! /bin/sh\n', '#BSUB -u shafquat@hsph.harvard.edu\n','#BSUB -q normal_serial\n','#BSUB -g /bowtie2\n']

#read_files = os.listdir('/n/huttenhower_lab_nobackup/downloads/HMP/HMASM/reads')
try:
	os.mkdir('asmbls_bsub')
except:
	pass
try: 
	os.mkdir('asmbls_err')
except:
	pass

try:
	os.mkdir('asmbls_out')
except:
	pass
for i in indexer_files:
	name = i.split('.')[0]
	tar_gz = name+'.tar.gz' in read_files
	foo = open('asmbls_bsub/'+name+'_bsub.bsub','w')
	foo.writelines(head)
	head_ii = ['#BSUB -J asmbls_'+name+'\n', '#BSUB -o asmbls_out/bwt2_'+name+'.out\n', '#BSUB -e asmbls_err/bwt2_'+name+'.err\n']
	foo.writelines(head_ii)
	foo.writelines(['\n','module load bio/bowtie2-2.0.0-beta6\n'])
	foo.writelines(['\n','module load bio/samtools-0.1.18\n'])
	foo.writelines(['\nbowtie2-build --quiet '+sys.argv[1]+'/'+i+' '+name+'_index\n'])
	if tar_gz:
		foo.writelines(['\ntar -xOvf '+sys.argv[2]+'/'+name+'.tar.gz | bowtie2 -x '+name+'_index -U - --no-unal --sensitive | samtools view -bS - > '+name+'.bam\n'])
	else:
		foo.writelines(['\ntar -xOvf '+sys.argv[2]+'/'+name+'.tar.bz2 | bowtie2 -x '+name+'_index -U - --no-unal --sensitive | samtools view -bS - > '+name+'.bam\n'])
	#foo.writelines(['\nrm '+name+'*.bt2\n'])
	foo.close()


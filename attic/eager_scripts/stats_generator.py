import os
import sys
indexer_files = os.listdir(sys.argv[1])
#indexer_files = os.listdir('/n/huttenhower_lab_nobackup/downloads/HMP/HMGI/stool/scripts/output')

head = ['#! /bin/sh\n', '#BSUB -u shafquat@hsph.harvard.edu\n','#BSUB -q normal_serial\n','#BSUB -g /samtools\n']

try:
	os.mkdir('sorter_bsubs')
except:
	pass

try:
	os.mkdir('sorter_err')
except:
	pass
try:
	os.mkdir('sorter_out')
except:
	pass
try:
        os.mkdir('sorted_bams')
except:
        pass
try:
        os.mkdir('bwt2_stats')
except:
        pass

for i in indexer_files:
	name = i.split('.')[0]
	foo = open('sorter_bsubs/'+name+'_bsub.bsub','w')
	foo.writelines(head)
	head_ii = ['#BSUB -J st2_'+name+'\n', '#BSUB -o sorter_out/st2_'+name+'.out\n', '#BSUB -e sorter_err/st2_'+name+'.err\n']
	foo.writelines(head_ii)
	foo.writelines(['\n','module load bio/samtools-0.1.18\n'])
	foo.writelines(['\nsamtools sort '+sys.argv[1]+'/'+i+' sorted_bams/'+i+'.sorted\n'])
	foo.writelines(['\nsamtools index sorted_bams/'+i+'.sorted.bam\n'])
	foo.writelines(['\nsamtools idxstats sorted_bams/'+i+'.sorted.bam > bwt2_stats/stats_'+i+'.txt\n'])
	#foo.writelines(['\nif [-f '+i+'.sorted.bam ];\nthen\nrm output/'+i+'\nfi\n'])
	foo.close()

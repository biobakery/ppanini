ASSEMBLIES=/n/huttenhower_lab_nobackup/data/hmp/anares_annotation/fna_only
READS=/n/huttenhower_lab_nobackup/downloads/HMP/HMASM/reads

#python read_mapper.py $ASSEMBLIES $READS

#ASMBLS_BSUB=asmbls_bsub

#files=$(ls $ASMBLS_BSUB)
#read -a array_files <<< $files
#for element in "${array_files[@]}"
#do bsub < $ASMBLS_BSUB/$element
#done

python stats_generator.py bams

files=$(ls sorter_bsubs)
read -a array_files <<< $files
for element in "${array_files[@]}"
do bsub < sorter_bsubs/$element
done



BODY_SITE="anteriornares"
DATA_FOLDER=/n/huttenhower_lab_nobackup/data/hmp/anares_annotation/data_files
RESULTS_U90=/n/huttenhower_lab_nobackup/data/hmp/anares_annotation/results/results_u90
SCRIPTS_FOLDER=/n/huttenhower_lab_nobackup/data/hmp/anares_annotation/src
ext_fasta=$DATA_FOLDER/extracted_$BODY_SITE/extracted_$BODY_SITE.txt

U50_DB=/n/huttenhower_lab_nobackup/data/hmp/stool_annotation/data_files/db_files/uniref50_udb

#cat $RESULTS_U90/*_u90.* | cut -f1 | sort| uniq > $DATA_FOLDER/pab_centroids_$BODY_SITE.extracted.txt
#cat $ext_fasta | grep '>' -n > $DATA_FOLDER/inds_headers_$BODY_SITE.extracted.txt
#cat $DATA_FOLDER/inds_headers_$BODY_SITE.extracted.txt | cut -f1 -d: > $DATA_FOLDER/inds_all_extracted_$BODY_SITE.txt 
#python $SCRIPTS_FOLDER/select_ids.py $DATA_FOLDER/pab_centroids_$BODY_SITE.extracted.txt $DATA_FOLDER/inds_headers_$BODY_SITE.extracted.txt $DATA_FOLDER/u90_$BODY_SITE.id.txt

#python $SCRIPTS_FOLDER/id_u50.py $DATA_FOLDER/u90_$BODY_SITE.id.txt $DATA_FOLDER/inds_all_extracted_$BODY_SITE.txt $DATA_FOLDER/inds_u50.txt
#mkdir $DATA_FOLDER/extracted_u50
#python $SCRIPTS_FOLDER/extract_ind.py $ext_fasta ind $DATA_FOLDER/inds_u50.txt $DATA_FOLDER/inds_all_extracted_$BODY_SITE.txt $DATA_FOLDER/extracted_u50 $BODY_SITE
#python $SCRIPTS_FOLDER/commander.py $DATA_FOLDER/extracted_u50 results/results_u50 $U50_DB u50


if [ ! -d err ]; then mkdir err; mkdir out; fi
date
files=$(ls bsub_u50)
read -a array_files <<< $files
for element in "${array_files[@]}"
do bsub < bsub_u50/$element
done

date
echo Loading files

####################################################
# Intializing variables
####################################################

BODY_SITE="anteriornares"
DATA_FOLDER=/n/huttenhower_lab_nobackup/data/hmp/anares_annotation/data_files
BAM_FOLDER=/n/huttenhower_lab_nobackup/data/hmp/anares_annotation/bwt2_stats
GFF3_FOLDER=/n/huttenhower_lab_nobackup/downloads/HMP/HMGI/anterior_nares
CLUSTER_FILE=/n/huttenhower_lab/data/hmp/index_to_cluster_mapping/anterior_nares.uc
STOOL_FASTA=/n/huttenhower_lab_nobackup/downloads/HMP/HMGC/anterior_nares_nr.fasta
SCRIPTS_FOLDER=/n/huttenhower_lab_nobackup/data/hmp/anares_annotation/src
U90_DB=/n/huttenhower_lab_nobackup/data/hmp/stool_annotation/data_files/db_files/uniref90_udb
U50_DB=/n/huttenhower_lab_nobackup/data/hmp/stool_annotation/data_files/db_files/uniref50_udb


echo "Creating gene only gff3 files to $DATA_FOLDER/gff3_gene"

#files=$(ls $GFF3_FOLDER)
#read -a array_files <<< $files

#if [ ! -d $DATA_FOLDER/gff3_gene ];
#then mkdir $DATA_FOLDER/gff3_gene
#fi

#for element in "${array_files[@]}"
#do cat $GFF3_FOLDER/$element | grep -w gene > $DATA_FOLDER/gff3_gene/$element
#done

echo "Cropping $BODY_SITE cluster file"
#cat $CLUSTER_FILE | grep -w H > $DATA_FOLDER/cropped_$BODY_SITE.uc

################################################
# Extracting prevalent and abundant centroids
################################################

echo "Processing data"
#python $SCRIPTS_FOLDER/abund_hmgc.py $DATA_FOLDER/gff3_gene $DATA_FOLDER/cropped_$BODY_SITE.uc $BAM_FOLDER $DATA_FOLDER


##########################################################
# Running annotations on prevalent and abundant centroids#
##########################################################

cat $DATA_FOLDER/prevalent_abundant_centroids_$BODY_SITE.txt | cut -f1 > $DATA_FOLDER/pab_centroids_$BODY_SITE.txt
cat $STOOL_FASTA | grep '>' -n > $DATA_FOLDER/inds_headers_$BODY_SITE.txt
cat $DATA_FOLDER/inds_headers_$BODY_SITE.txt | cut -f1 -d: > $DATA_FOLDER/inds_$BODY_SITE.txt
python $SCRIPTS_FOLDER/select_ids.py $DATA_FOLDER/pab_centroids_$BODY_SITE.txt $DATA_FOLDER/inds_headers_$BODY_SITE.txt $DATA_FOLDER/prevabundant_inds_$BODY_SITE.txt

echo "Indices found"

python $SCRIPTS_FOLDER/extract_ind.py $STOOL_FASTA ind $DATA_FOLDER/prevabundant_inds_$BODY_SITE.txt $DATA_FOLDER/inds_$BODY_SITE.txt $DATA_FOLDER $BODY_SITE

no=$(cat data_files/prevabundant_inds_$BODY_SITE.txt | wc -l)
if [ ! -d $DATA_FOLDER/extracted_$BODY_SITE ]; then mkdir $DATA_FOLDER/extracted_$BODY_SITE; fi
mv $DATA_FOLDER/extracted_$BODY_SITE.txt $DATA_FOLDER/extracted_$BODY_SITE/


if [ ! -d results ]; then mkdir results; mkdir results/results_u50; mkdir results/results_u90; fi

if [ $no -gt 10000 ]; then 
if [ ! -d $DATA_FOLDER/extracted_$BODY_SITE/divided_$BODY_SITE ]; then mkdir $DATA_FOLDER/extracted_$BODY_SITE/divided_$BODY_SITE; fi
cat $DATA_FOLDER/extracted_$BODY_SITE/extracted_$BODY_SITE.txt | grep '>' -n | cut -f1 -d:> $DATA_FOLDER/extracted_$BODY_SITE.inds.txt
python $SCRIPTS_FOLDER/extract_ind.py $DATA_FOLDER/extracted_$BODY_SITE/extracted_$BODY_SITE.txt n 10000 $DATA_FOLDER/extracted_$BODY_SITE.inds.txt $DATA_FOLDER/extracted_$BODY_SITE/divided_$BODY_SITE $BODY_SITE
python $SCRIPTS_FOLDER/commander.py $DATA_FOLDER/extracted_$BODY_SITE/divided_$BODY_SITE results/results_u90 $U90_DB u90
else
python $SCRIPTS_FOLDER/commander.py $DATA_FOLDER/extracted_$BODY_SITE results/results_u90 $U90_DB u90
fi

#files=$(ls bsub_u90)
#read -a array_files <<< $files
#for element in "${array_files[@]}"
#do bsub < bsub_u90/$element
#done

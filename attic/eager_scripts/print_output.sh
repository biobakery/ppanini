
BODY_SITE="anteriornares"
DATA_FOLDER=/n/huttenhower_lab_nobackup/data/hmp/anares_annotation/data_files
RESULTS_U50=/n/huttenhower_lab_nobackup/data/hmp/anares_annotation/results/results_u50
RESULTS_U90=/n/huttenhower_lab_nobackup/data/hmp/anares_annotation/results/results_u90
SCRIPTS_FOLDER=/n/huttenhower_lab_nobackup/data/hmp/anares_annotation/src
HMGI_FOLDER=/n/huttenhower_lab_nobackup/downloads/HMP/HMGI/anterior_nares

files=$(ls $HMGI_FOLDER)
read -a array_files <<< $files
mkdir $DATA_FOLDER/gff3_text
for ele in "${array_files[@]}"
do
cat $HMGI_FOLDER/$ele | grep product > $DATA_FOLDER/gff3_text/$ele
done

echo 'gff3 text completed'

cat $RESULTS_U50/*_u50.* $RESULTS_U90/*_u90.* | cut -f3 -d'|' > $DATA_FOLDER/annotation_results_$BODY_SITE.txt

python $SCRIPTS_FOLDER/print_table.py $DATA_FOLDER/prevalent_abundant_centroids_$BODY_SITE.txt $DATA_FOLDER/annotation_results_$BODY_SITE.txt $DATA_FOLDER/prevabundannotate_centroids_$BODY_SITE.tmp.txt $DATA_FOLDER/gff3_text $DATA_FOLDER

#nos=$(cat data_files/prevabundannotate_centroids_$BODY_SITE.txt | cut -f5 |grep -v . -n| cut -f1 -d:)
#read -a array_nos <<< $nos
#for element in "${array_nos[@]}";do x=$(echo $element'!d'); sed $x data_files/prevabundannotate_centroids_$BODY_SITE.txt; done > data_files/unannotated.centroids.txt


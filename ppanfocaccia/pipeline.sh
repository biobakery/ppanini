#!/bin/sh
echo "starting pipeline..."
nohup python selectGOTerms.py $1 $2 aterms
GOTERM=$(tail -1 aterms.sorted | cut -f 1)
rm -rf aterms.sorted  
nohup python goTermToTrainingAndTestSets.py $GOTERM $1 $2 $3 $4
nohup python train.py
./classify.sh 
nohup python mlEvaluationPlots.py 
./outputFolderOrganization.sh 
echo "pipeline succesfull!" 
   

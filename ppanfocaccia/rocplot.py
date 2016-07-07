print(__doc__)

import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.cross_validation import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from scipy import interp
import sys
global dist,label

dist = []
label = []

def genVectors(fname1, fname2):
	global dist,label
	labels = open(fname1)
	dists = open(fname2)
	label2 = []
	dist2 = []
	for lines in labels:
		label2.append(int(lines.rstrip().split(" ")[0]))
	for lines in dists:
		dist2.append(float(lines.rstrip()))
	dist.append(dist2)
	label.append(label2)

text = []

genVectors("samplePipelineOutput/svmtestfiles/fold1HeldOut","samplePipelineOutput/svmpredictions/fold1HeldOutOutput")
text.append("Fold 1 Held Out") 
genVectors("samplePipelineOutput/svmtrainingfiles/fold1","samplePipelineOutput/svmpredictions/fold1TrainingOutput")
text.append("Fold 1 Training")
genVectors("samplePipelineOutput/svmtestfiles/fold2HeldOut","samplePipelineOutput/svmpredictions/fold2HeldOutOutput")
text.append("Fold 2 Held Out")
genVectors("samplePipelineOutput/svmtrainingfiles/fold2","samplePipelineOutput/svmpredictions/fold2TrainingOutput")
text.append("Fold 2 Training")
genVectors("samplePipelineOutput/svmtestfiles/fold3HeldOut","samplePipelineOutput/svmpredictions/fold3HeldOutOutput")
text.append("Fold 3 Held Out")
genVectors("samplePipelineOutput/svmtrainingfiles/fold3","samplePipelineOutput/svmpredictions/fold3TrainingOutput")
text.append("Fold 3 Training")
genVectors("samplePipelineOutput/svmtestfiles/fold4HeldOut","samplePipelineOutput/svmpredictions/fold4HeldOutOutput")
text.append("Fold 4 Held Out")
genVectors("samplePipelineOutput/svmtrainingfiles/fold4","samplePipelineOutput/svmpredictions/fold4TrainingOutput")
text.append("Fold 4 Training")
genVectors("samplePipelineOutput/svmtrainingfiles/full","samplePipelineOutput/svmpredictions/fullOutput")
text.append("Full")

fpr = dict()
tpr = dict()
roc_auc = dict()
for i in range(9):
    fpr[i], tpr[i], _ = roc_curve(label[i], dist[i])
    roc_auc[i] = auc(fpr[i], tpr[i])

for i in range(9):
	plt.plot(fpr[i], tpr[i], label=text[i]+' (area = %0.2f)' % roc_auc[i])
plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic Curves')
plt.legend(loc="lower right",prop={'size':10})
plt.show()
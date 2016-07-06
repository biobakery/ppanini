#!/bin/sh
#runs joachims svm_light classifer (http://www.cs.cornell.edu/people/tj/) on all folds
svm_classify full fullModel fullOutput
svm_classify fold1HeldOut fold1Model fold1HeldOutOutput
svm_classify fold2HeldOut fold2Model fold2HeldOutOutput
svm_classify fold3HeldOut fold3Model fold3HeldOutOutput
svm_classify fold4HeldOut fold4Model fold4HeldOutOutput
svm_classify fold1 fold1Model fold1TrainingOutput
svm_classify fold2 fold2Model fold2TrainingOutput
svm_classify fold3 fold3Model fold3TrainingOutput
svm_classify fold4 fold4Model fold4TrainingOutput

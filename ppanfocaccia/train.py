import multiprocessing 
import os

#function to run Joachims(http://www.cs.cornell.edu/people/tj/) svmLight classifier by jumping out to a shell
def learn((inputf,outputf)):
	os.system("nohup svm_learn "+inputf+" "+outputf+" &")
	return True 	

#doing this in parallel (because finding best margin sep is O(n^3)) by setting up a multiprocessing pool of 5 processors and having 4 train folds and 1 train full model
p = multiprocessing.Pool(5)
p.map(learn,[("fold1","fold1Model"),("fold2","fold2Model"),("fold3","fold3Model"),("fold4","fold4Model"),("full","fullModel")])


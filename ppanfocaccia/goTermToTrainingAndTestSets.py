import sys
import random
import argparse

#help if user selected wrong command line arguments
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write("python goTermToTrainingAndTestSets.py GOTerm PPANINIOutputFile SwissProtFile POS/NEGfeatureRatio RandomnessSeed")
        sys.exit(2)
parser=MyParser()
parser.add_argument('foo', nargs='+')
args=parser.parse_args()

#handling command line arguments 
term = sys.argv[1]
ppaniniOutput = open(sys.argv[2],'r')
swissProt = open(sys.argv[3],'r')
pos2NegFeatureRatio = int(sys.argv[4])
randomnessSeed = int(sys.argv[5])

#instantiating data structures 
random.seed(randomnessSeed)
posExamples = set()
negExamples = set()
clusters = set()
possibleNegExamples = []

#adding known clusters from uncharacterized ones from PPANINI Output. known ones -> set; uncharacterized -> array 
for lines in ppaniniOutput:
        prekey = lines.rstrip().split('\t')[0]
        if prekey[0] == 'U':
                clusters.add(prekey.split('_')[1])
	else:
		possibleNegExamples.append(prekey)
		
#if priortized characterized cluster also in swissprot and is annotated to term binary classifier is working on add to pos examples
for lines in swissProt:
	if term in lines:
		candidateCluster = lines.rstrip().split('\t')[0]
		if candidateCluster in clusters:
			posExamples.add(candidateCluster)

#randomly select n negative examples from uncharacterized set of prioritized clusters & add to set
for i in range(0,pos2NegFeatureRatio*len(posExamples)):
	negExamples.add(random.choice(possibleNegExamples))

ppaniniOutput = open(sys.argv[2],'r')

negativeFeatures = []
positiveFeatures = []

#format positive and negative training examples and store in array 
for lines in ppaniniOutput:
	processedLine = lines.rstrip().split('\t')
        prekey = processedLine[0]
        if prekey[0] == 'U':
		key = prekey.split('_')[1]
                if key in posExamples:
			positiveFeatures.append("1 1:"+processedLine[1]+" 2:"+processedLine[2]+" 3:"+processedLine[3]+" 4:"+processedLine[4]+" #"+key)
        else:
                if prekey in negExamples:
			negativeFeatures.append("-1 1:"+processedLine[1]+" 2:"+processedLine[2]+" 3:"+processedLine[3]+" 4:"+processedLine[4]+" #"+prekey)

fold1 = []
fold2 = []
fold3 = []
fold4 = []

f = open("full",'w')

#function to assign 1 labeled Example to 1 of 4 folds randomly 
def assignFold(labeledExample):
	fold = random.randint(1,4)
	if fold == 1:
		fold1.append(labeledExample)
	if fold == 2:
		fold2.append(labeledExample)
	if fold == 3:
		fold3.append(labeledExample)
	if fold == 4:	
		fold4.append(labeledExample)

#assign all positive features to folds 
for elems in positiveFeatures:
	assignFold(elems)
	f.write(elems+'\n')

#assign all negative features to folds 
for elems in negativeFeatures:
	assignFold(elems) 
	f.write(elems+'\n')


#functions to write training & held out files 
def writeTrainingFold(fa,fb,fc,fileName):
	fileHandle = open(fileName,"w")
	for elems in fa:
		fileHandle.write(elems+'\n')	
	for elems in fb:
                fileHandle.write(elems+'\n')
	for elems in fc:
                fileHandle.write(elems+'\n')
	fileHandle.close()
	
def writeHeldOutFold(fa,fileName):
	fileHandle = open(fileName,"w")
	for elems in fa:
                fileHandle.write(elems+'\n')
	fileHandle.close()

#actually writing training & held out files
writeTrainingFold(fold2,fold3,fold4,"fold1")
writeTrainingFold(fold1,fold3,fold4,"fold2")
writeTrainingFold(fold2,fold1,fold4,"fold3")
writeTrainingFold(fold2,fold3,fold1,"fold4")
writeHeldOutFold(fold1,"fold1HeldOut")
writeHeldOutFold(fold2,"fold2HeldOut")
writeHeldOutFold(fold3,"fold3HeldOut")
writeHeldOutFold(fold4,"fold4HeldOut")
		
ppaniniOutput.close()
f.close()




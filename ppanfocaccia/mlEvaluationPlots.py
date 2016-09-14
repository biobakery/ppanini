#this program will output PR curves, convexed hulled PR curves, ROC curves, and data files to help in evaluating machine learning model performance. Needs some refactoring. 
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import time
import sys
import numpy as np
import matplotlib.cm as mplcm
import matplotlib.colors as colors
NUM_COLORS = 9 

cm = plt.get_cmap('gist_rainbow')
cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)

fontP = FontProperties()
fontP.set_size('small')
def generator(cutoff):
        global count,origdata,recalla,precisiona,fpra,tpra,fp,fn,tp,tn
	fp = 0
	fn = 0
	tp = 0 
	tn = 0 
	for y in range (0,count):
        	if origarray[y] > cutoff and origdata[y] == '1':
                	tp = tp + 1
                if origarray[y] > cutoff and origdata[y] == '-':
                        fp = fp + 1
                if origarray[y] <= cutoff and origdata[y] == '1':
                        fn = fn + 1
        tn = count - fn - fp - tp
	recall = (float(tp)/(float(tp)+float(fn)))
	precision = float(tp)/(float(tp)+float(fp))
	fpr = float(fp)/(float(tn)+float(fp))
	tpr = float(tp)/(float(tp)+float(tn))
	fpra.append(fpr)
	recalla.append(recall)
	precisiona.append(precision)
	tpra.append(tpr)
        lst = [str(cutoff),str(tp),str(tn),str(fp),str(fn),str(recall),str(precision),str(fpr),str(tpr)]
        out = " ".join(lst)
        gd.write(out)
	gd.write('\n')

def generator2(cutofftup):
	global fp,tp,tn,fn,fpra,recalla,precisiona,tpra
	if cutofftup[1] == '-':
		#if fp == 0:
		#	fp = fp
		#else:
		tn = tn + 1
		fp = fp - 1
	else:
		tp = tp - 1
		fn = fn + 1
 
        recall = (float(tp)/(float(tp)+float(fn)))
        precision = np.divide(float(tp),float(tp)+float(fp))
	#if str(precision) == 'nan':
	#	precision = float(np.divide(int(float(tp)),int(float(tp)+float(fp))))		
	fpr = float(fp)/(float(tn)+float(fp))
        tpr = recall
	fpra.append(fpr)
        recalla.append(recall)
        precisiona.append(precision)
        tpra.append(tpr)
        lst = [str(cutofftup[0]),str(tp),str(tn),str(fp),str(fn),str(recall),str(precision),str(fpr),str(tpr)]
        out = " ".join(lst)
        gd.write(out)
        gd.write('\n')


def convexHulling():						
	for x in range(1, count):
		if chprecisiona[x] < chprecisiona[x-1]:
			chprecisiona[x] = chprecisiona[x-1]
		

def auprc(xval,yval):
	ar = 0
	xs = 0
	for x in range (1,count):
		xsv = xval[x] - xs
		ar1 = xsv * yval[x]
		ar = ar + ar1
		xs = xs+xsv
	return ar

def clear():
	global chprecisiona,chrecalla,recalla,precisiona,fpra,tpra
	chprecisiona = []
	chrecalla = []
	precisiona = []
	fpra = []
	tpra = []
	
def integration():
	global auprcv,aucv,auprcs,aucs
	auprcv = auprc(chprecisiona,chrecalla)
        aucv = auprc(fpra[::-1],tpra[::-1])
        auprcs = "auprc: " + str(auprcv)
        aucs = "auc: " + str(aucv)	

			
def graph():
	global auprcv,aucv,auprcs,aucs	
	plt.subplots_adjust(hspace=1.2)
	plt.figure(1)
	ax1 = plt.subplot(312)
	ax1.set_color_cycle([scalarMap.to_rgba(i) for i in range(NUM_COLORS)])
	i = 0
	plt.text(0.5, 1.5, 'Precision-Recall Curves From 5-Fold Validation', horizontalalignment='center',fontsize=12,transform = ax1.transAxes)
	work("fullOutput","full","fullData")
	ffpra = fpra
        ftpra = tpra
        fchprecisiona = chprecisiona
        fchrecalla = chrecalla
        faucs = aucs
        fauprcs = auprcs
	plt.plot(recalla,precisiona)
	work("fold1HeldOutOutput","fold1HeldOut","fold1HeldOutData")
	fa1fpra = fpra
	fa1tpra = tpra
	fa1chprecisiona = chprecisiona
	fa1chrecalla = chrecalla
	fa1aucs = aucs
	fa1auprcs = auprcs
        plt.plot(recalla,precisiona)
	work("fold2HeldOutOutput","fold2HeldOut","fold2HeldOutData")
	fa2fpra = fpra
        fa2tpra = tpra
        fa2chprecisiona = chprecisiona
        fa2chrecalla = chrecalla
        fa2aucs = aucs
        fa2auprcs = auprcs
        plt.plot(recalla,precisiona)
	work("fold3HeldOutOutput","fold3HeldOut","fold3HeldOutData")
	fa3fpra = fpra
        fa3tpra = tpra
        fa3chprecisiona = chprecisiona
        fa3chrecalla = chrecalla
        fa3aucs = aucs
        fa3auprcs = auprcs
        plt.plot(recalla,precisiona)
	work("fold4HeldOutOutput","fold4HeldOut","fold4HeldOutData")
	fa4fpra = fpra
        fa4tpra = tpra
        fa4chprecisiona = chprecisiona
        fa4chrecalla = chrecalla
        fa4aucs = aucs
        fa4auprcs = auprcs
        plt.plot(recalla,precisiona)
	work("fold1TrainingOutput","fold1","fold1TrainingData")
	fa1tfpra = fpra
        fa1ttpra = tpra
        fa1tchprecisiona = chprecisiona
        fa1tchrecalla = chrecalla
        fa1taucs = aucs
        fa1tauprcs = auprcs
	plt.plot(recalla,precisiona)
	work("fold2TrainingOutput","fold2","fold2TrainingData")
	fa2tfpra = fpra
        fa2ttpra = tpra
        fa2tchprecisiona = chprecisiona
        fa2tchrecalla = chrecalla
        fa2taucs = aucs
        fa2tauprcs = auprcs
        plt.plot(recalla,precisiona)
	work("fold3TrainingOutput","fold3","fold3TrainingData")
	fa3tfpra = fpra
        fa3ttpra = tpra
        fa3tchprecisiona = chprecisiona
        fa3tchrecalla = chrecalla
        fa3taucs = aucs
        fa3tauprcs = auprcs
        plt.plot(recalla,precisiona)
	work("fold4TrainingOutput","fold4","fold4TrainingData")
	fa4tfpra = fpra
        fa4ttpra = tpra
        fa4tchprecisiona = chprecisiona
        fa4tchrecalla = chrecalla
        fa4taucs = aucs
        fa4tauprcs = auprcs
        plt.plot(recalla,precisiona)
	ax1.legend( ('full', 'fold 1 test', 'fold 2 test', 'fold 3 test', 'fold 4 test', 'fold 1', 'fold 2', 'fold 3', 'fold 4'),fancybox=True,shadow=True, loc='upper center',prop = fontP,bbox_to_anchor=(0.5, 1.43),ncol=3)
	plt.axis([0,1,0,1.4])
	plt.xlabel('recall')
	plt.ylabel('precision')
	ax2 = plt.subplot(311)
	ax2.set_color_cycle([scalarMap.to_rgba(i) for i in range(NUM_COLORS)])
	plt.text(0.5, 1.50, 'ROC Curves From 5-Fold Validation', horizontalalignment='center',fontsize=12,transform = ax2.transAxes)
	#work("fullOutput","SVMData.txt","fulldata")
	#a1 = aucs
        plt.plot(ffpra,ftpra)
        #work("fa1output","fa1t","fa1data")
	#a2 = aucs
	plt.plot(fa1fpra,fa1tpra)
        #work("fa2output","fa2t","fa2data")
	#a3 = aucs
        plt.plot(fa2fpra,fa2tpra)
        #work("fa3output","fa3t","fa3data")
	#a4 = aucs
        plt.plot(fa3fpra,fa3tpra)
        #work("fa4output","fa4t","fa4data")
        plt.plot(fa4fpra,fa4tpra)
	plt.plot(fa1tfpra,fa1ttpra)
	plt.plot(fa2tfpra,fa2ttpra)
	plt.plot(fa3tfpra,fa3ttpra)
	plt.plot(fa4tfpra,fa4ttpra)
	plt.axis([0,1,0,1.8])
	plt.xlabel('fpr')
	plt.ylabel('tpr')
	ax2.legend( ('full;'+str(faucs)[0:10], 'fold 1 test;'+str(fa1aucs)[0:10], 'fold 2 test;'+str(fa2aucs)[0:10], 'fold 3 test;'+str(fa3aucs)[0:10], 'fold 4 test;'+str(fa4aucs)[0:10],'fold 1;'+str(fa1taucs)[0:10],'fold 2;'+str(fa2taucs)[0:10],'fold 3;'+str(fa3taucs)[0:10],'fold 4;'+str(fa4taucs)[0:10]),fancybox=True,shadow=True, prop = fontP,bbox_to_anchor=(.5, 1.43), loc='upper center', ncol=3)
	ax3 = plt.subplot(313)
	ax3.set_color_cycle([scalarMap.to_rgba(i) for i in range(NUM_COLORS)])
	plt.text(0.5, 1.50, 'Convex Hulls Of Precision-Recall Curves From 5-Fold Validation', horizontalalignment='center',fontsize=12,transform = ax3.transAxes)
	#work("fullOutput","SVMData.txt","fulldata")
	#b1 = auprcs
        plt.plot(fchrecalla,fchprecisiona)
        #work("fa1output","fa1t","fa1data")
	#b2 = auprcs
        plt.plot(fa1chrecalla,fa1chprecisiona)
        #work("fa2output","fa2t","fa2data")
	#b3 = auprcs
        plt.plot(fa2chrecalla,fa2chprecisiona)
        #b4 = auprcs
	#work("fa3output","fa3t","fa3data")
        plt.plot(fa3chrecalla,fa3chprecisiona)
        #work("fa4output","fa4t","fa4data")
        plt.plot(fa4chrecalla,fa4chprecisiona)
	plt.plot(fa1tchrecalla,fa1tchprecisiona)
	plt.plot(fa2tchrecalla,fa2tchprecisiona)
	plt.plot(fa3tchrecalla,fa3tchprecisiona)
	plt.plot(fa4tchrecalla,fa4tchprecisiona)
        plt.axis([0,1,0,1.4])
	plt.xlabel('recall')
	plt.ylabel('precision')
	#ax3.legend( ('fold 1 test;'+str(fa1auprcs)[0:12], 'fold 1 training;'+str(fa1tauprcs)[0:12]), fancybox=True,shadow=True, loc='upper center',prop = fontP,bbox_to_anchor=(.5, 1.43),ncol =3)
	ax3.legend( ('full;'+str(fauprcs)[0:12], 'fold 1 test;'+str(fa1auprcs)[0:12], 'fold 2 test;'+str(fa2auprcs)[0:12], 'fold 3 test;'+str(fa3auprcs)[0:12], 'fold 4 test;'+str(fa4auprcs)[0:12], 'fold 1;'+str(fa1tauprcs)[0:12], 'fold 2;'+str(fa2tauprcs)[0:12], 'fold 3;'+str(fa3tauprcs)[0:12], 'fold 4;'+str(fa4tauprcs)[0:12]),fancybox=True,shadow=True, loc='upper center',prop = fontP,bbox_to_anchor=(.5, 1.43),ncol =3)
	#plt.show()
	plt.savefig('graphs.png')

def work(outputfile,actualfile,gfile):
	ins = open (outputfile,"r")
        data = open(actualfile,"r")
        global origarray,count,origdata,recalla,precisiona,fpra,tpra,fp,fn,tn,tp,tupa,chprecisiona,chrecalla,aucs,auprc,aucv,auprv
        origarray = []
        array = []
        origdata = []
        recalla = []
        precisiona = []
        fpra = []
        tpra=[]
        tupa = []
        for line in ins:
                origarray.append(float(line[:-1]))
        ins.close()
        array = origarray[:]
        for line in data:
                origdata.append(line[0:1])
        count = len(origarray)
        for x in range (0, count):
                new = array[x],origdata[x]
                tupa.append(new)
        data.close
        tupa.sort()
	executeSingleFold(tupa,gfile)
	convexHulling()        
	integration()
	aucs = "auc: " + str(aucv)
	#
	#length = len(chprecisiona)
	#for x in range (0,length-1):
	#	if chprecisiona[x] > 1:
	#		chprecisiona[x] = 1
	#
def executeSingleFold(array,outputname):
	global count,gd,tpra,fpra,chprecisiona,chrecalla,precisiona,recalla
	gd = open(outputname,"w+")
	gd.write("Cutoff TP TN FP FN Recall Precision FPR TPR")
	gd.write('\n')
	count = len(array)
	generator(array[0][0]-1)
	for x in range (0,count):
		generator2(array[x])
	gd.close()
	chprecisiona = precisiona[:]
	chrecalla = recalla[:]
	tpra = [1] + tpra
	fpra = [1] + fpra

graph()



		
	


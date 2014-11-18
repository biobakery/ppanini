import os
import sys
import pdb
import re

def read_abund_folder(foldername):
	files = os.listdir(foldername)
	id_abund = {}
	for fname in files:
		id_abund[fname.split('.')[0].strip()] = {}
		foo = open(foldername+'/'+fname)
		foo = foo.readlines()
		for line in foo:
			id_abund[fname.split('.')[0].strip()][re.sub('[\r\t\n]','', line.split('\t')[0].strip())] = float(re.sub('[\n\t\r]','',line.split('\t')[1].strip()))
	return id_abund

def read_gi_folder(foldername):
	files = os.listdir(foldername)
	gi_reads = {}
	for fname in files:
		gi_reads[fname.split('.')[0].strip()] = {}
		foo = open(foldername+'/'+fname)
		foo = foo.readlines()
		for line in foo:
			gi_reads[fname.split('.')[0].strip()][re.sub('[\r\t\n]','', line.split('\t')[0].strip())] = re.sub('[\n\t\r]','',line.split('\t')[1].strip())
	return gi_reads

def read_prev_gc_file(filename):
	foo = open(filename)
	foo = foo.readlines()
	prev_gc = {}
	for i in foo:
		prev_gc[i.split('\t')[0].strip()] = float(re.sub('[\r\t\n]','',i.split('\t')[1].strip()))
	return prev_gc

def read_gc_file(filename):
	foo = open(filename)
	foo = foo.readlines()
	gc_gi = {}
	for i in foo:
		split_i = i.split('\t')
		gc_gi[split_i[0].strip()] = [k.strip() for k in split_i[1:]]
	return gc_gi


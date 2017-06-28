#! /usr/bin/env python

import sys, re, argparse
import subprocess

# argument parsing (python argparse)
parser = argparse.ArgumentParser()
parser.add_argument( 
    '-i', '--input', 
    help='fasta file', 
    )
parser.add_argument(
    '-f', '--ok_file',
    help='like grep -f behavior', 
    )
args = parser.parse_args()

ok = {}
with open( args.ok_file ) as fh:
    for line in fh:
        if line[0] != ">":
            line = ">"+line
        ok[line.replace('\r\n', '')] = 1
#print (ok)
# process the fasta file
total = 0
count = 0
skipping = True
with open( args.input ) as fh:
    for line in fh:
        if line[0] == ">":
            total += 1
            if line.split(" ")[0]in ok:
                skipping = False
                count += 1
            else:
                skipping = True
        if not skipping:
            print line.strip()

print >>sys.stderr, "ok seqs:", len( ok ), "total seqs:", total, "found:", count

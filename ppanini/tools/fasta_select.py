#! /usr/bin/env python

import sys, re, argparse
import subprocess

def  get_args():
    parser = argparse.ArgumentParser(
        description="PPANINI selecting sequence of genes form fasta and mapping file",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument( 
        '-i', '--input', 
        help='fasta file',
        required = True, 
        )
    parser.add_argument(
        '-f', '--ok_file',
        help='like grep -f behavior', 
        required = True,
        )
    parser.add_argument( 
        '-o', '--output', 
        help='fasta file',
        required = True, 
        )
    args = parser.parse_args()
    return args

def main():
    args = get_args()
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
    f1=open(args.output, 'w')
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
                f1.write(line)
                #print line.strip()
    
    print >>sys.stderr, "Hits:", len( ok ), "total sequences:", total, "found hits sequences:", count

if __name__ == '__main__':
    main()

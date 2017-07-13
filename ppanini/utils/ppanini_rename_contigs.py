#! /usr/bin/env python

import sys, re, argparse, os
import subprocess

def get_args():
    # argument parsing (python argparse)
    parser = argparse.ArgumentParser()
    parser.add_argument( 
        '-i', '--input', 
        help='fasta file', 
        )
    parser.add_argument( 
        '-o', '--output', 
        help='fasta file', 
        )
    args = parser.parse_args()
    return args
def main():
    args = get_args()
    f1=open(args.output, 'w')
    base=os.path.basename( args.input)
    sample_name = os.path.splitext(base)[0]
    with open( args.input ) as fh:
        for line in fh:
            if line[0] == ">":
                line = line.replace(">", (">"+sample_name+"_"))
                line = line.replace(".", "_", 1)
            f1.write(line)
if __name__ == '__main__':
    main()
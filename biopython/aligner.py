#!/usr/bin/env python3.4

'''
Name: aligner.py
Date: 2014-03-30
Author: Dan Shea <shea.d@husky.neu.edu>
Description: aligner.py utilizes BioPython to provide a convenient wrapper to
             various alignment programs.  This provides a standardization to
             the naming convention used for the input/output files, allowing us
             to automate various aspects of the data processing pipeline.
'''

from Bio import SeqIO
from Bio import Align
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
import os.path
import sys
import argparse

def runClustalw():
    pass

def runMuscle():
    pass

def runMafft():
    pass

def main():
    # create an argparse parser
    parser = argparse.ArgumentParser()
    # add positional arguments to the parser
    # alignment program to use
    parser.add_argument('aligner', nargs='?', choices=['clustal', 'muscle', 'mafft'], default='clustal', help='aligner to use')
    # if infile is specified, open that file for reading, otherwise use stdin
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='input file')
    # if outfile is specified, open that file for writing, otherwise use stdout
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='output file')
    args = parser.parse_args()
    return(0)

if __name__ == '__main__':
    sys.exit(main())
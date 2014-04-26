#!/usr/bin/env python

################################################################################
#
# Name: rna.py
# Date: 2014-04-26
# Author: dshea
# Description: transcribe DNA to RNA
#
# Problem
#
# An RNA string is a string formed from the alphabet containing 'A', 'C', 'G',
# and 'U'.
#
# Given a DNA string t corresponding to a coding strand, its transcribed RNA
# string u is formed by replacing all occurrences of 'T' in t with 'U' in u.
#
# Given: A DNA string t having length at most 1000 nt.
#
# Return: The transcribed RNA string of t.
# Sample Dataset
#
# GATGGAACTTGACTACGTAAATT
#
# Sample Output
#
# GAUGGAACUUGACUACGUAAAUU
#
################################################################################

import sys
from Bio import Seq

def main():
    if len(sys.argv) != 2:
        print 'usage: {0:s} infile'.format(sys.argv[0])
        sys.exit(1)
    else:
        with open(sys.argv[1], 'r') as fh:
            dna = fh.readline().strip()
            print str(Seq.Seq(dna).transcribe())
    sys.exit(0)

if __name__ == '__main__':
    main()


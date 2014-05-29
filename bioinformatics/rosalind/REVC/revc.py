#!/usr/bin/env python

################################################################################
#
# Name: revc.py
# Date: 2014-04-26
# Author: dshea
# Description: reverse complement
#
# Problem
#
# In DNA strings, symbols 'A' and 'T' are complements of each other, as are 'C'
# and 'G'.
#
# The reverse complement of a DNA string s is the string sc formed by reversing
# the symbols of s, then taking the complement of each symbol (e.g., the reverse
# complement of "GTCA" is "TGAC").
#
# Given: A DNA string s of length at most 1000 bp.
#
# Return: The reverse complement sc of s.
# Sample Dataset
#
# AAAACCCGGT
#
# Sample Output
#
# ACCGGGTTTT
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
            print str(Seq.Seq(dna).reverse_complement())
    sys.exit(0)

if __name__ == '__main__':
    main()


#!/usr/bin/env python

################################################################################
#
# Name: rvco.py
# Date: 2014-04-26
# Author: dshea
# Description:
#
# Problem
#
# Recall that in a DNA string s, 'A' and 'T' are complements of each other, as
# are 'C' and 'G'. Furthermore, the reverse complement of s is the string sc
# formed by reversing the symbols of s and then taking the complement of each
# symbol (e.g., the reverse complement of "GTCA" is "TGAC").
#
# The Reverse Complement program from the SMS 2 package can be run online here.
#
# Given: A collection of n (n<=10) DNA strings.
#
# Return: The number of given strings that match their reverse complements.
# Sample Dataset
#
# >Rosalind_64
# ATAT
# >Rosalind_48
# GCATA
#
# Sample Output
#
# 1
#
################################################################################

import sys
from Bio import SeqIO

def rvco_matches(infile):
    total = 0
    for seq_record in SeqIO.parse(infile, format='fasta'):
        if str(seq_record.seq) == str(seq_record.reverse_complement().seq):
            total += 1
    return(total)

def main():
    if len(sys.argv) != 2:
        print 'usage: {0:s} infile'.format(sys.argv[0])
        sys.exit(1)
    else:
        infile = sys.argv[1]
        print rvco_matches(infile)
        sys.exit(0)

if __name__ == '__main__':
    main()


#!/usr/bin/env python

################################################################################
#
# Name: hamm.py
# Date: 2014-04-27
# Author: dshea
# Description:
#
# Problem
#
# Given two strings s and t of equal length, the Hamming distance between s and t,
# denoted dH(s,t), is the number of corresponding symbols that differ in s and t.
# 
# Given: Two DNA strings s and t of equal length (not exceeding 1 kbp).
#
# Return: The Hamming distance dH(s,t).
# Sample Dataset
#
# GAGCCTACTAACGGGAT
# CATCGTAATGACGGCCT
#
# Sample Output
#
# 7
#
################################################################################

import sys

def hamm(seq_a, seq_b):
    total = 0
    for a, b in zip(seq_a, seq_b):
        if a != b:
            total += 1
    return(total)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'usage: {0:s} infile'.format(sys.argv[0])
        sys.exit(1)
    infile = sys.argv[1]
    with open(infile, 'r') as fh:
        seq_a = fh.readline().strip()
        seq_b = fh.readline().strip()
        print hamm(seq_a, seq_b)
    sys.exit(0)
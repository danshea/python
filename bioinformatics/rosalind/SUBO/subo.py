#!/usr/bin/env python

################################################################################
# Name: subo.py
# Date: 2014-04-26
# Author: Dan Shea <shea.d@husky.neu.edu>
# Description:
# Problem
#
# The Lalign program for finding multiple alternative matches via suboptimal
# alignment is available here.
#
# Given: Two DNA strings s and t in FASTA format that share some short inexact
# repeat r of 32-40 bp. By "inexact" we mean that r may appear with slight
# modifications (each repeat differ by <=3 changes/indels).
#
# Return: The total number of occurrences of r as a substring of s, followed by
# the total number of occurrences of r as a substring of t.  Sample Dataset
#
# >Rosalind_12
# GACTCCTTTGTTTGCCTTAAATAGATACATATTTACTCTTGACTCTTTTGTTGGCCTTAAATAGATACATATTTGTGCGACTCCACGAGTGATTCGTA
# >Rosalind_37
# ATGGACTCCTTTGTTTGCCTTAAATAGATACATATTCAACAAGTGTGCACTTAGCCTTGCCGACTCCTTTGTTTGCCTTAAATAGATACATATTTG
#
# Sample Output
#
# 2 2
#
################################################################################

import sys
import subprocess
import shlex
import regex

def lalign(A, B):
    count = 0
    command = shlex.split('lalign36 -r 3 -m 10 {0:s} {1:s}'.format(A, B))
    output = subprocess.check_output(command)
    for line in output.split('\n'):
        if regex.search('lsw_overlap', line):
            if 32 <= int(line.split(':')[1].strip()) <= 40:
                count += 1
    return(count)

def main():
    if len(sys.argv) != 3:
        print 'usage: {0:s} input_a input_b'.format(sys.argv[0])
        sys.exit(1)
    else:
        a = sys.argv[1]
        b = sys.argv[2]
        print '{0:d} {1:d}'.format(lalign(a, b), lalign(b, a))
    sys.exit(0)
    
if __name__ == '__main__':
    main()
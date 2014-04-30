#!/usr/bin/env python

################################################################################
#
# Name: subs.py
# Date: 2014-04-30
# Author: dshea
# Description:
#
# Problem
#
# Given two strings s and t, t is a substring of s if t is contained as a
# contiguous collection of symbols in s (as a result, t must be no longer than s).
#
# The position of a symbol in a string is the total number of symbols found to its
# left, including itself (e.g., the positions of all occurrences of 'U' in
# "AUGCUUCAGAAAGGUCUUACG" are 2, 5, 6, 15, 17, and 18). The symbol at position i
# of s is denoted by s[i].
#
# A substring of s can be represented as s[j:k], where j and k represent the
# starting and ending positions of the substring in s; for example, if s =
# "AUGCUUCAGAAAGGUCUUACG", then s[2:5] = "UGCU".
#
# The location of a substring s[j:k] is its beginning position j; note that t will
# have multiple locations in s if it occurs more than once as a substring of s
# (see the Sample below).
#
# Given: Two DNA strings s and t (each of length at most 1 kbp).
#
# Return: All locations of t as a substring of s. Sample Dataset
#
# GATATATGCATATACTT
# ATAT
#
# Sample Output
#
# 2 4 10
#
################################################################################

import sys

def main():
    if len(sys.argv) != 2:
        print 'usage: {0:s} infile'.format(sys.argv[0])
        sys.exit(1)
    else:
        infile = sys.argv[1]
        with open(infile, 'r') as fh:
            results = list()
            s = fh.readline().strip()
            s_len = len(s)
            t = fh.readline().strip()
            t_len = len(t)
            for i in range(s_len - t_len):
                if s[i:i+t_len] == t:
                    results.append(i+1)
            for result in results:
                print result,
            print
        sys.exit(0)

if __name__ == '__main__':
    main()


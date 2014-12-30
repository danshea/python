#!/usr/bin/env python3.4
"""
Filename:     SSEQ.py
Date Created: 2014-12-31 01:22
Author:       dshea <danshea@iastate.edu>
Description:  Given:  Two DNA strings s and t (each of length at most 1 kbp) in
                      FASTA format.
              Return: One collection of indices of s in which the symbols of t
                      appear as a subsequence of s. If multiple solutions exist,
                      you may return any one.
"""
import argparse
import sys
from Bio import SeqIO

def readFasta(infile):
    """ Reads in record(s) from FASTA file returning a list of SeqRecord objects """
    records = list()
    for seq in SeqIO.parse(infile, "fasta"):
        records.append(seq.upper())
    return(records)

def main():
    """ The main() """
    parser = argparse.ArgumentParser(description='Description goes here.')
    parser.add_argument('-i', '--infile', nargs='?',
                        type=argparse.FileType('r'), default=sys.stdin)
    args = parser.parse_args()
    s,t = readFasta(args.infile.name)
    s = str(str(s.seq))
    t = str(str(t.seq))
    previous_offset = 0
    offset = 0
    solution = list()
    for nt in t:
        previous_offset = offset
        offset = s[offset:].find(nt)+1
        offset += previous_offset
        solution.append(offset)
    print('{}'.format(' '.join([str(i) for i in solution])))

if __name__ == '__main__':
    main()

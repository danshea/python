#!/usr/bin/env python3.4
"""
Filename:     REVP.py
Date Created: 2014-12-28 23:45
Author:       dshea <danshea@iastate.edu>
Description:  The file description goes here.
"""
import argparse
import sys
from Bio import SeqIO

def findReversePalindrome(kmers):
    """ Take a list of kmers and find all reverse palindromes. """
    # A kmer is a reverse palindrome if it equals its reverse complement
    position = 1
    for kmer in kmers:
        #print('{}\t{}\t{}'.format(position, str(kmer.seq),str(kmer.seq.reverse_complement())))
        if str(kmer.seq) == str(kmer.seq.reverse_complement()):
            print('{}\t{}'.format(position, len(kmer)))
        position = position + 1

def makeKmers(seq, k):
    """ Takes a Seq object and returns a list of all kmers of size k """
    return([seq[i:i+k] for i in range(0,len(seq)-k+1)])

def readFasta(infile):
    """ Reads in record from FASTA file returning a Seq object """
    record = SeqIO.read(infile, "fasta")
    return(record)

def main():
    """ The main() """
    parser = argparse.ArgumentParser(description='Description goes here.')
    parser.add_argument('-i', '--infile', nargs='?',
                        type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('-o', '--outfile', nargs='?',
                        type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()

    # Print out the argument values
    #print('infile {}'.format(args.infile.name))
    #print('outfile {}'.format(args.outfile.name))
    for k in range(4,13):
        findReversePalindrome(makeKmers(readFasta(args.infile.name), k))

if __name__ == '__main__':
    main()

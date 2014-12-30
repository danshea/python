#!/usr/bin/env python3.4
"""
Filename:     KMER.py
Date Created: 2014-12-30 23:43
Author:       dshea <danshea@iastate.edu>
Description:  Given:  A DNA string s in FASTA format (having length at most 100
                      kbp).
              Return: The 4-mer composition of s.


"""
import argparse
import sys
import numpy as np
from Bio import SeqIO

def makeAlphabet(symbols):
    """ Return a dictionary with an integer key for each of the (up to) 10
        possible symbols
    """
    alphabet = dict()
    i = 0
    for symbol in symbols:
        alphabet[i] = symbol
        i += 1
    return(alphabet)

def enumerateStrings(alphabet, n):
    """ Given:  A collection of at most 10 symbols defining an ordered
                alphabet, and a positive integer n (n≤10).
        Return: All strings of length n that can be formed from the
                alphabet, ordered lexicographically.
    """
    base = len(alphabet.keys())
    enumeration = [('{0:0'+str(n)+'d}').format(int(np.base_repr(i, base))) for i in range(base**n)]
    return(enumeration)

def translateEnumeration(enumerations,alphabet):
    """ Given:  an enumeration of possible strings and a dictionary providing
                the symbol lookup
        Return: The possible strings
    """
    words = list()
    for enumeration in enumerations:
        word = list()
        for symbol in list(enumeration):
            word.append(alphabet[int(symbol)])
        word = ''.join(word)
        words.append(word)
    return(words)

def readFasta(infile):
    """ Reads in record from FASTA file returning a SeqRecord object """
    seq = SeqIO.read(infile, "fasta")
    return(seq.upper())

def countMatches(word, seq):
    """ Given:  a word and a sequence
        Return: the number of times the word is found in the sequence
    """
    count = 0
    for kmer in [seq[i:i+4] for i in range(0,len(seq)-3)]:
        if word == kmer:
            count += 1
    return(count)

def main():
    """ The main() """
    parser = argparse.ArgumentParser(description='Description goes here.')
    parser.add_argument('-i', '--infile', nargs='?',
                        type=argparse.FileType('r'), default=sys.stdin)
    args = parser.parse_args()
    alphabet = makeAlphabet(['A','C','G','T'])
    words = translateEnumeration(enumerateStrings(alphabet, 4), alphabet)
    seq = readFasta(args.infile.name)
    kcomp = list()
    for word in words:
        kcomp.append(countMatches(word, str(seq.seq)))
    print(' '.join([str(i) for i in kcomp]))

if __name__ == '__main__':
    main()

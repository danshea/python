#!/usr/bin/env python3.4
"""
Filename:     LEXV.py
Date Created: 2015-01-12 00:25
Author:       dshea <danshea@iastate.edu>
Description:  Given:  A permutation of at most 12 symbols defining an ordered
                      alphabet 𝒜 and a positive integer n (n≤4).
              Return: All strings of length at most n formed from 𝒜, ordered
                      lexicographically. (Note: As in “Enumerating k-mers
                      Lexicographically”, alphabet order is based on the order
                      in which the symbols are given.)
"""
import argparse
import sys
import numpy as np

def makeAlphabet(symbols):
    """ Return a dictionary with an integer key for each of the (up to) 10
        possible symbols
    """
    alphabet = dict()
    i = 0
    for symbol in symbols:
        alphabet['{0:x}'.format(i)] = symbol
        i += 1
    return(alphabet)

def enumerateStrings(alphabet, n):
    """ Given:  A collection of at most 12 symbols defining an ordered
                alphabet, and a positive integer n (n≤4).
        Return: All strings of length at most n that can be formed from the
                alphabet, ordered lexicographically.
    """
    base = len(alphabet.keys())
    enumeration = list()
    for i in range(0,n+1):
        enumeration += [('{0:0'+str(i)+'x}').format(int(np.base_repr(j, base),16)) for j in range(base**n)]
    enumeration = sorted(list(set(enumeration)))
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
            word.append(alphabet[symbol])
        word = ''.join(word)
        words.append(word)
    return(words)

def main():
    """ The main() """
    parser = argparse.ArgumentParser(description='Description goes here.')
    parser.add_argument('-i', '--infile', nargs='?',
                        type=argparse.FileType('r'), default=sys.stdin)
    args = parser.parse_args()
    symbols = args.infile.readline()
    n = args.infile.readline()
    symbols = list(symbols.strip().replace(' ',''))
    n = int(n.strip())
    alphabet = makeAlphabet(symbols)
    words = translateEnumeration(enumerateStrings(alphabet, n), alphabet)
    for word in words:
        print('{}'.format(word))

if __name__ == '__main__':
    main()

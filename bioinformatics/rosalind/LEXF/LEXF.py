#!/usr/bin/env python3.4
"""
Filename:     LEXF.py
Date Created: 2014-12-30 22:08
Author:       dshea <danshea@iastate.edu>
Description:  Given: A collection of at most 10 symbols defining an ordered
              alphabet, and a positive integer n (n≤10).
              Return: All strings of length n that can be formed from the
              alphabet, ordered lexicographically.
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

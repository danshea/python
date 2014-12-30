#!/usr/bin/env python3.4
"""
Filename:     SIGN.py
Date Created: 2014-12-30 19:57
Author:       dshea <danshea@iastate.edu>
Description:  Given: A positive integer n≤6.
              Return: The total number of signed permutations of length n,
              followed by a list of all such permutations (you may list the
              signed permutations in any order).
"""
import argparse
import sys
import itertools

def main():
    """ The main() """
    parser = argparse.ArgumentParser(description='Description goes here.')
    parser.add_argument("n", type=int)
    args = parser.parse_args()
    n = args.n
    positive = list(range(1,n+1))
    negative = [-i for i in positive]
    perms = list(itertools.permutations(positive+negative, n))
    # Filter permutations where duplicates exist (i.e. -1 and 1 in same output)
    results = []
    for perm in perms:
        tmp = [abs(i) for i in perm]
        tmp = set(tmp)
        if len(tmp) == len(perm):
            results.append(perm)
    print('{}'.format(len(results)))
    for result in results:
        print('{}'.format(' '.join([str(i) for i in result])))

if __name__ == '__main__':
    main()

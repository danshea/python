#!/usr/bin/env python3.4
"""
Filename:     PERM.py
Date Created: 2014-12-30 19:29
Author:       dshea <danshea@iastate.edu>
Description:  Given: A positive integer n≤7.
              Return: The total number of permutations of length n, followed by
              a list of all such permutations (in any order).
"""
import argparse
import sys
import itertools

def main():
    """ The main() """
    parser = argparse.ArgumentParser(description='Description goes here.')
    parser.add_argument("n", type=int)
    args = parser.parse_args()
    perms = list(itertools.permutations(range(1,args.n+1)))
    print('{}'.format(len(perms)))
    for perm in perms:
        print('{}'.format(' '.join([str(s) for s in perm])))

if __name__ == '__main__':
    main()

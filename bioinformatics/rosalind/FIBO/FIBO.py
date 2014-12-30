#!/usr/bin/env python3.4
"""
Filename:     FIBO.py
Date Created: 2014-12-31 02:10
Author:       dshea <danshea@iastate.edu>
Description:  Return n-th number in the Fibonacci Sequence
"""
import argparse
import sys

def fib():
    a = 0
    b = 1
    while(True):
        yield(a)
        a,b = b, a+b

def main():
    """ The main() """
    parser = argparse.ArgumentParser(description='Description goes here.')
    parser.add_argument('n', type=int)
    args = parser.parse_args()
    n = args.n
    f = fib()
    for i,j in zip(range(n+1), f):
        print('{} {}'.format(i, j))

if __name__ == '__main__':
    main()

#!/usr/bin/env python3.4
"""
Filename:     PPER.py
Date Created: 2014-12-30 19:43
Author:       dshea <danshea@iastate.edu>
Description:  Given: Positive integers n and k such that 100≥n>0 and 10≥k>0.
              Return: The total number of partial permutations P(n,k),
              modulo 1,000,000.
"""
import argparse
import sys
import math
def main():
    """ The main() """
    parser = argparse.ArgumentParser(description='Description goes here.')
    parser.add_argument("n", type=int)
    parser.add_argument("k", type=int)
    args = parser.parse_args()
    n = args.n
    k = args.k
    print('{}'.format((math.factorial(n)/math.factorial(n-k)) % 1000000))

if __name__ == '__main__':
    main()

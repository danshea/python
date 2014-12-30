#!/usr/bin/env python3.4
"""
Filename:     INS.py
Date Created: 2014-12-31 03:00
Author:       dshea <danshea@iastate.edu>
Description:  Given:  A positive integer n≤10**3 and an array A[1..n] of
                      integers.
              Return: The number of swaps performed by insertion sort algorithm
                      on A[1..n].
"""
import argparse
import sys

def insertionSort(A):
    swap_count = 0
    for i in range(1,len(A)):
        k = i
        while k > 0 and A[k] < A[k-1]:
            A[k-1],A[k] = A[k],A[k-1]
            swap_count += 1
            k = k - 1
    return(swap_count, A)

def main():
    """ The main() """
    parser = argparse.ArgumentParser(description='Description goes here.')
    parser.add_argument('-i', '--infile', nargs='?',
                        type=argparse.FileType('r'), default=sys.stdin)
    args = parser.parse_args()
    n = int(args.infile.readline().strip())
    A = [int(i) for i in args.infile.readline().strip().split()]
    swap_count, sortedA = insertionSort(A[:n+1])
    print(swap_count)
    print(sortedA)

if __name__ == '__main__':
    main()

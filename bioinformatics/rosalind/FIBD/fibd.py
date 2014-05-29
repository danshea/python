#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
#
# Name: fibd.py
# Date: 2014-04-26
# Author: dshea
# Description: Mortal Fibonacci Rabbits
#
# Problem
#
# Figure 4. A figure illustrating the propagation of Fibonacci's rabbits if they
# die after three months.
#
# Recall the definition of the Fibonacci numbers from “Rabbits and Recurrence
# Relations”, which followed the recurrence relation Fn=Fn−1+Fn−2 and assumed
# that each pair of rabbits reaches maturity in one month and produces a single
# pair of offspring (one male, one female) each subsequent month.
#
# Our aim is to somehow modify this recurrence relation to achieve a dynamic
# programming solution in the case that all rabbits die out after a fixed number
# of months. See Figure 4 for a depiction of a rabbit tree in which rabbits live
# for three months (meaning that they reproduce only twice before dying).
#
# Given: Positive integers n≤100 and m≤20.
#
# Return: The total number of pairs of rabbits that will remain after the n-th
# month if all rabbits live for m months. Sample Dataset
#
# 6 3
#
# Sample Output
#
# 4
#
################################################################################

import sys

def fibd(n, m):
    '''
    B[n] = sum(b[n-1]...b[n-m-1])
    b[n] = B[n-1]
    F[n] = B[n] + b[n]
    '''
    B = [0 for i in range(n+1)]
    b = [0 for i in range(n+1)]
    F = [0 for i in range(n+1)]
    b[1] = 1
    F[1] = B[1] + b[1]
    for i in range(2,n+1):
            #print b, B, F
            b[i] = B[i-1]
            if i-(m-1) < 0:
                B[i] = sum(b[0:i])
            else:
                B[i] = sum(b[i-(m-1):i])
            F[i] = b[i] + B[i]
    return(F[n])

def main():
    if len(sys.argv) != 3:
        print 'usage: {0:s} n m'.format(sys.argv[0])
        sys.exit(1)
    else:
        n = int(sys.argv[1])
        m = int(sys.argv[2])
        print fibd(n, m)
    sys.exit(0)

if __name__ == '__main__':
    main()


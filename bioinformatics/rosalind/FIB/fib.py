#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
#
# Name: fib.py
# Date: 2014-04-26
# Author: dshea
# Description: Rabbits and Recurrence Relations
#
# Problem
#
# A sequence is an ordered collection of objects (usually numbers), which are
# allowed to repeat. Sequences can be finite or infinite. Two examples are the
# finite sequence (π,−2‾‾√,0,π) and the infinite sequence of odd numbers
# (1,3,5,7,9,…). We use the notation an to represent the n-th term of a sequence.
#
# A recurrence relation is a way of defining the terms of a sequence with respect
# to the values of previous terms. In the case of Fibonacci's rabbits from the
# introduction, any given month will contain the rabbits that were alive the
# previous month, plus any new offspring. A key observation is that the number of
# offspring in any month is equal to the number of rabbits that were alive two
# months prior. As a result, if Fn represents the number of rabbit pairs alive
# after the n-th month, then we obtain the Fibonacci sequence having terms Fn that
# are defined by the recurrence relation Fn=Fn−1+Fn−2 (with F1=F2=1 to initiate
# the sequence). Although the sequence bears Fibonacci's name, it was known to
# Indian mathematicians over two millennia ago.
#
# When finding the n-th term of a sequence defined by a recurrence relation, we
# can simply use the recurrence relation to generate terms for progressively
# larger values of n. This problem introduces us to the computational technique of
# dynamic programming, which successively builds up solutions by using the answers
# to smaller cases.
#
# Given: Positive integers n≤40 and k≤5.
#
# Return: The total number of rabbit pairs that will be present after n months if
# we begin with 1 pair and in each generation, every pair of reproduction-age
# rabbits produces a litter of k rabbit pairs (instead of only 1 pair).
# Sample Dataset
#
# 5 3
#
# Sample Output
#
# 19
#
################################################################################

import sys

def fib(n , k):
    if n < 0:
        return(0)
    if n == 0:
        return(1)
    if n == 1:
        return(1)
    else:
        return(fib(n-1, k)+ k*fib(n-2, k))

def main():
    if len(sys.argv) != 3:
        print 'usage: {0:s} n k'.format(sys.argv[0])
        sys.exit(1)
    else:
        n = int(sys.argv[1]) - 1
        k = int(sys.argv[2])
        print fib(n, k)
    sys.exit(0)

if __name__ == '__main__':
    main()


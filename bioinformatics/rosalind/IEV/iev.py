#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
# Name: iev.py
# Date: 2014-04-28
# Author: Dan Shea <shea.d@husky.neu.edu>
# Description:
#
# Problem
#
# For a random variable X taking integer values between 1 and n, the expected
# value of X is E(X)=∑nk=1k×Pr(X=k). The expected value offers us a way of taking
# the long-term average of a random variable over a large number of trials.
#
# As a motivating example, let X be the number on a six-sided die. Over a large
# number of rolls, we should expect to obtain an average of 3.5 on the die (even
# though it's not possible to roll a 3.5). The formula for expected value confirms
# that E(X)=∑6k=1k×Pr(X=k)=3.5.
#
# More generally, a random variable for which every one of a number of equally
# spaced outcomes has the same probability is called a uniform random variable (in
# the die example, this "equal spacing" is equal to 1). We can generalize our die
# example to find that if X is a uniform random variable with minimum possible
# value a and maximum possible value b, then E(X)=a+b2. You may also wish to
# verify that for the dice example, if Y is the random variable associated with
# the outcome of a second die roll, then E(X+Y)=7.
#
# Given: Six positive integers, each of which does not exceed 20,000. The integers
# correspond to the number of couples in a population possessing each genotype
# pairing for a given factor. In order, the six given integers represent the
# number of couples having the following genotypes:
#
#    AA-AA AA-Aa AA-aa Aa-Aa Aa-aa aa-aa
#    1.00  1.00  1.00  0.75  0.50  0.00
# Return: The expected number of offspring displaying the dominant phenotype in
# the next generation, under the assumption that every couple has exactly two
# offspring. Sample Dataset
#
# 1 0 0 1 0 1
#
# Sample Output
#
# 3.5
#
###############################################################################

import pandas
import sys

def calculate_probability(phenotypes):
    probabilities = pandas.Series([1.0, 1.0, 1.0, 0.75, 0.50, 0.00])
    return(sum(2*probabilities*phenotypes))

def main():
    if len(sys.argv) != 2:
        print 'usage {0:s} infile'.format(sys.argv[0])
        sys.exit(1)
    with open(sys.argv[1], 'r') as fh:
        phenotypes = pandas.Series([float(i) for i in fh.readline().split()])
        print calculate_probability(phenotypes)
    sys.exit(0)

if __name__ == '__main__':
    main()
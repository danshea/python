#!/usr/bin/env python

'''
Author: Dan Shea
Date: 2013-11-08
Description: Illustrative examples of computing sequence similarity

hamming_distance(a, b):
The Hamming distance between two strings of equal length is the number of
positions at which the corresponding symbols are different. 

levenshtein_distance(a, b):
The Levenshtein distance is a string metric for measuring the difference between
two sequences. Informally, the Levenshtein distance between two words is the
minimum number of single-character edits (insertion, deletion, substitution)
required to change one word into the other.

Computing the Levenshtein distance is based on the observation that if we
reserve a matrix to hold the Levenshtein distances between all prefixes of the
first string and all prefixes of the second, then we can compute the values in
the matrix in a dynamic programming fashion, and thus find the distance between
the two full strings as the last value computed. It turns out that only two rows
of the table are needed for the construction: the previous row and the current
row (the one being calculated). The Levenshtein distance may be calculated
iteratively using the following algorithm

Example: Saturday and Sunday have a levenshtein-distance of 3
(The lower right matrix value is the answer)

                    __Saturday
                    _012345678
                    S101234567
                    u211223456
                    n322233456
                    d433334345
                    a543444434
                    y654455543
'''

def hamming_distance(a, b):
    '''hamming_distance(a, b): return the hamming distance between two sequences of equal length'''
    # Are the two sequences a and b of equal length?
    if len(a) != len(b):
        raise ValueError('Arguments a and b must be the same length')
    else:
        return(sum([0 if (ord(i)^ord(j) == 0) else 1 for (i,j) in zip(list(a),list(b))]))
    
def levenshtein_distance(a, b):
    '''levenshtein_distance(a, b): compute the levenshtein-distance between two sequences'''
    # if a is longer than b, swap them
    if len(a) > len(b):
        a,b = b,a
    # degenerate cases
    if (a == b):
        return 0
    if (len(a) == 0):
        return len(b)
    if (len(b) == 0):
        return len(a)
    
    # create two working lists of integer distances
    # initialize list_0 as the edit distance for an empty a
    list_0 = [i for i in xrange(len(b)+1)]
    # initialize list_1 to be the same size as list_0
    list_1 = [None for i in xrange(len(b)+1)]
    
    for i in xrange(len(a)):
        # edit distance is delete (i+1) chars from a to match empty b
        list_1[i] = i+1
        #use formula to fill in the rest of the row
        for j in xrange(len(b)):
            cost = 0 if a[i] == b[j] else 1
            list_1[j+1] = min(list_1[j]+1, list_0[j+1]+1, list_0[j]+cost)
        # copy list_1 (current row) to list_0 (previous row) for next iteration
        list_0 = list_1[:]
    return(list_1[len(b)])

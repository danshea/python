#!/usr/bin/env python

import collections
import itertools

'''
Name: PyBLAST
Date: 2013-12-15
Desc: Me having a go at implementing BLAST in python, as I can't find a pure
      python implementation of the algorithm and I want to really understand
      what is going on at a low level.
Author: Dan Shea daniel_shea2@hms.harvard.edu
Notes:
I am a graduate student of Bioinformatics at Northeastern University.
I make no claims as to the reliability of this program with respect to the
actual NCBI BLAST programs.  This is a learning exercise for myself and
hopefully can be used pedagogically to illustrate to others what blast is
doing.  Use at your own risk!  If you're trying to do real work I would
recommend using NCBI BLAST.
'''

# Overview of what BLAST does
'''
 1. Remove low-complexity region or sequence repeats in the query sequence.
 2. Make a k-letter word list of the query sequence.
 3. List the possible matching words.
 4. Organize the remaining high-scoring words into an efficient search tree.
 5. Repeat steps 3 and 4 for each k-letter word in the query sequence.
 6. Scan the database sequences for exact matches with the remaining high-scoring words.
 7. Extend the matches to high-scoring segment pair (HSP)
 8. List all the HSPs in the database whose score is high enough to be considered.
 9. Evaluate the significance of the HSP score.
10. Make two or more HSP regions into a longer alignment
11. Show the gapped Smith-Waterman local alignments of the query and each matched db sequence.
12. Report every match whose expect score is lower than a threshold parameter E.
'''
# Entries for the BLOSUM62 matrix at a scale of ln(2)/2.0.
#BLOSUM62 = collections.OrderedDict({
#    'A':[4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,-2,-1,-1,-1],
#    'R':[-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3,-1,-2,0,-1],
#    'N':[-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3,4,-3,0,-1],
#    'D':[-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3,4,-3,1,-1],
#    'C':[0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-1,-3,-1],
#    'Q':[-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2,0,-2,4,-1],
#    'E':[-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2,1,-3,4,-1],
#    'G':[0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3,-1,-4,-2,-1],
#    'H':[-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3,0,-3,0,-1],
#    'I':[-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3,-3,3,-3,-1],
#    'L':[-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1,-4,3,-3,-1],
#    'K':[-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2,0,-3,1,-1],
#    'M':[-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1,-3,2,-1,-1],
#    'F':[-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1,-3,0,-3,-1],
#    'P':[-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2,-2,-3,-1,-1],
#    'S':[1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2,0,-2,0,-1],
#    'T':[0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0,-1,-1,-1,-1],
#    'W':[-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3,-4,-2,-2,-1],
#    'Y':[-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1,-3,-1,-2,-1],
#    'V':[0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4,-3,2,-2,-1],
#    'B':[-2,-1,4,4,-3,0,1,-1,0,-3,-4,0,-3,-3,-2,0,-1,-4,-3,-3,4,-3,0,-1],
#    'J':[-1,-2,-3,-3,-1,-2,-3,-4,-3,3,3,-3,2,0,-3,-2,-1,-2,-1,2,-3,3,-3,-1],
#    'Z':[-1,0,0,1,-3,4,4,-2,0,-3,-3,1,-1,-3,-1,0,-1,-2,-2,-2,0,-3,4,-1],
#    'X':[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
#})

BLOSUM62 = collections.OrderedDict({
    'A': collections.OrderedDict({'A': 4, 'R': -1, 'N': -2, 'D': -2, 'C': 0, 'Q': -1, 'E': -1, 'G': 0, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 0, 'W': -3, 'Y': -2, 'V': 0, 'B': -2, 'J': -1, 'Z': -1, 'X': -1}),
    'R': collections.OrderedDict({'A': -1, 'R': 5, 'N': 0, 'D': -2, 'C': -3, 'Q': 1, 'E': 0, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 2, 'M': -1, 'F': -3, 'P': -2, 'S': -1, 'T': -1, 'W': -3, 'Y': -2, 'V': -3, 'B': -1, 'J': -2, 'Z': 0, 'X': -1}),
    'N': collections.OrderedDict({'A': -2, 'R': 0, 'N': 6, 'D': 1, 'C': -3, 'Q': 0, 'E': 0, 'G': 0, 'H': 1, 'I': -3, 'L': -3, 'K': 0, 'M': -2, 'F': -3, 'P': -2, 'S': 1, 'T': 0, 'W': -4, 'Y': -2, 'V': -3, 'B': 4, 'J': -3, 'Z': 0, 'X': -1}), 
    'D': collections.OrderedDict({'A': -2, 'R': -2, 'N': 1, 'D': 6, 'C': -3, 'Q': 0, 'E': 2, 'G': -1, 'H': -1, 'I': -3, 'L': -4, 'K': -1, 'M': -3, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': 4, 'J': -3, 'Z': 1, 'X': -1}),
    'C': collections.OrderedDict({'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': 9, 'Q': -3, 'E': -4, 'G': -3, 'H': -3, 'I': -1, 'L': -1, 'K': -3, 'M': -1, 'F': -2, 'P': -3, 'S': -1, 'T': -1, 'W': -2, 'Y': -2, 'V': -1, 'B': -3, 'J': -1, 'Z': -3, 'X': -1}),
    'Q': collections.OrderedDict({'A': -1, 'R': 1, 'N': 0, 'D': 0, 'C': -3, 'Q': 5, 'E': 2, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 1, 'M': 0, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -2, 'Y': -1, 'V': -2, 'B': 0, 'J': -2, 'Z': 4, 'X': -1}),
    'E': collections.OrderedDict({'A': -1, 'R': 0, 'N': 0, 'D': 2, 'C': -4, 'Q': 2, 'E': 5, 'G': -2, 'H': 0, 'I': -3, 'L': -3, 'K': 1, 'M': -2, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 1, 'J': -3, 'Z': 4, 'X': -1}),
    'G': collections.OrderedDict({'A': 0, 'R': -2, 'N': 0, 'D': -1, 'C': -3, 'Q': -2, 'E': -2, 'G': 6, 'H': -2, 'I': -4, 'L': -4, 'K': -2, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -2, 'W': -2, 'Y': -3, 'V': -3, 'B': -1, 'J': -4, 'Z': -2, 'X': -1}),
    'H': collections.OrderedDict({'A': -2, 'R': 0, 'N': 1, 'D': -1, 'C': -3, 'Q': 0, 'E': 0, 'G': -2, 'H': 8, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -1, 'P': -2, 'S': -1, 'T': -2, 'W': -2, 'Y': 2, 'V': -3, 'B': 0, 'J': -3, 'Z': 0, 'X': -1}),
    'I': collections.OrderedDict({'A': -1, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -3, 'E': -3, 'G': -4, 'H': -3, 'I': 4, 'L': 2, 'K': -3, 'M': 1, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -3, 'Y': -1, 'V': 3, 'B': -3, 'J': 3, 'Z': -3, 'X': -1}),
    'L': collections.OrderedDict({'A': -1, 'R': -2, 'N': -3, 'D': -4, 'C': -1, 'Q': -2, 'E': -3, 'G': -4, 'H': -3, 'I': 2, 'L': 4, 'K': -2, 'M': 2, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -2, 'Y': -1, 'V': 1, 'B': -4, 'J': 3, 'Z': -3, 'X': -1}),
    'K': collections.OrderedDict({'A': -1, 'R': 2, 'N': 0, 'D': -1, 'C': -3, 'Q': 1, 'E': 1, 'G': -2, 'H': -1, 'I': -3, 'L': -2, 'K': 5, 'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 0, 'J': -3, 'Z': 1, 'X': -1}),
    'M': collections.OrderedDict({'A': -1, 'R': -1, 'N': -2, 'D': -3, 'C': -1, 'Q': 0, 'E': -2, 'G': -3, 'H': -2, 'I': 1, 'L': 2, 'K': -1, 'M': 5, 'F': 0, 'P': -2, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V': 1, 'B': -3, 'J': 2, 'Z': -1, 'X': -1}),
    'F': collections.OrderedDict({'A': -2, 'R': -3, 'N': -3, 'D': -3, 'C': -2, 'Q': -3, 'E': -3, 'G': -3, 'H': -1, 'I': 0, 'L': 0, 'K': -3, 'M': 0, 'F': 6, 'P': -4, 'S': -2, 'T': -2, 'W': 1, 'Y': 3, 'V': -1, 'B': -3, 'J': 0, 'Z': -3, 'X': -1}),
    'P': collections.OrderedDict({'A': -1, 'R': -2, 'N': -2, 'D': -1, 'C': -3, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -4, 'P': 7, 'S': -1, 'T': -1, 'W': -4, 'Y': -3, 'V': -2, 'B': -2, 'J': -3, 'Z': -1, 'X': -1}),
    'S': collections.OrderedDict({'A': 1, 'R': -1, 'N': 1, 'D': 0, 'C': -1, 'Q': 0, 'E': 0, 'G': 0, 'H': -1, 'I': -2, 'L': -2, 'K': 0, 'M': -1, 'F': -2, 'P': -1, 'S': 4, 'T': 1, 'W': -3, 'Y': -2, 'V': -2, 'B': 0, 'J': -2, 'Z': 0, 'X': -1}),
    'T': collections.OrderedDict({'A': 0, 'R': -1, 'N': 0, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 5, 'W': -2, 'Y': -2, 'V': 0, 'B': -1, 'J': -1, 'Z': -1, 'X': -1}),
    'W': collections.OrderedDict({'A': -3, 'R': -3, 'N': -4, 'D': -4, 'C': -2, 'Q': -2, 'E': -3, 'G': -2, 'H': -2, 'I': -3, 'L': -2, 'K': -3, 'M': -1, 'F': 1, 'P': -4, 'S': -3, 'T': -2, 'W': 11, 'Y': 2, 'V': -3, 'B': -4, 'J': -2, 'Z': -2, 'X': -1}),
    'Y': collections.OrderedDict({'A': -2, 'R': -2, 'N': -2, 'D': -3, 'C': -2, 'Q': -1, 'E': -2, 'G': -3, 'H': 2, 'I': -1, 'L': -1, 'K': -2, 'M': -1, 'F': 3, 'P': -3, 'S': -2, 'T': -2, 'W': 2, 'Y': 7, 'V': -1, 'B': -3, 'J': -1, 'Z': -2, 'X': -1}),
    'V': collections.OrderedDict({'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -2, 'E': -2, 'G': -3, 'H': -3, 'I': 3, 'L': 1, 'K': -2, 'M': 1, 'F': -1, 'P': -2, 'S': -2, 'T': 0, 'W': -3, 'Y': -1, 'V': 4, 'B': -3, 'J': 2, 'Z': -2, 'X': -1}),
    'B': collections.OrderedDict({'A': -2, 'R': -1, 'N': 4, 'D': 4, 'C': -3, 'Q': 0, 'E': 1, 'G': -1, 'H': 0, 'I': -3, 'L': -4, 'K': 0, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': 4, 'J': -3, 'Z': 0, 'X': -1}),
    'J': collections.OrderedDict({'A': -1, 'R': -2, 'N': -3, 'D': -3, 'C': -1, 'Q': -2, 'E': -3, 'G': -4, 'H': -3, 'I': 3, 'L': 3, 'K': -3, 'M': 2, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -2, 'Y': -1, 'V': 2, 'B': -3, 'J': 3, 'Z': -3, 'X': -1}),
    'Z': collections.OrderedDict({'A': -1, 'R': 0, 'N': 0, 'D': 1, 'C': -3, 'Q': 4, 'E': 4, 'G': -2, 'H': 0, 'I': -3, 'L': -3, 'K': 1, 'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -2, 'Y': -2, 'V': -2, 'B': 0, 'J': -3, 'Z': 4, 'X': -1}),
    'X': collections.OrderedDict({'A': -1, 'R': -1, 'N': -1, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -1, 'H': -1, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -1, 'P': -1, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V': -1, 'B': -1, 'J': -1, 'Z': -1, 'X': -1}),
})

AMINO_ACIDS = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','J','Z','X']

def remove_low_complexity_regions(query_sequence):
    '''
    remove_low_complexity_regions(query_sequence): mask out regions in the query sequence
    of low complexity
    '''
    pass

def make_kmer_list(query_sequence, k=3):
    '''
    make_kmer_list(query_sequence): create list of (N-k)+1 words for a sequence of length N
    and a word size of k
    '''
    index = 0
    word_list = list()
    while (index <= len(query_sequence)-k):
        word_list.append((index, index+k, query_sequence[index:index+k]))
        index += 1
    return word_list

def compute_neighbor_words(query_sequence, scoring_matrix=BLOSUM62, k=3, t=15):
    '''
    compute_neighbor_words(query_sequence, scoring_matrix, k, t):
    Given a query_sequence create a word_list of all kmers of length k,
    compute possible nieghbor words that score above threshold t.
    Return word_list and dictionary of kmers containing dictionary of neighbor
    words with score that meets the threshold.
    '''
    # Generate the word_list of kmers of length k
    word_list = make_kmer_list(query_sequence,k)    
    
    # Generate a list of all possible permutations of Amino Acids of length k
    unique_permutations = list(set(list(itertools.permutations(AMINO_ACIDS*3, k))))
    
    # Store results
    neighbor_words = collections.OrderedDict()
    
    # for every word in the word list calculate the score with respect to possible
    # permutations and add it to the results if the score is above threshold.
    for (start, stop, word) in word_list:
        key = ''.join(word)
        for permutation in unique_permutations:
            score = sum([BLOSUM62[a][b] for a, b in zip(word, permutation)])
            if  score >= t:
                pkey = ''.join(permutation)
                if neighbor_words.has_key(key):
                    neighbor_words[key][pkey] = score
                else:
                    neighbor_words[key] = {pkey: score}
    return (word_list, neighbor_words)
    
def compare_kmers_to_target_sequence(query_sequence, target_sequence, scoring_matrix=BLOSUM62, k=3, t=15):
    '''
    compare_kmers_to_target_sequence(query_sequence, target_sequence, k):
    given two sequences, a query and a target compare k-mers of length k from the query
    sequence to matches in the target sequence.  Return word_list, neigh_words and
    matches as a dict of dicts pointing to lists of matched k-mers each element is a pair of tuples
    consisting of ((query slice indices), (target slice indices))
    '''
    # Get the word list from the query_sequence and compute the neighbor_words
    word_list, neighbor_words = compute_neighbor_words(query_sequence, scoring_matrix, k, t)
    
    matches = collections.OrderedDict()
    for (start, stop, word) in word_list:
        index = 0
        key = ''.join(word)
        while (index <= len(target_sequence)-k):
            tkey = ''.join(target_sequence[index:index+k])
            if neighbor_words[word].has_key(tkey):
                if matches.has_key(key) and matches[key].has_key(tkey):
                    matches[key][tkey].append(((start,stop),(index,index+k)))
                else:
                    matches[key] = collections.OrderedDict({tkey: list(((start,stop),(index,index+k)))})
            index += 1
    return (word_list, neighbor_words, matches)

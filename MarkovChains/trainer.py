#!/usr/bin/env python

'''
Author: djs
Date: 2012-07-03
Description: Using the text of Ulysses by James Joyce attempt to train a Markov
chain to construct sentences.  The results have been, shall we say, interesting?

>>> construct_sentence('This','world',M)
'This thing! Poor collar. Roman of purse silken where and say: hog, you Do all
are they oysters of fire Looking lines. Hard knock. short very the dipped were
softly he can There cigarette. a with be wouldnt it and her through coming her
Flowers Lethargy. climate. the person, in light one, one, bright the for lard
No Keyes. on he biscuit the all in prescribed quantity the then But Instinct.
Cream. cones.'
>>> construct_sentence('This','world',M)
"This mouth. his to going he's and in Harrison's of traditions metaphysical the
wielding son, my In came. he where busy perch, the by them knows me. kissed
open, orifice rectal for lie than rational less increasingly purchase of boys
the land.)_"

>>> construct_sentence('Rock','castle',M)
'Rock Three Israel. for Good possible. as shallow as ignorant Ha ho! Ho
_(Laughing)_ WHORES: THREE'

Obviously the algorithm needs to be fine tuned to strip offending punctuation
and the like, but I figured I'd give it a few runs to see what it comes up with.
'''

import Markov

def train(filename, freq = dict()):
    '''
    read the contents of the file and store the state transitions from word to
    word.  Pass this information on to the analyzer to createa Markov chain for
    generating sentences.  Can be run on multiple files to train the Markov
    chain before passing to analyze to create the Markov object.
    '''
    with open(filename, 'r') as fh:
        for line in fh:
            words = line.strip().split()
            if words == []:
                continue
            prev = words.pop(0)
            while words:
                curr = words.pop()
                if prev in freq:
                    if curr in freq[prev]:
                        freq[prev][curr] += 1.0
                    else:
                        freq[prev][curr] = 1.0
                else:
                    freq[prev] = dict()
                    freq[prev][curr] = 1.0
                prev = curr
    return freq

def analyze(freq):
    for k in freq:
        count = sum([i for i in freq[k].values()])
        for k2 in freq[k]:
            freq[k][k2] = freq[k][k2] / count
    return freq

def construct_sentence(start, stop, markov):
    sentence = start
    current_word = start
    while current_word != stop:
        try:
            next_word = markov.next(current_word)
            sentence += ' ' + next_word
            current_word = next_word
        except KeyError:
            break
    return sentence

if __name__ == '__main__':
    freq = dict()
    filenames = ['Ulysses.txt']
    for f in filenames:
        freq = train(f, freq)
    M = Markov.Markov(analyze, freq)
    
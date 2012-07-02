#!/usr/bin/env python

'''
Author: djs
Date: 2012-07-02
Description: Using Markov chains attempt to analyze the frequencies a letter
will transition from one to the next in the English language, using this
information, attempt to generate words using the given transition table.
'''

import random

def analyze(filename, freq=dict()):
    '''
    Perform statistical analysis on the file returning a dictionary of the
    probabilities of state transitions from one character to another.
    '''
    with open(filename, 'r') as fh:
        prev = fh.read(1)
        if not prev:
            return freq
        count = 1
        while True:
            c = fh.read(1)
            if not c:
                break
            if freq.has_key(prev):
                if freq[prev].has_key(c):
                    freq[prev][c] += 1.0
                else:
                    freq[prev][c] = 1.0
            else:
                freq[prev] = dict()
                freq[prev][c] = 1.0
            count += 1
            prev = c
        for k in freq:
            total = sum([i for i in freq[k].values()])
            for k2 in freq[k]:
                freq[k][k2] = freq[k][k2] / total
    return freq

def sanity_check(freq):
    '''
    This function should return a value close to 1.0 as it is the sum of all the
    probabilities in the dict of dicts, that makes up the transitions table.
    '''
    probability_totals = [sum(v.values()) for v in [val for val in freq.values()]]
    return sum(probability_totals)/len(probability_totals)

def build_transition_table(freq):
    '''
    Given a dictionary of dictionaries detailing the probability of transitions
    from one state to another, construct the transition dictionary.
    '''
    t_table = dict()
    for k in freq:
        t_table[k] = dict()
        trans = zip(freq[k].keys(), freq[k].values())
        prev = trans.pop()
        t_table[k][prev[0]] = (0, prev[1])
        prev = t_table[k][prev[0]]
        while trans:
            curr = trans.pop()
            t_table[k][curr[0]] = (prev[1], prev[1] + curr[1])
            prev = t_table[k][curr[0]]
    return t_table

def next_state(state, t_table):
    '''
    Given a state and a table of transition probabilities, compute the next
    state.
    '''
    partitions = [(k, t_table[state][k][0], t_table[state][k][1]) for k in t_table[state]]
    rval = random.random()
    for (next_state, lower_bound, upper_bound) in partitions:
        if lower_bound < rval <= upper_bound:
            return next_state
    
def build_word(start, stop, t_table):
    '''
    Given a start state, a stop state and a table of transition probabilities,
    construct a word and output that word.
    '''
    state = start
    word = start
    while state != stop:
        state = next_state(state, t_table)
        word += state
    return word
#!/usr/bin/env python

'''
Author: djs
Date: 2012-07-02
Description: Encapsulate the idea of a Markov chain into a class.  Pass the class
a function that analyzes some data set and returns a dictionary of dictionaries
that details the state transitions and the associated probabilities.  Using that
information construct a transition table and use it to calculate how the chain
would move from one state to the next based on the given probabilities.
'''

import random

class Markov(object):
    '''
    Markov class requires a function that returns a dictionary of dictionaries
    that details the probability of transition from one state to another.
    {'a': {'a': 0.25, 'b': 0.25, 'c': 0.25, 'd': 0.25}}
    '''
    def __init__(self, analyzer, *args, **kwargs):
        self.analyzer = analyzer
        self.frequencies = self.analyzer(*args, **kwargs)
        self.transitions = dict()
        self.build_transition_table()
        with open('/dev/random', 'r') as fh:
            random.seed(fh.read(4))
        
    def build_transition_table(self):
        '''
        Construct the transition dictionary.
        '''
        for k in self.frequencies:
            self.transitions[k] = dict()
            trans = zip(self.frequencies[k].keys(), self.frequencies[k].values())
            prev = trans.pop()
            self.transitions[k][prev[0]] = (0, prev[1])
            prev = self.transitions[k][prev[0]]
            while trans:
                curr = trans.pop()
                self.transitions[k][curr[0]] = (prev[1], prev[1] + curr[1])
                prev = self.transitions[k][curr[0]]

    def next(self, state):
        '''
        Given the current state, compute the next state.
        '''
        partitions = [(k, self.transitions[state][k][0], self.transitions[state][k][1]) for k in self.transitions[state]]
        rval = random.random()
        for (next_state, lower_bound, upper_bound) in partitions:
            if lower_bound < rval <= upper_bound:
                return next_state
    
if __name__ == '__main__':
    '''
    Example of a very basic analyzer function that just returns a dictionary of
    dictionaries.
    table = fn()
    table['a']['a'] = 0.25, meaning there is a 25% probably of transitioning
    from state 'a' to state 'a'.
    '''
    def fn():
        return {'a': {'a': 0.25, 'b': 0.25, 'c': 0.25, 'd': 0.25},
                'b': {'a': 0.25, 'b': 0.25, 'c': 0.25, 'd': 0.25},
                'c': {'a': 0.25, 'b': 0.25, 'c': 0.25, 'd': 0.25},
                'd': {'a': 0.25, 'b': 0.25, 'c': 0.25, 'd': 0.25},}
    M = Markov(fn)
    
    
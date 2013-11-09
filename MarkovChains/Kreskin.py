#!/usr/bin/env python

'''
Author: djs
Date: 2012-07-02
Description: Quick analysis of historical market data using Markov chains.
Interesting to note after analyzing several stocks going back 10 years the
statistics in the Markov chain for a given stock show the likelihood of moving
up or down in price based on the previous days close only is close to chance.
There is a slight upward pressure on the market as evidenced by the probabilities
given in the frequencies calculated.
'''

import Markov
import urllib2
import csv

def analyzer(symbol='IBM'):
    '''
    Grab historical pricing information from Yahoo Finance for a given ticker
    symbol.
    '''
    # Set up arguments to pass to the Yahoo Finance Website
    # In our case we will examine 10 years of historical data
    kwargs = { 'symbol': symbol,
               'start_month': '01',
               'start_day': '01',
               'start_year': '2002',
               'end_month': '12',
               'end_day': '31',
               'end_year': '2012',
             }
    # Format the url
    urlstring = 'http://ichart.finance.yahoo.com/table.csv?s={symbol}&a={start_month}&b={start_day}&c={start_year}&d={end_month}&e={end_day}&f={end_year}&g=d&ignore=.csv'.format(**kwargs)
    # grab the data as a csv and place it into a list
    data = [row for row in csv.reader(urllib2.urlopen(urlstring))]
    # pop off the csv column headings
    headings = data.pop(0)
    
    # initialize the frequencies dictionary
    freq = {'+': {'+': 0.0, '-': 0.0}, '-': {'+': 0.0, '-': 0.0}}
    
    # Analyze the data
    total = len(data)
    prev = data.pop()
    freq['-']['-'] += 1
    up = False
    while data:
        curr = data.pop()
        if prev[-1] < curr[-1] and up:
            freq['+']['+'] += 1.0
            up = True
        elif prev[-1] > curr[-1] and up:
            freq['+']['-'] += 1.0
            up = False
        elif prev[-1] < curr[-1]:
            freq['-']['+'] += 1.0
            up = True
        else:
            freq['-']['-'] += 1.0
            up = False
        prev = curr

    for k in freq:
            count = sum([i for i in freq[k].values()])
            for k2 in freq[k]:
                freq[k][k2] = freq[k][k2] / count
    
    return freq

if __name__ == '__main__':
    mchains = [(symbol, Markov.Markov(analyzer, symbol=symbol)) for symbol in ['IBM','YHOO','GOOG']]
    
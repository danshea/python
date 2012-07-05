#!/usr/bin/env python

import csv
import datetime
from matplotlib import pyplot as plt
import urllib2

def fetch_data(symbol='IBM'):
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
    return (headings, data)

if __name__ == '__main__':
    headings, data = fetch_data()
    plt.plot([datetime.datetime.strptime(d[0], '%Y-%m-%d') for d in data],[d[-1] for d in data])
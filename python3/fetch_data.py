#!/usr/bin/env python3

import codecs
import csv
import datetime
import urllib.request

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

    # since python3 returns bytes instead of strings, we can no longer use our previous python2 code that
    # simply tossed the response to the csv.reader, we must first map the HttpResponse to decode with the proper
    # encoding and then pass that to the csv.reader

    # It is important to note map is mapping bytes.decode to each __next__ in the
    # HttpResponse and the encoding is passed via a generator that never terminates
    # iter(int, 1) will never terminate since int() will always return 0, therefore,
    # map returns once the input HttpResponse object is exhausted.

    # data = [row for row in csv.reader(map(bytes.decode, urllib.request.urlopen(urlstring), ('iso-8859-1' for i in iter(int,1))))]
    # can be more clearly written as 
    # data = [row for row in csv.reader(map(bytes.decode, urllib.request.urlopen(urlstring), ('iso-8859-1' for i in iter(lambda:0,1))))]

    # turns out the codecs library offers us this functionality in a cleaner, more pythonic fashion
    data = [row for row in csv.reader(codecs.iterdecode(urllib.request.urlopen(urlstring), 'iso-8859-1'))]
    
    # since we define a lambda statement that always returns 0, one might not know int() always returns 0 and what if that behaviour
    # were to change in later versions of python?

    # pop off the csv column headings
    headings = data.pop(0)
    return (headings, data)

def rolling_avg(data, size):
    '''
    rolling_avg computes the rolling average for samples in data
    using a window size of size
    '''
    return [sum(data[start:end])/size for start, end in zip(range(0,len(data)-size),range(size,len(data)-size))]

if __name__ == '__main__':
    headings, data = fetch_data('IBM')
    avgs = rolling_avg([float(d[-1]) for d in sorted(data)], 30)

#!/usr/bin/env python
'''
Author: djs
Date: 2014-05-29
Description: The venerable FizzBuzz

Given a number, determine if it is divisiable by either of two other numbers
Or simply print the number if it is evenly divisble by neither number.
'''

def fizzbuzz(number, a, b):
    '''def fizzbuzz(number, a, b): print fizz if evenly divisible by a
                                   print buzz if evenly divisible by b
    '''
    divisible = False
    if number % a == 0:
        print 'Fizz',
        divisible = True
    if number % b == 0:
        print 'Buzz'
        divisible = True
    if not divisible:
        print '{0:d}'.format(number)
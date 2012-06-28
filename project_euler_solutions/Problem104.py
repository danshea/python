#!/usr/bin/env python

# Project Euler Problem 104 - http://projecteuler.net/problem=104
'''
The Fibonacci sequence is defined by the recurrence relation:

    Fn = Fn-1 + Fn-2, where F1 = 1 and F2 = 1.

It turns out that F541, which contains 113 digits, is the first Fibonacci number for which the last nine digits are 1-9 pandigital (contain all the digits 1 to 9, but not necessarily in order). And F2749, which contains 575 digits, is the first Fibonacci number for which the first nine digits are 1-9 pandigital.

Given that Fk is the first Fibonacci number for which the first nine digits AND the last nine digits are 1-9 pandigital, find k.
'''

# This fails due to overflow errors - numerical computing dictates a simpler appraoch!
# Using the closed form solution (Binet's formula) for generating fibonacci numbers we can define a function to generate the sequence
# Fn = (phi**n - psi**n) / sqrt(5.0)
#import math
#def fib(n):
#    phi = (1.0 + math.sqrt(5.0)) / 2.0 # The golden ratio
#    psi = (1.0 - math.sqrt(5.0)) / 2.0
#    return (phi**n - psi**n) / math.sqrt(5.0)

# Iterative computation, computing fib(n) in n steps, still too costly
#def fib(n):
#    a, b = 0, 1
#    for i in range(n):
#        a, b = b, a + b
#    return a

# Using Djikstra's algorithm which I am still working through his paper to understand (The man was a genius)
fibs = {0: 0, 1: 1}
def fib(n):
    if n in fibs: return fibs[n]
    if n % 2 == 0:
        fibs[n] = ((2 * fib((n / 2) - 1)) + fib(n / 2)) * fib(n / 2)
        return fibs[n]
    else:
        fibs[n] = (fib((n - 1) / 2) ** 2) + (fib((n+1) / 2) ** 2)
        return fibs[n]

# we need to examine if the first 9 and last 9 digits of the fibonacci number are pandigital, let's write a helper function to verify this
def is_pandigital(n):
    sn = str(n)
    head = [int(i) for i in sn[0:9]]
    tail = [int(i) for i in sn[-9:]]
    if sum(head) == 45 and sum(tail) == 45:
        try:
            if head.index(0) or tail.index(0):
                return False
        except ValueError:
            match = [0 for i in range(9)]
            for val in head:
                match[val-1] = 1
            if sum(match) == 9:
                match = [0 for i in range(9)]
                for val in tail:
                    match[val -1] = 1
                if sum(match) == 9:
                    return True
    return False

def main():
    # Now, since from the problem we are told the first fibonacci number with the first nine digits being pandigital is F2749
    # let's start at F2750 and work our way up via brute force
    k = 2750
    while is_pandigital(fib(k)) == False:
        k += 1
    print 'Solution: k = {0}'.format(k)

if __name__ == '__main__':
    main()
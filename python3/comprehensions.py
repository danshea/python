#!/usr/bin/env python3

'''
Author: djs
Date: 2012-07-09
Description: Brief examples of comprehensions
Comprehensions are a powerful tool for creating complex data structures
in a relatively straightforward and elegant fashion.

Below are some simple examples.
'''

# List comprehension
L = [i for i in range(10)]
'''
>>> L
[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
>>> type(L)
<class 'list'>
'''

# Generator comprehension
G = (i for i in range(10))
'''
>>> G
<generator object <genexpr> at 0xb57ac504>
>>> type(G)
<class 'generator'>
'''

# Set comprehension
S = {i for i in range(10)}
'''
>>> S
{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
>>> type(S)
<class 'set'>
'''

# Dictionary comprehension
D = {key:value for key, value in zip(range(10), range(10))}
'''
>>> D
{0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9}
>>> type(D)
<class 'dict'>
'''


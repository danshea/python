#!/usr/bin/env python
'''
Author: djs
Date: 2012-06-23
Description: Python argument unpacking, or how to write functions that take an
arbitrary number of arguments.
'''

# So let's say you have a simple function
def f(x, y):
    return x+y

'''
If you have a bunch of arguyments you would like to pass to this function as
written, you have to do something like this

>>> total=0
>>> for a, b in zip(range(0,5),range(5,10)):
...     total += f(a,b)
... 
>>> total
45

However, we can utilize argument unpacking and a bit of a rewrite to create
a function that will sum up an arbitrary number of arguments.
'''

def g(*args):
    total = 0
    for arg in args:
        total += arg
    return total

'''
>>> g(1,2,3,4,5,6,7,8,9)
45

Or less verbosely

>>> g(*range(10))
45

So why a * when using range?  In Python 2.x range returns a list.

>>> range(10)
[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

The * unpacks the list which is then passed to the function.  The function takes
arguments into it as a tuple.  So in our case we want g to be given
(1,2,3,4,5,6,7,8,9)

As in the first example.  That is the meaning behind unpacking.  If we define
*args in our function definition Python will unpack the tuple it is given.  It
is still up to use to handle how those arguments are parsed, used appropriately
inside our function, but now we can have an arbitrary number of arguments.

Python also supports this functionality for arguments that are named.  Named
arguments have an associated keyword.  They are passed in as a dictionary.

Let's take a look at how that works.
'''

def h(*args, **kwargs):
    for arg in args:
        print 'Unamed argument: {0}'.format(arg)
    for kw in kwargs:
        print 'Named  argument: {0} = {1}'.format(kw, kwargs[kw])

'''
>>> h(1,2,3,foo='bar')
Unamed argument: 1
Unamed argument: 2
Unamed argument: 3
Named  argument: foo = bar
'''


#!/usr/bin/env python

'''
Author: djs
Date: 2012-06-29
Description: Closures - In Python, functions are first class objects
So much as we utilize this concept to create decorators, we can further extend
it to allow nested functions to better exploit the power of a closure.

Values inside the nested functions are bound at definition time to the arguments
passed in, they are said to be "closed over."
'''

def tracer(func):
    def _wrap(fh):
        def _deco(*args, **kwargs):
            fh.write('Calling function {0}\n'.format(func.__name__))
            return func(*args, **kwargs)
        return _deco
    return _wrap

if __name__ == '__main__':
    import sys
    # define a function to be wrapped
    def add_one(x):
        return x + 1
    # call the function inside the closure and tell it to use sys.stdout for output
    wrapped_add_one = tracer(add_one)(sys.stdout)
    result = wrapped_add_one(1)
    print result
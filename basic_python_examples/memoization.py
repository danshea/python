#!/usr/bin/env python
'''
Author: djs
Date: 2012-07-07
Description: Example of memoization, take a function and store the arguments to
the function as a key into hash table that stores the function's return value.
If the function is called again  with the same arguments you then simply return
the pre-computed value that is now stored in the hash table.
'''

# basic example of a decorator
def memoize(fn):
    _hash = dict()
    def memoized(*args):
        if args not in _hash:
            _hash[args] = fn(*args)
        return _hash[args]
    return memoized

# example call
@memoize
def square1(x):
    return x * x

# example that allows you to optionally pass the cache in as an argument
def memoize(_hash = None):
    if _hash is None:
        _hash = dict()
    def decorator(fn):
        def decorated(*args):
            if args not in _hash:
                _hash[args] = fn(*args)
            return _hash[args]
        return decorated
    return decorator

# make square2(1) return a bogus value by pre-caching
@memoize({(1,): -256})
def square2(x):
    return x * x
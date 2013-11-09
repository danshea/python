#!/usr/bin/env python

'''
Author: djs
Date: 2012-06-23
Description: Python decorators
'''

# In Python a function can be passed around like any other object
# Here is a simple function f(x) that returns x
def f(x):
    return x
# set g equal to f
g = f

'''
So here, we see that both f and g point to the same code

>>> def f(x):
...     return x
... 
>>> f(1)
1
>>> g=f
>>> g(1)
1
'''

'''
Additionally, in Python, a function can be defined inside of another function
Here we have defined the function b inside of the function a
>>> def a(x):
...     def b(x):
...             return x
...     return x * b(x)
... 
>>> a(2)
4
>>> a(3)
9

So what does this have to do with decorators?
A decorator is a function that takes a function as an argument and modifies that
function, returning the newly modified function.
'''

# Here is a simple decorator
def all_caps(f):
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs).upper()
    return wrapper

# And here is a simple function that uses that decorator
@all_caps
def print_string(s):
    return s

'''
So we can see that the print_string function is wrapped by the all_caps decorator
all_caps defines a function wrapper and wrapper calls f with its arguments.
Since f returns a string the upper() method of the string object is called.

>>> def all_caps(f):
...     def wrapper(*args, **kwargs):
...         return f(*args, **kwargs).upper()
...     return wrapper
... 
>>> @all_caps
... def print_string(s):
...     return s
... 
>>> print_string('Hello')
'HELLO'

***Notice all_caps is returning the function wrapper, not the results of calling
upper.  This is because we can define some other function that returns a string
and use this same decorator
'''

@all_caps
def join_two(word, another_word):
    return word + another_word

'''
>>> @all_caps
... def join_two(word, another_word):
...     return word + another_word
... 
>>> join_two('foo','bar')
'FOOBAR'
'''

'''
You are probably asking yourself now, "But what is a good real world example of why I want to use this
crazy stuff?"

Debugging - The ability to wrap function calls with debugging statements
Logging, likewise, you can write your code and add the logging later
Error Handling - wrapping your function calls in try/except blocks
Connection handling - opening and closing connections by wrapping the function in a decorator that adds
the code for you

Let's look at the try/except example for a more concrete reason to use decorators
'''

def handler(func):
    def _handle_errors(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            print 'Houston we have a problem. {0} encountered.'.format(e)
    return _handle_errors

@handler
def risky_func():
    1/0
    
'''
>>> def handler(func):
...     def _handle_errors(*args, **kwargs):
...         try:
...             return func(*args, **kwargs)
...         except Exception as e:
...             print 'Houston we have a problem. {0} encountered.'.format(e)
...     return _handle_errors
... 
>>> @handler
... def risky_func():
...     1/0
... 
>>> risky_func()
Houston we have a problem. integer division or modulo by zero encountered.
'''

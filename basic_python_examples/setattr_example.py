#!/usr/bin/python

'''
Author: djs
Date: 2012-06-29
Description:  Was recently asked an interesting question regarding how we can
ensure that class attribute values are within a valid range and of a valid type.
Here is a simple example using __setattr__ to override default behaviour and to
ensure we get and use only positive integers, We raise an exception otherwise.

'''

class PositiveInteger:
    def __init__(self, i):
        self.val = i
        
    def __setattr__(self, name, value):
        # ensure value is of type int
        if type(value) == int:
            # ensure value is positive
            if value < 0:
                raise ValueError
            else:
                # set the value, must use self.__dict__[name] = value to
                # avoid a recursive call to itself (i.e. - do not use self.val = value)
                self.__dict__[name] = value
        else:
            raise ValueError
        

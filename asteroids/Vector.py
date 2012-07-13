#!/usr/bin/env python

import math

class Vector(object):
    '''Vector class.'''
    def __init__(self, r, theta):
        self.r = r
        self.theta = theta
        self.x = r * math.cos(math.radians(theta))
        self.y = r * math.sin(math.radians(theta))

    def __add__(self, v):
        x = self.x + v.x
        y = self.y + v.y
        r = math.sqrt(x**2 + y**2)
        theta = math.degrees(math.atan2(y, x))
        return Vector(r, theta)

    def __setattr__(self, name, value):
        '''
        Ensure we are dealing with floats.
        Coerce integers to float and raise ValueError otherwise.
        '''
        if name in ['r', 'theta', 'x', 'y']:
            if type(value) is int:
                super(Vector, self).__setattr__(name, float(value))
            elif type(value) is float:
                super(Vector, self).__setattr__(name, value)
            else:
                raise ValueError('{0} must be float or int'.format(name))
        else:
            raise AttributeError('{0} object has no attribute {1}'.format(Vector.__name__, name))

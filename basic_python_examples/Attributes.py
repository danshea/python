#!/usr/bin/env python

'''
Author: djs
Date: 2012-07-01
Description: Example of how __getattribute__, __getattr__ and __setattr__ overloading works

__getattribute__ is called whenever an attribute is accessed
__getattr__ is called when an attribute can not be found via the usual methods
__setattr__ is called when an attribute value is to be set or modified

Here we have a class Page which is derived from the object class

1. The object's create attribute is modified with the current time when the object is created
2. Whenever the data attribute is accessed, the object's access attribute is modified with the current time
3. Whenever the data attribute is modified, the object's modify attribute is modified with the current time
4. If someone attempts to set one of these three attributes to a non-float value, a ValueError is raised
'''

import time
class Page(object):
    def __init__(self, data=None):
        epoch = time.time()
        self.data = data
        self.create = epoch
        self.access = epoch
        self.modify = epoch
    
    def __getattribute__(self, name):
        if name == 'data':
            super(Page, self).__getattribute__('__dict__')['access'] = time.time()
        return super(Page, self).__getattribute__(name)
    
    def __getattr__(self, name):
        if name == 'data':
            epoch = time.time()
            self.__dict__['create'] = epoch
            self.__dict__['modify'] = epoch
            self.__dict__['access'] = epoch
            return self.__dict__[name]
        else:
            return self.__dict__[name]
    
    def __setattr__(self, name, value):
        if name in ['create', 'access', 'modify']:
            if type(value) is float:
                self.__dict__[name] = value
            else:
                raise ValueError
        else:
            self.__dict__['modify'] = time.time()
            self.__dict__[name] = value

    def stat(self):
        print 'Create: {0}'.format(self.create)
        print 'Access: {0}'.format(self.access)
        print 'Modify: {0}'.format(self.modify)
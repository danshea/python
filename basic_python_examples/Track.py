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

Taking this and extending it a bit, we can define a decorator tracked that wraps an object in the Track class
and stores the object in the data attribute of the Track class, we can now keep track of when the object is
accessed or modified and when the object was created.
'''
import datetime
import time

class Track(object):
    def __init__(self, data=None):
        epoch = time.time()
        self.data = data
        self.create = epoch
        self.access = epoch
        self.modify = epoch
    
    def __getattribute__(self, name):
        if name == 'data':
            super(Track, self).__getattribute__('__dict__')['access'] = time.time()
        return super(Track, self).__getattribute__(name)
    
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
        print 'Create: {0}'.format(datetime.datetime.fromtimestamp(self.create))
        print 'Access: {0}'.format(datetime.datetime.fromtimestamp(self.access))
        print 'Modify: {0}'.format(datetime.datetime.fromtimestamp(self.modify))

def tracked(fn):
    def wrap(*args, **kwargs):
        return Track(data=fn(*args, **kwargs))
    return wrap

if __name__ == '__main__':
    @tracked
    def tracked_list(*args):
        return list(*args)
    
    tlist = tracked_list(range(10))
    for i in tlist.data:
        print 'Value: {0}'.format(i)
        tlist.stat()

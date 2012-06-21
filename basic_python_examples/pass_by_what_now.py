#!/usr/bin/env python
'''
Author: djs
Date: 2012-06-22
Description: Python, pass-by-value or pass-by-reference? Answer: Yes please.

Python passes parameters by value, however some values being passed are actually
references.  To muddy the water even more, some of those references are mutable,
while others are immutable.  This leads to some confusion when dealing with
scope and whether or not the outer scope's variable had its contents modifed.

Let's go through some examples to see what this means.
'''

def add_one(x):
    print 'I was passed x={0:d}'.format(x)
    x = x + 1
    print 'Inside the function it is now x={0:d}'.format(x)

def append_one(l):
    print 'I was passed l={0:s}'.format(l)
    l.append(len(l))
    print 'Inside the function it is now l={0:s}'.format(l)

def reassign(l):
    print 'I was passed l={0:s}'.format(l)
    l = range(3,6)
    print 'Inside the function it is now l={0:s}'.format(l)

def reassign2(l):
    print 'I was passed l={0:s} with an object id of {1:d}'.format(l, id(l))
    l = range(3,6)
    print 'Inside the function it is now l={0:s} with an object id of {1:d}'.format(l, id(l))

if __name__ == '__main__':
    # First we look at basic variables
    x = 0
    print 'Initially, x={0:d}'.format(x)
    add_one(x)
    print 'After the function call x={0:d}'.format(x)
    
    # Now, we'll look at a mutable object like a list
    l = range(3)
    print 'Initially, l={0:s}'.format(l)
    append_one(l)
    print 'After the function call l={0:s}'.format(l)
    # The list has been modified!  This is because the value being passed
    # to the function is the address of the list object, if you are familiar
    # with C this is similar to how we pass by value in C a pointer to an
    # object and modify the contents of that memory address.  The address
    # is not changed, but the contents it points to have been.
    
    # Aha! but now let's see what happens when I try to reassign a list
    l=range(3)
    print 'Initially, l={0:s}'.format(l)
    reassign(l)
    print 'After the function call l={0:s}'.format(l)
    # So initially, we passed a reference to some object and then reaplced
    # the reference to that object with another one, this does not change
    # the outer scope's reference.
    
    # Let's look at the above example again, this time keeping track of
    # of the object id as well as the value.
    l=range(3)
    print 'Initially, l={0:s} and its object id is {1:d}'.format(l, id(l))
    reassign2(l)
    print 'After the function call l={0:s} and its object id is {1:d}'.format(l, id(l))
    # Notice that the object id changed inside the function when we reassigned l
    # But outside the function the reference still points to the original object
    
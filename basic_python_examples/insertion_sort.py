#!/usr/bin/env python

'''
Author: djs
Date: 2012-06-26
Description:
simple implementation of insertion sort using a second list
O(n**2)
'''

def insertion_sort(l):
    result = []
    for i in l:
        inserted = False
        for j in result:
            if i <= j:
                result.insert(result.index(j),i)
                inserted = True
                break
        if not inserted:
            result.append(i)
    return result

'''
re-written to sort in place.
'''

def inplace_isort(l):
    for i in range(1,len(l)):
        val = l.pop(i)
        inserted = False
        for j in range(i-1, -1, -1):
            if val > l[j]:
                l.insert(j+1,val)
                inserted = True
                break
        if not inserted:
            l.insert(0,val)
    return l

if __name__ == '__main__':
    for i in xrange(10000):
        inplace_isort([random.randint(0,100) for i in xrange(100)])

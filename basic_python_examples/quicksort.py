#!/usr/bin/env python

# Example implementation of quicksort
# This is laughably slower than the built in sort algorithm for the list data type and far uglier
# Purely intended as an exercise to ensure my brain didn't fully atrophy decades after university ;-)

def quicksort(l):
    '''quicksort(l) - sort a list l using the infamous quicksort algorithm'''
    # if the list is of length less than 2 you're done, so return the list
    if len(l) < 2:
        return l
    # if the list is of length 2, determine the proper order of the elements and return the list
    if len(l) == 2:
        if l[0] > l[1]:
            return [l[1],l[0]]
        else:
            return l
    # if the list is of length 3 or more, pick the pivot and recurse
    pivot = len(l) / 2
    # everything less than the pivot value goes left
    pivot_val = l[pivot]
    left = list()
    right = list()
    tmp = l[0:pivot]
    tmp.extend(l[pivot+1:])
    for i in tmp:
        if i < pivot_val:
            left.append(i)
        else:
            right.append(i)
    l=quicksort(left)
    l.append(pivot_val)
    l.extend(quicksort(right))
    return l
#!/usr/bin/env python3.4
"""
Filename:     LGIS.py
Date Created: 2015-01-13 22:37
Author:       dshea <danshea@iastate.edu>
Description:  Given:  A positive integer n≤10000 followed by a permutation π of
                      length n.
              Return: A longest increasing subsequence of π, followed by a
                      longest decreasing subsequence of π.
"""
import argparse
import sys

class Item:
    def __init__(self, val):
        self.val = val
        self.link = None

def findStack(stacks, s):
    """ Find leftmost stack where s.val < top of stack return index"""
    index = 0
    for stack in stacks:
        if stack == []:
            return(index)
        if stack[-1].val > s.val:
            return(index)
        else:
            index += 1
    return(index)

def patienceSort(n, seq):
    """ Patience sort to determine longest increasing and decreasing subsequences """
    # Initialize an array of n empty stacks
    stacks = [list() for i in range(n)]
    for s in seq:
        # Find the leftmost stack which fits Item and put it there, keeping a
        # reference to the top-most Item on the stack to the left of this stack
        index = findStack(stacks, s)
        if index != 0:
            s.link = stacks[index-1][-1]
            stacks[index].append(s)
        else:
            stacks[index].append(s)
    # Reconstruct the subsequence starting with the top-most Item on the
    # right-most occupied stack and work back from there.
    index = n-1
    while stacks[index] == []:
        index -= 1
    # First Item in lgis
    current = stacks[index][-1]
    lgis = [current]
    while current.link != None:
        current = current.link
        lgis.append(current)
    return(reversed(lgis))

def main():
    """ The main() """
    parser = argparse.ArgumentParser(description='Description goes here.')
    parser.add_argument('-i', '--infile', nargs='?',
                        type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('-o', '--outfile', nargs='?',
                        type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()

    # Print out the argument values
    #print('infile {}'.format(args.infile.name))
    #print('outfile {}'.format(args.outfile.name))

    n   = int(args.infile.readline().strip())
    seq = [Item(int(i)) for i in args.infile.readline().strip().split()]
    print('For n = {}'.format(n))
    print('The sequence {} yields:'.format([i.val for i in seq]))
    print([i.val for i in patienceSort(n, seq)])

if __name__ == '__main__':
    main()

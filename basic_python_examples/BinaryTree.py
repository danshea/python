#!/usr/bin/env python

'''
Author: djs
Date: 2012-06-21
Description: Binary trees stored as arrays

A simple binary tree, can be stored as an array.
If we observe the structure of the tree we see the following, labeling each node
in order and mapping that as the index into the array for that node:

Tree Structure  Index into an array
    0           i=0                  
   / \
  1   2         i=1 i=2              
 / \ / \
3  4 5  6       i=3 i=4 i=5 i=6

Properties:
2*i+1 will yield the left  child of the node located at index i
2*i+2 will yield the right child of the node located at index i

Left children  (i-1)/2 will yield the parent
Right children (i-2)/2 will yield the parent

In this way we can navigate a tree structure stored as a flat list of objects.
The only information we need to know about a given node is if it is a left or
a right child.  This can be determined by the index itself, as all left children
will have an odd numbered index, and all right children will have an even
numbered index.

The number of leaves of the binary tree can be determined by the property that
at a depth n, the tree will have 2**n-1 leaves.  The total size of the array
required to store a binary tree in this fashion can be computed as the
sum of 2**k k->[0,n-1] i.e.- a tree of depth 3 (n=3) requires an array of size
2**0 + 2**1 + 2**2 = 1+2+4 = 7

A breadth first traversal of a tree stored in this way is simple.  We simply
iterate through the list.

A depth first traversal of a tree stored in this fashion requires slightly
more work however because we know the total length of the array and how to
navigate both forwards (from parent->child) and backwards (from child->parent)
we can implement a function to produce the desired results.
'''

class BinaryTree:
    def __init__(self):
        self.tree = list()
        self.last = None
        self.size = len(self.tree)
    
    def _updatesize(self):
        self.size = len(self.tree)
    
    def _updatelast(self, index):
        self.last = index
    
    def children(self, i):
        '''return indices of all children of node located at i'''
        children = []
        if 2*i+2 > self.size:
            return children
        elif 2*i+2 == self.size:
            children.append(2*i+1)
            return children
        else:
            children.extend([2*i+1, 2*i+2])        
        for i in children:
            if 2*i+2 > self.size:
                return children
            elif 2*i+2 == self.size:
                children.append(2*i+1)
                return children
            else:
                children.extend([2*i+1, 2*i+2])

if __name__ == '__main__':
    btree = BinaryTree()
    btree.tree = range(7)
    btree.size = len(btree.tree)
    print btree.children(0)
    btree.tree = range(8)
    btree.size = len(btree.tree)
    print btree.children(1)

'''
Output from above given the following tree(s):

      0
     / \
    1   2
   / \ / \
  3  45   6
 /
7
>>> execfile('BinaryTree.py')
[1, 2, 3, 4, 5, 6]
[3, 4, 7]

Note the actual values could be anything, they were merely set to be equal to
the index into the array for ease of illustration.
'''

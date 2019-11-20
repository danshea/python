#!/usr/bin/env python
# Author: Dan Shea
# Date: 2019.11.20

# Simple example of binary tree in python
# I don't use tree data structure much, so I thought it would be a good exercise
# to just write up an example.

class Node:
    def __init__(self):
        self.left = None
        self.right = None
        self.val = None

class Tree:
    def __init__(self):
        self.root = None

    def parse(self, values):
        nodestack = list() # nodes to be added to the tree
        treestack = list() # breadth first traversal of nodes in the tree
        ptr = None # pointer to the current node in the tree
        # Make a list of nodes to be added to the tree
        for val in values:
            node = Node()
            node.val = val
            nodestack.append(node)
        # If that list is not empty the first one becomes the root
        if nodestack != []:
            self.root = nodestack.pop(0)
            ptr = self.root # We now point to the current root of the tree

        # While we still have nodes to add to the list
        while nodestack != []:
            # If the current node in the tree has no left child
            if ptr.left == None:
                # Set the left child to the node at the top of the nodestack
                ptr.left = nodestack.pop(0)
                # Add this new node in the tree to the list of nodes to be visited
                treestack.append(ptr.left)
            # Otherwise if the current node has no right child
            elif ptr.right == None:
                # Set the right child to the node at the top of the nodestack
                ptr.right = nodestack.pop(0)
                # Add this new node in the tree to the list of nodes to be visited
                treestack.append(ptr.right)
            # If the node has a left and a right child, set the current pointer
            # to the next node in the list of nodes to be visited
            else:
                ptr = treestack.pop(0)

    def bfprint(self):
        # Print the nodes in the tree in a breadth-first traversal
        ptr = None
        treestack = list()
        if self.root != None:
            treestack.append(self.root)
        while treestack != []:
            ptr = treestack.pop(0)
            print(ptr.val, end=' ', flush=True)
            if ptr.left != None:
                treestack.append(ptr.left)
            if ptr.right != None:
                treestack.append(ptr.right)
        print('', flush=True)

    def dfprinter(node):
        print(node.val, end=' ', flush=True)
        if node.left != None:
            Tree.dfprinter(node.left)
        if node.right != None:
            Tree.dfprinter(node.right)

    def dfprint(self):
        if self.root != None:
            Tree.dfprinter(self.root)
        print('', flush=True)

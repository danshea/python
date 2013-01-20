#!/usr/bin/env python
#
# Name: sequencer.py
# Date: 2013-01-19
# Author: Dan Shea
# Description: Based on reads in some file, construct a De Bruijn graph and
# show the Eulerian walks to result in possible sequences of the reconstructed
# reads.  This is my attempt at writing a general purpose, dumb sequencer to
# ensure I understand the practice behind the theory.
#
# Some properties of directed graphs we need to consider when examing the graph
# for the existence of an Eulerian cycle or trail.
#
# A directed graph has an Eulerian cycle if and only if every vertex has equal
# in degree and out degree, and all of its vertices with nonzero degree belong
# to a single strongly connected component. Equivalently, a directed graph has
# an Eulerian cycle if and only if it can be decomposed into edge-disjoint
# directed cycles and all of its vertices with nonzero degree belong to a single
# strongly connected component.
#
# A directed graph has an Eulerian trail if and only if at most one vertex has
# (out-degree) - (in-degree) = 1, at most one vertex has
# (in-degree) - (out-degree) = 1, every other vertex has equal in-degree and
# out-degree, and all of its vertices with nonzero degree belong to a single
# connected component of the underlying undirected graph.

class Vertice(object):
    def __init__(self, key, in_edges=None, out_edges=None):
        '''Vertice(key, in_edges, out_edges) - given a key and a list of
           in_edges and a list of out_edges, create a Vertice object'''
        self.key = key
        if in_edges:
            self.in_edges = in_edges
        else:
            self.in_edges = list()
        if out_edges:
            self.out_edges = out_edges
        else:
            self.out_edges = list()
        self.in_degree = len(self.in_edges)
        self.out_degree = len(self.out_edges)

    def append_out_edge(self, edge):
        self.out_edges.append(edge)
        self.out_degree += 1

    def append_in_edge(self, edge):
        self.in_edges.append(edge)
        self.in_degree += 1

class DeBruijn(object):
    def __init__(self, filename):
        '''DeBruijn(filename) - open filename and construct a DeBruijn graph
           based on the reads in the file.'''
        # initialize variables
        self.kmers = list()
        self.alphabet = set()
        self.vertices = dict()

        # determine the alphabet
        self._alphabet(filename)
        # construct the graph
        self._construct_graph(filename)

    def _alphabet(self, filename):
        '''_alphabet(filename) - determine the alphabet used by examining all
           reads in the file and constructing the set of unique symbols.'''
        with open(filename, 'r') as fh:
            for read in fh:
                # The union f distinct elements in the read and the existing
                # set of known symbols gives us the new alphabet
                self.alphabet = self.alphabet.union(set(read.strip()))

    def _overlap(self, a, b):
        '''_overlap(a, b) - given to n-length sequences determine if they
           overlap and return the largest overlap offset or -1 if no overlap
           exists.'''
        offset = 0
        length = len(a)
        while offset < length:
            csum = sum([ord(i)^ord(j) for i,j in zip(a[offset:],b[0:length-offset])])
            if csum == 0:
                return offset
            offset += 1
        return -1

    def _construct_graph(self, filename):
        '''_construct_graph(filename) - construct the graph representation of
           the k-mer reads in filename'''
        # NOTE: This constructs a directed graph containing all reads as
        # vertices with edges pointing to the overlapping vertices.
        with open(filename, 'r') as fh:
            # load every read into memory
            for read in fh:
                self.kmers.append(read.strip())
        # Compare every read to every other read, determining overlaps
        for a in self.kmers:
            if not self.vertices.has_key(a):
                self.vertices[a] = Vertice(a)
            for b in self.kmers:
                # each vertices edges are stored in a dict with the offset of the
                # overlap
                if not self.vertices.has_key(b):
                    self.vertices[b] = Vertice(b)
                # compute overlap, if any
                overlap = self._overlap(a,b)
                # if there is overlap, add an out_edge to vertice a
                # and an in_edge to vertice b
                if overlap > -1:
                    self.vertices[a].append_out_edge((b,overlap))
                    self.vertices[b].append_in_edge((a,overlap))

    def _is_eulerian(self):
        '''_is_eulerian() - return True if and only if every vertex has equal in
           degree and out degree'''
        for v in self.vertices:
            if self.vertices[v].in_degree != self.vertices[v].out_degree:
                return False
        return True

    def _is_eulerian_trail(self):
        '''_is_eulerian_trail() - return True if and only if at most one vertex
           has (out-degree) - (in-degree) = 1, at most one vertex has
           (in-degree) - (out-degree) = 1, every other vertex has equal
           in-degree and out-degree'''
        out_in_flag = True
        in_out_flag = True
        for v in self.vertices:
            if self.vertices[v].in_degree == self.vertices[v].out_degree:
                continue
            elif self.vertices[v].out_degree - self.vertices[v].in_degree == 1 \
                 and out_in_flag:
                out_in_flag = False
                continue
            elif self.vertices[v].in_degree - self.vertices[v].out_degree == 1 \
                 and in_out_flag:
                in_out_flag = False
                continue
            else:
                return False
        return True

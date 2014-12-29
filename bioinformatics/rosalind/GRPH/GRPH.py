#!/usr/bin/env python3.4
"""
Filename:     GRPH.py
Date Created: 2014-12-29 23:37
Author:       dshea <danshea@iastate.edu>
Description:  The file description goes here.
"""
import argparse
import sys
from Bio import SeqIO

def buildGraph(records, k=3):
    """ Given a list of SeqRecord objects, construct a dictionary detailing
        directed edges from s to t such that s !=t with overlap of k nucleotides """
    # Dictionary to store the graph
    graph = dict()
    # Track the current position in the list
    current = 0
    for s in records:
        # Make a copy of all the edges sans the current edge
        nodes = records[:current] + records[current+1:]
        for t in nodes:
            #print('{}\t{}'.format(s.seq[len(s.seq)-k:], t.seq[:k]))
            if(s.seq[len(s.seq)-k:] == t.seq[:k]):
                if(s.id in graph.keys()):
                    graph[s.id].append(t)
                else:
                    graph[s.id] = [t]
        current = current + 1
    return(graph)

def readFasta(infile):
    """ Reads in record(s) from FASTA file returning a list of SeqRecord objects """
    records = list()
    for seq in SeqIO.parse(infile, "fasta"):
        records.append(seq.upper())
    return(records)

def printGraph(graph):
    """ Given a dictionary representing the directed edges, display the pairs """
    for s in graph.keys():
        for t in graph[s]:
            print('{} {}'.format(s, t.id))

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
    printGraph(buildGraph(readFasta(args.infile.name)))

if __name__ == '__main__':
    main()

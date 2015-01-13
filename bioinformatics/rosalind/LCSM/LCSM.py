#!/usr/bin/env python3.4
"""
Filename:     LCSM.py
Date Created: 2015-01-13 14:41
Author:       dshea <danshea@iastate.edu>
Description:  Given:  A collection of k (k≤100) DNA strings of length at most
                      1 kbp each in FASTA format.
              Return: A longest common substring of the collection. (If multiple
                      solutions exist, you may return any single solution.)
"""
import argparse
import sys
from Bio import SeqIO

def readFasta(infile):
    """ Reads in record(s) from FASTA file returning a list of SeqRecord objects """
    records = list()
    for seq in SeqIO.parse(infile, "fasta"):
        records.append(seq.upper())
    return(records)

def findLongestCommonSubstring(records):
    """ Return the longest common substring from a list of strings """
    # The shortest string could be the longest possible common substring,
    # therefore, compare it to the others
    strings = records[:]
    shortest = list(map(len,strings))
    shortest = strings.pop(shortest.index(min(shortest)))
    possible_subs = set([shortest[i:j] for i in range(0,len(shortest)) for j in range(i+1,len(shortest)+1)])
    substrings = list()
    for s in strings:
        substrings.append(set())
        for possible_sub in possible_subs:
            if possible_sub in s:
                substrings[-1].add(possible_sub)
    solution = possible_subs
    for sols in substrings:
        solution = set.intersection(solution, sols)
    solution = list(solution)
    sol_lens = list(map(len,solution))
    return(solution[sol_lens.index(max(sol_lens))])

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

    records = [str(seqrecord.seq) for seqrecord in readFasta(args.infile.name)]
    print(findLongestCommonSubstring(records))

if __name__ == '__main__':
    main()

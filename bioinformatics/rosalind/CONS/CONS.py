#!/usr/bin/env python3.4
"""
Filename:     CONS.py
Date Created: 2014-12-29 23:37
Author:       dshea <danshea@iastate.edu>
Description:  Given: A collection of at most 10 DNA strings of equal length
              (at most 1 kbp) in FASTA format.
              Return: A consensus string and profile matrix for the collection.
              (If several possible consensus strings exist, then you may return
              any one of them.)
"""
import argparse
import sys
from Bio import SeqIO
import pandas

def buildMatrix(records):
    """ Construct the profile matrix from a list of SeqRecords
        Return the resulting data frame object
    """
    # Create a data frame
    df = pandas.DataFrame(data=records)
    # transpose it so columns are compriosed of sequences
    df = df.transpose()
    return(df)

def countNucleotides(sequence, nt):
    """ Count the total number of a given nucleotide, returning the count """
    total = 0
    for base in sequence:
        if base == nt:
            total = total + 1
    return(total)

def buildProfile(df):
    """ Construct a data frame that consists of the profile matrix for the given
        sequences
    """
    A = df.apply(countNucleotides, axis=1, args=('A'))
    C = df.apply(countNucleotides, axis=1, args=('C'))
    G = df.apply(countNucleotides, axis=1, args=('G'))
    T = df.apply(countNucleotides, axis=1, args=('T'))
    return(pandas.DataFrame(data=[A,C,G,T], index=['A','C','G','T']))

def readFasta(infile):
    """ Reads in record(s) from FASTA file returning a list of SeqRecord objects """
    records = list()
    for seq in SeqIO.parse(infile, "fasta"):
        records.append(seq.upper())
    return(records)

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
    profile = buildProfile(buildMatrix(readFasta(args.infile.name)))
    # Display the consensus string
    print('{}'.format(''.join(profile.idxmax())))
    # Display the profile matrix
    print('A: {}'.format(' '.join([str(x) for x in profile.transpose()['A']])))
    print('C: {}'.format(' '.join([str(x) for x in profile.transpose()['C']])))
    print('G: {}'.format(' '.join([str(x) for x in profile.transpose()['G']])))
    print('T: {}'.format(' '.join([str(x) for x in profile.transpose()['T']])))

if __name__ == '__main__':
    main()

#!/usr/bin/env python3.4
"""
Filename:     SPLC.py
Date Created: 2014-12-30 20:14
Author:       dshea <danshea@iastate.edu>
Description:  Given: A DNA string s (of length at most 1 kbp) and a collection
              of substrings of s acting as introns. All strings are given in
              FASTA format.
              Return: A protein string resulting from transcribing and
              translating the exons of s. (Note: Only one solution will exist
              for the dataset provided.)
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

def removeIntrons(records):
    """ Remove the introns from the DNA string
        The first SeqRecord is the DNA sequence
        The other SeqRecords in the list of are the defined introns
    """
    # Copy the records to a new list, we're going to destructively modify it
    introns = records[:]
    # First SeqRecord is the DNA sequence
    DNA = introns.pop(0)
    for intron in introns:
        # Remove each intron from the DNA sequence
        # Where the intron starts
        start = str(DNA.seq).index(str(intron.seq))
        DNA = DNA[:start] + DNA[start+len(intron):]
    # Return the coding sequence (DNA with introns removed)
    return(DNA)

def translateSeq(seq):
    """ Take the SeqRecord and return the protein sequence as a string """
    # Drop the '*' that BioPython uses to denote the end of the sequence
    # and return the translation
    return(str(seq.seq.translate())[:len(seq.seq.translate())-1])

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
    print(translateSeq(removeIntrons(readFasta(args.infile.name))))

if __name__ == '__main__':
    main()

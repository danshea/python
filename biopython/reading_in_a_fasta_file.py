#!/usr/bin/env python
from Bio import SeqIO
import sys

'''
Using BioPython to read in a fasta file.
'''

def readFasta(filename):
    '''readFasta - read in fasta formatted file filename and print each record id'''
    with open(filename, 'r') as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            print record.id

def readEntireFasta(filename):
    '''readEntireFasta - read in the entire fasta formatted file filename and return a list of the records'''
    try:
        with open(filename, 'r') as fh:
            records = list(SeqIO.parse(fh, 'fasta'))
            return(records)
    except IOError as err:
        sys.stderr.write(repr(err))
        return(None)
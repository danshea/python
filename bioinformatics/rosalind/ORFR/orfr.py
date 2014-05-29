#!/usr/bin/env python

################################################################################
#
# Name: orfr.py
# Date: 2014-04-26
# Author: dshea
# Description:
#
# Problem
#
# An ORF begins with a start codon and ends either at a stop codon or at the end
# of the string. We will assume the standard genetic code for translating an RNA
# string into a protein string (i.e., see the standard RNA codon table).
#
# ORF finder from the SMS 2 package can be run online here.
#
# Given: A DNA string s of length at most 1 kbp.
#
# Return: The longest protein string that can be translated from an ORF of s. If
# more than one protein string of maximum length exists, then you may output any
# solution.
#
# Sample Dataset
#
# AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG
#
# Sample Output
#
# MLLGSFRLIPKETLIQVAGSSPCNLS
#
# Programming Shortcut
#
# We can also find ORFs using the EMBOSS program getorf. It can be downloaded and
# run locally. The documentation can be found here. To find ORFs using Biopython,
# it may be useful to recall the translate() and reverse_complement() methods from
# the Bio.Seq module.
#
################################################################################

import sys
from Bio import Seq
import regex

def pad_seq(seq):
    if len(str(seq)) % 3 != 0:
        return (seq + 'N'*(3 - (len(str(seq)) % 3)))
    else:
        return(seq)

def find_ORF(dna_sequence):
    # Create the sequence
    seq = Seq.Seq(dna_sequence)
    # Translate all 6 reading frames
    trans = [pad_seq(seq).translate(),
             pad_seq(seq[1:]).translate(),
             pad_seq(seq[2:]).translate(),
             pad_seq(seq.reverse_complement()).translate(),
             pad_seq(seq.reverse_complement()[1:]).translate(),
             pad_seq(seq.reverse_complement()[2:]).translate()]
    orfs = list()
    for t in trans:
        m = regex.findall('M\w*\*', str(t))
        if m != []:
            m = max(regex.findall('M\w*\*', str(t)), key=lambda a: len(a))
            orfs.append(m)
    return(max(orfs, key=lambda a: len(a)))
    #return(orfs)

def main():
    if len(sys.argv) != 2:
        print 'usage: {0:s} infile'.format(sys.argv[0])
        sys.exit(1)
    else:
        infile = sys.argv[1]
        with open(infile, 'r') as fh:
            print find_ORF(fh.readline().strip())
        sys.exit(0)

if __name__ == '__main__':
    main()


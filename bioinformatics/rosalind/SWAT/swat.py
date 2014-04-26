#!/usr/bin/env python

################################################################################
#
# Name: swat.py
# Date: 2014-04-26
# Author: dshea
# Description:
# Problem
#
# The EBI portal for Water, a program for local alignment in the EMBOSS suite, can
# be accessed here.
#
# Use:
#
#    The BLOSUM62 scoring matrix. Gap opening penalty of 10. Gap extension
#    penalty of 1.
#
# Given: Two UniProt ID's corresponding to two protein strings s and t.
#
# Return: The maximum score of any local alignment of s and t. Sample Dataset
#
# B3ET80 Q78PG9
#
# Sample Output
#
# 35
#
################################################################################

import os.path
import sys
from Bio.Emboss.Applications import WaterCommandline
from Bio import Entrez
from Bio import SeqIO

# GLOBAL
Entrez.email = "shea.d@husky.neu.edu"

def getFasta(uniprot_id):
    filename = '{0:s}.fa'.format(uniprot_id)
    if not os.path.isfile(filename):
        handle = Entrez.efetch(db='nucleotide', rettype='fasta', retmode='text', id=uniprot_id)
        seq_record = SeqIO.read(handle, 'fasta')
        handle.close()
        SeqIO.write(seq_record, filename, format='fasta')
    return(filename)

def main():
    if len(sys.argv) != 3:
        print 'usage {0:s} uniprot_id1 uniprot_id2'.format(sys.argv[0])
        sys.exit(1)
    uniprot_a = sys.argv[1]
    uniprot_b = sys.argv[2]
    sequence_a = getFasta(uniprot_a)
    sequence_b = getFasta(uniprot_b)
    water_cline = WaterCommandline('/usr/local/bin/water',
                                     asequence=sequence_a, bsequence=sequence_b,
                                     gapopen=10, gapextend=1, outfile='{0:s}_{1:s}_water.txt'.format(uniprot_a,uniprot_b))
    stdout, stderr = water_cline()
    sys.exit(0)

if __name__ == '__main__':
    main()


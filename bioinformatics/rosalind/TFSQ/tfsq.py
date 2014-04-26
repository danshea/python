#!/usr/bin/env python

################################################################################
#
# Name: tfsq.py
# Date: 2014-04-25
# Author: dshea
# Description:
# Problem
#
# Sometimes it's necessary to convert data from FASTQ format to FASTA format. For
# example, you may want to perform a BLAST search using reads in FASTQ format
# obtained from your brand new Illumina Genome Analyzer.
#
# Links:
#
#     A FASTQ to FASTA converter can be accessed from the Sequence conversion
# website
#
#     A free GUI converter developed by BlastStation is available here for
# download or as an add-on to Google Chrome.
#
#     There is a FASTQ to FASTA converter in the Galaxy web platform. Note that
# you should register in the Galaxy and upload your file prior to using this
# tool.
#
# Given: FASTQ file
#
# Return: Corresponding FASTA records Sample Dataset
#
# @SEQ_ID GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT +
# !*((((***+))%%%++)(%%%%).1***-+*****))**55CCF>>>>>>CCCCCCC65
#
# Sample Output
#
# >SEQ_ID
# GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
#
################################################################################

import os.path
import sys
from Bio import SeqIO

def main():
    if len(sys.argv) != 2:
        print 'usage: {0:s} fastq_file'.format(sys.argv[0])
        sys.exit(1)
    infile = sys.argv[1]
    outfile = infile.split(os.path.extsep)
    outfile = os.path.extsep.join(outfile[0:len(outfile)-1]+['fa'])
    SeqIO.write(SeqIO.parse(infile, format='fastq'),outfile, format='fasta')
    sys.exit(0)

if __name__ == '__main__':
    main()


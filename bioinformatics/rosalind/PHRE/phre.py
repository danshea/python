#!/usr/bin/env python

################################################################################
#
# Name: phre.py
# Date: 2014-04-26
# Author: dshea
# Description:
# Problem
#
#     A version of FastQC can be downloaded here and run locally on any operating
# system with a suitable Java Runtime Environment (JRE) installed.
#
#     An online version of FastQC is also available here in the "Andromeda"
# Galaxy instance.
#
# Given: A quality threshold, along with FASTQ entries for multiple reads.
#
# Return: The number of reads whose average quality is below the threshold.
# Sample Dataset
#
# 28
# @Rosalind_0041
# GGCCGGTCTATTTACGTTCTCACCCGACGTGACGTACGGTCC
# +
# 6.3536354;.151<211/0?::6/-2051)-*"40/.,+%)
# @Rosalind_0041
# TCGTATGCGTAGCACTTGGTACAGGAAGTGAACATCCAGGAT
# +
# AH@FGGGJ<GB<<9:GD=D@GG9=?A@DC=;:?>839/4856
# @Rosalind_0041
# ATTCGGTAATTGGCGTGAATCTGTTCTGACTGATAGAGACAA
# +
# @DJEJEA?JHJ@8?F?IA3=;8@C95=;=?;>D/:;74792.
#
# Sample Output
#
# 1
#
#     Programming Shortcut
#
#     When used to read FASTQ data, BioPython's function SeqIO.parse returns a
# SeqRecord object containing the Phred quality scores corresponding to each base
# of the sequence. The scores are found in the .letter_annotations attribute,
# which is a Python dictionary having, in this case, a single key,
# 'phred_quality.'
#
#     >>> print record.letter_annotations.keys()
#     ['phred_quality']
#     >>> print record.letter_annotations["phred_quality"]
#     [40, 39, 38, 37, 36, 35, 34, 13, 12, 11, 10, 9, 8]
#
#
################################################################################

import sys
from Bio import SeqIO

def average_quality(seq_record):
    q_scores = seq_record.letter_annotations['phred_quality']
    num_q_scores = len(q_scores)
    return(sum(q_scores)/float(num_q_scores))

def report_q_scores(infile, threshold):
    total = 0
    seqio = SeqIO.parse(infile, format='fastq')
    for seq_record in seqio:
        avg_q_score = average_quality(seq_record)
        if avg_q_score < threshold:
            total += 1
    return(total)

def main():
    if len(sys.argv) != 3:
        print 'usage: {0:s} threshold infile'.format(sys.argv[0])
        sys.exit(1)
    else:
        threshold = float(sys.argv[1])
        infile = sys.argv[2]
        print report_q_scores(infile, threshold)
    sys.exit(0)

if __name__ == '__main__':
    main()


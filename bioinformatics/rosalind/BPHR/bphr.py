#!/usr/bin/env python

################################################################################
#
# Name: bphr.py
# Date: 2014-04-26
# Author: dshea
# Description:
# Problem
#
# Quality of the bases can vary depends on position in read due to nature of the
# sequencing procedure. One can check this quality distribution using "Per Base
# Sequence Quality" module of the FastQC program.
#
# Average accepted quality values is a 10 for the lower quartile and 25 for
# median. If the values falls below this limit the module rises a warning.
#
# Note that for the reads >50bp long FastQC will group the bases. To show data
# for every base in the read use "--nogroup" option.
#
# Given: FASTQ file, quality threshold q
#
# Return: Number of positions where mean base quality falls below given threshold
# Sample Dataset
#
# 26
# @Rosalind_0029
# GCCCCAGGGAACCCTCCGACCGAGGATCGT
# +
# >?F?@6<C<HF?<85486B;85:8488/2/
# @Rosalind_0029
# TGTGATGGCTCTCTGAATGGTTCAGGCAGT
# +
# @J@H@>B9:B;<D==:<;:,<::?463-,,
# @Rosalind_0029
# CACTCTTACTCCCTAGCCGAACTCCTTTTT
# +
# =88;99637@5,4664-65)/?4-2+)$)$
# @Rosalind_0029
# GATTATGATATCAGTTGGCTCCGAGAGCGT
# +
# <@BGE@8C9=B9:B<>>>7?B>7:02+33.
#
# Sample Output
#
# 17
#
################################################################################

import sys
from Bio import SeqIO
import pandas

def report_q_scores(infile, threshold):
    total = 0
    seqio = SeqIO.parse(infile, format='fastq')
    q_scores = list()
    for seq_record in seqio:
        q_scores.append(seq_record.letter_annotations['phred_quality'])
    df = pandas.DataFrame(q_scores)
    means = df.mean()
    for m in means:
        if m < threshold:
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
#!/bin/bash

################################################################################
#
# Name: filt.sh
# Date: 2014-04-26
# Author: Dan Shea <shea.d@husky.neu.edu>
# Description:
# Problem
#
# Poor-quality reads can be filtered out using the FASTQ Quality Filter tool from
# the FASTX toolkit. A command-line version of FASTX can be downloaded for Linux
# or MacOS from its website. An online interface for the FASTQ Quality Filter is
# also available here within the Galaxy web platform.
#
# Given: A quality threshold value q, percentage of bases p, and set of FASTQ
# entries.
#
# Return: Number of reads in filtered FASTQ entries
# Sample Dataset
#
# 20 90
# @Rosalind_0049_1
# GCAGAGACCAGTAGATGTGTTTGCGGACGGTCGGGCTCCATGTGACACAG
# +
# FD@@;C<AI?4BA:=>C<G=:AE=><A??>764A8B797@A:58:527+,
# @Rosalind_0049_2
# AATGGGGGGGGGAGACAAAATACGGCTAAGGCAGGGGTCCTTGATGTCAT
# +
# 1<<65:793967<4:92568-34:.>1;2752)24')*15;1,.3*3+*!
# @Rosalind_0049_3
# ACCCCATACGGCGAGCGTCAGCATCTGATATCCTCTTTCAATCCTAGCTA
# +
# B:EI>JDB5=>DA?E6B@@CA?C;=;@@C:6D:3=@49;@87;::;;?8+
#
# Sample Output
#
# 2
#
################################################################################

if [[ $# != 3 ]]; then
    echo "usage: ${0} q_threshold p_bases infile"
    exit 1
else
    q_threshold=$1
    p_bases=$2
    infile=$3
    outfile=$(echo ${infile} | cut -d. -f1)
    outfile="${outfile}_filt.fq"
    fastq_quality_filter -q ${q_threshold} -p ${p_bases} -i ${infile} -o ${outfile}
    grep -c -e '^@Rosalind_' ${outfile}
fi
exit 0

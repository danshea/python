#!/usr/bin/env python3.4
"""
Filename:     ORF.py
Date Created: 2015-01-13 06:39
Author:       dshea <danshea@iastate.edu>
Description:  Given: A DNA string s of length at most 1 kbp in FASTA format.
              Return: Every distinct candidate protein string that can be
              translated from ORFs of s. Strings can be returned in any order.
"""
import argparse
import sys
from Bio import SeqIO, Seq
import re

def readFasta(infile):
    """ Reads in record(s) from FASTA file returning a list of SeqRecord objects """
    records = list()
    for seq in SeqIO.parse(infile, "fasta"):
        records.append(seq.upper())
    return(records)

def translateSeq(seqs):
    """ Take the SeqRecord and return a list of possible translations to protein
        for that sequence in all 6 reading frames """
    start_codons = ['ATG']
    stop_codons  = ['TAG', 'TGA', 'TAA']
    for seq in seqs:
        seq = seq.seq
        # split the sequence into codons for the 3 fwd and 3 rev reading frames
        rfs = list()
        for j in range(3):
            tmp = [str(seq)[i:i+3] for i in range(j, len(str(seq)), 3)]
            rtmp = [str(seq.reverse_complement())[i:i+3] for i in range(j, len(str(seq)), 3)]
            rfs.append(tmp)
            rfs.append(rtmp)

        orfs = dict()
        rfindex = 1
        for rf in rfs:
            index = 0
            orfs[rfindex] = dict()
            for codon in rf:
                if(codon in start_codons):
                    if 'start' in orfs[rfindex].keys():
                        orfs[rfindex]['start'].append(index)
                    else:
                        orfs[rfindex]['start'] = [index]
                if(codon in stop_codons):
                    if 'stop' in orfs[rfindex].keys():
                        orfs[rfindex]['stop'].append(index)
                    else:
                        orfs[rfindex]['stop'] = [index]
                index += 1
            rfindex += 1
        solution = set()
        for rfindex in orfs.keys():
            if 'start' in orfs[rfindex].keys():
                for sindex in sorted(orfs[rfindex]['start']):
                    for stindex in sorted(orfs[rfindex]['stop']):
                        if stindex > sindex:
                            #print('(frame {}: {},{})'.format(rfindex, sindex, stindex))
                            #print(Seq.Seq(''.join(rfs[rfindex-1][sindex:stindex+1])).translate()[:-1])
                            solution.add(str(Seq.Seq(''.join(rfs[rfindex-1][sindex:stindex+1])).translate()[:-1]))
                            break
        for s in solution:
            print(s)

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
    translateSeq(readFasta(args.infile.name))

if __name__ == '__main__':
    main()

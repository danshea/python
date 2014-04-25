#!/usr/bin/env python

################################################################################
#
# Name: data_formats.py
# Date: 2014-04-24
# Author: dshea
# Description:
#
# GenBank can be accessed here. A detailed description of the GenBank format can
# be found here. A tool, from the SMS 2 package, for converting GenBank to FASTA
# can be found here.
#
# Given: A collection of n (n<=10) GenBank entry IDs.
#
# Return: The shortest of the strings associated with the IDs in FASTA format.
#
################################################################################

from Bio import Entrez
from Bio import SeqIO
Entrez.email = 'shea.d@husky.neu.edu'

def shortest(genbank_ids):
    handle = Entrez.efetch(db='nucleotide', id=genbank_ids, rettype="fasta")
    records = list(SeqIO.parse(handle, "fasta"))
    seqlens = [len(record.seq) for record in records]
    print '>'+records[seqlens.index(min(seqlens))].description
    print records[seqlens.index(min(seqlens))].seq

if __name__ == '__main__':
    shortest(['GU292427', 'JX472277', 'NM_001025158' ,'JX308819', 'JX317645','JQ796071', 'JF927163' ,'NM_001270868'])

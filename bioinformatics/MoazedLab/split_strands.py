#!/usr/bin/env python
'''
Name: split_strands.py
Date: 2014-08-13
Author: Dan Shea <daniel.john.shea@gmail.com>
Description:

    We would like to examine RNA-Seq information generated from tophat in IGV
    using the strand information that is present in the BAM files tophat
    creates.  IGV seems to not use the flag value for reverse complemented,
    which we may use to discern which strand the read came from.

    Fortunately, pysam can do this for us relatively simply.  Otherwise we would
    need to take the flag value and bitwise AND it with 2**5 to read the value
    of the 5th bit in the flag.
'''

import argparse
import pysam
import os
import sys

def pe_strandify(inputfile, foutputfile, routputfile):
    try:
        sfile = pysam.Samfile(inputfile,'rb')
        try:
            fout = pysam.Samfile(foutputfile, 'wb', template=sfile)
            rout = pysam.Samfile(routputfile, 'wb', template=sfile)
            for read in sfile:
                if not read.is_unmapped:
                    if read.is_paired and read.is_proper_pair and read.is_read1:
                        pos = sfile.tell()
                        #print pos
                        try:
                            #print 'calling mate()'
                            mate = sfile.mate(read)
                            if read.is_reverse:
                                rout.write(read)
                                rout.write(mate)
                            else:
                                fout.write(read)
                                fout.write(mate)
                        except ValueError:
                            continue
                        finally:
                            #print 'calling seek()'
                            sfile.seek(pos)
            sfile.close()
            fout.close()
            rout.close()
        except IOError:
            print 'Failed to create output files {} and {} for writing.'.format(foutputfile, routputfile)
            sys.exit(2)
    except IOError:
        print 'Failed to open {} for reading.'.format(inputfile)
        sys.exit(3)
    except ValueError:
        print 'Is {} in BAM format?'.format(inputfile)
        sys.exit(4)

def strandify(inputfile, foutputfile, routputfile):
    try:
        sfile = pysam.Samfile(inputfile,'rb')
        try:
            fout = pysam.Samfile(foutputfile, 'wb', template=sfile)
            rout = pysam.Samfile(routputfile, 'wb', template=sfile)
            for read in sfile:
                if not read.is_unmapped:
                    if read.is_reverse:
                        rout.write(read)
                    else:
                        fout.write(read)
            sfile.close()
            fout.close()
            rout.close()
        except IOError:
            print 'Failed to create output files {} and {} for writing.'.format(foutputfile, routputfile)
            sys.exit(2)
    except IOError:
        print 'Failed to open {} for reading.'.format(inputfile)
        sys.exit(3)
    except ValueError:
        print 'Is {} in BAM format?'.format(inputfile)
        sys.exit(4)

def main():
    parser = argparse.ArgumentParser(description='Split tophat output into two files based on the origin strand (sense or anti-sense)')
    parser.add_argument('inputfile', type=argparse.FileType('r'), help='input file in BAM format')
    parser.add_argument('--paired', help='flag for paired end data', action='store_true')
    args = parser.parse_args()
    inputfile = args.inputfile
    ifpath = os.path.dirname(inputfile.name)
    ifname = os.path.basename(inputfile.name)
    inputfile.close() # Close the file handle as we don't need it open here anymore.
    
    # Before we head off the try and parse the files, let's ensure there is a valid index so we don't segfault
    indexfile = os.path.join(ifpath,inputfile.name+'.bai')
    if not os.path.isfile(indexfile):
        print 'You need a BAM index file named {} to process {}'.format(indexfile, inputfile.name)
        sys.exit(1)

    if args.paired:
        foutputfile = os.path.join(ifpath,'Forward_pe_'+ifname)
        routputfile = os.path.join(ifpath,'Reverse_pe_'+ifname)
        print 'Processing as paired end reads.'
        print 'inputfile is {}'.format(inputfile.name)
        print 'forward strand hits will be in {}'.format(foutputfile)
        print 'reverse strand hits will be in {}'.format(routputfile)
        pe_strandify(inputfile.name, foutputfile, routputfile)
    else:
        foutputfile = os.path.join(ifpath,'Forward_'+ifname)
        routputfile = os.path.join(ifpath,'Reverse_'+ifname)
        print 'Processing as single end reads.'
        print 'inputfile is {}'.format(inputfile.name)
        print 'forward strand hits will be in {}'.format(foutputfile)
        print 'reverse strand hits will be in {}'.format(routputfile)
        strandify(inputfile.name, foutputfile, routputfile)
    sys.exit(0)

if __name__ == '__main__':
    main()
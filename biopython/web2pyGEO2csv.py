#!/usr/bin/env python
'''
Name: web2pyGEO2csv.py
Date: 2014-03-08
Author: Dan Shea
Description:
    BioPython provides a very basic parser for Gene Omnibus Expression files.
    I have attempted to extend some additional functionality to it in order to
    read in GEO files and dump a csv file for use with web2py.
    
    The GEO soft file format is as follows:
    Lines that begin with ^ identify the entity_type
    Valid types are:
    DATABASE
    PLATFORM
    SERIES
    SAMPLE
    
    entity_attribute lines begin with !
    
    For SAMPLE entity_types, we want to construct a TABLE for each and insert
    the table_rows.  col_defs.keys() or table_rows[0] will give is the column
    names for the table.
'''

from Bio import Geo
import sys
import gzip
import os
import os.path
import re
import getopt
import csv

def parseGeo(filename):
    '''parseGeo(filename): open a GEO soft file and return the records it contains.'''
    try:
            fh = gzip.open(filename)
            records = Geo.parse(fh)
            return(records)
    except:
        try:
                fh = open(filename)
                records = Geo.parse(fh)
                return(records)
        except:
            print 'Could not open filename'
            sys.exit(3)

def createCsv(outfile, records):
    '''createCsv(records): create csv file based on the records contained in the GEO file'''
    with open(outfile, 'wb') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',', quotechar="|", quoting=csv.QUOTE_MINIMAL)
        series = ''
        for record in records:
            if re.match('DATABASE', record.entity_type, re.IGNORECASE) != None:
                pass
            elif re.match('PLATFORM', record.entity_type, re.IGNORECASE) != None:
                pass
            elif re.match('SERIES', record.entity_type, re.IGNORECASE) != None:
                series = record.entity_id
            elif re.match('SAMPLE', record.entity_type, re.IGNORECASE) != None:            
                #csvwriter.writerow(record.col_defs.keys())
                for row in record.table_rows[1:]:
                    csvwriter.writerow([series, record.entity_id] + row)
    return 0

def usage(argv):
    print 'USAGE: {0:s}'.format(argv[0])
    print 'Options:'
    print '-h, --help\t\tDisplay this message.'
    print '-i, --infile\t\tThe input GEO soft file (can be gzip compressed.)'
    print '-o, --outfile\t\tThe output csv file.'
    return 1

def main(argv=None):
    infile = None
    outfile = None
    if argv is None:
        argv = sys.argv
    try:
        optlist, args = getopt.getopt(argv[1:], 'hi:o:', ['help','infile=','outfile='])
    except:
        return(usage(argv))
    for o, a in optlist:
        if o in ('-h', '--help'):
            usage(argv)
        elif o in ('-i', '--infile'):
            infile = a
            if not os.path.isfile(a):
                print '{0:s} is not a file.'.format(a)
                return 3
        elif o in ('-o', '--outfile'):
            if os.path.isfile(a):
                print '{0:s} exists.'.format(a)
                print 'Exiting without modifying file.'
                return 4
            else:
                outfile = a
    if infile is None or outfile is None:
        usage(argv)        
    else:
        records = parseGeo(infile)
        return(createCsv(outfile, records))

if __name__ == '__main__':
    sys.exit(main())
#!/usr/bin/env python

from Bio import Geo
import sys
import gzip
import os
import os.path

def parseGeo(filename):
    try:
        with gzip.open(filename) as fh:
            records = Geo.parse(fh)
            for record in records:
                print record
    except:
        try:
            with open(filename) as fh:
                records = Geo.parse(fh)
                for record in records:
                    print record
        except:
            print 'Could not open filename'
            sys.exit(3)

def main(argv=None):
    if argv is None:
        argv = sys.argv
    if len(argv) != 2:
        print 'Please provide a filename.'
        return 2
    elif not os.path.isfile(argv[1]):
        print '{0:s} is not a file.'.format(argv[1])
        return 3
    else:
        parseGeo(argv[1])
        return 0

if __name__ == '__main__':
    sys.exit(main())
#!/usr/bin/env python

import getopt
import sys
import os

def usage():
    print 'usage: {0} -a -b value [args]'.format(os.path.basename(sys.argv[0]))
    sys.exit(1)

def main(*args):
    try:
        optlist, args = getopt.getopt(args, 'ab:')
        for opt, val in optlist:
            print 'Option {0} has value {1}'.format(opt, val)
        for arg in args:
            print 'Argument {0}'.format(arg)
    except getopt.GetoptError:
        usage()

if __name__ == '__main__':
    main(*sys.argv[1:])
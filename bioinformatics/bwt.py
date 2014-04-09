#!/usr/bin/env python

################################################################################
#
# Name: bwt.py
# Date: 2014-04-07
# Author: dshea
# Description: Burrows-Wheeler Transformation on an input string
#
################################################################################

from collections import OrderedDict
import sys

class bwtString(object):

    def build_BWTmatrix(self, input_string):
        '''build_BWTmatrix(input_string): construct a matrix of permutations for bwt'''
        input_string = input_string + '\0'
        od = OrderedDict()
        bwt = list()
        od[input_string] = 0
        bwt.append(input_string)
        for i in range(1,len(input_string)):
            od[input_string[i:]+input_string[0:i]] = i
            bwt.append(input_string[i:]+input_string[0:i])
        bwt.sort()
        index = [od[''.join(s)] for s in bwt]
        string = ''.join([input_string[-1] for input_string in bwt])
        return(string,index)

    def __init__(self, input_string):
        self.string, self.index = self.build_BWTmatrix(input_string)
    
    def decode(self):
        '''decode(): return the original input string'''
        orig = range(len(self.string))
        for i in range(len(self.string)):
            orig[self.index[i]-1] = self.string[i]
        return(''.join(orig)[0:-1])

def main():
    return(0)

if __name__ == '__main__':
    sys.exit(main())


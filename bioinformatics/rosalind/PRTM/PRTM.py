#!/usr/bin/env python3.4
"""
Filename:     PRTM.py
Date Created: 2014-12-30 18:24
Author:       dshea <danshea@iastate.edu>
Description:  Given: A protein string P of length at most 1000 aa.
              Return: The total weight of P. Consult the monoisotopic mass
              table.
"""
import argparse
import sys

MONOISOTOPIC_MASS = {
'A':   71.03711,
'C':   103.00919,
'D':   115.02694,
'E':   129.04259,
'F':   147.06841,
'G':   57.02146,
'H':   137.05891,
'I':   113.08406,
'K':   128.09496,
'L':   113.08406,
'M':   131.04049,
'N':   114.04293,
'P':   97.05276,
'Q':   128.05858,
'R':   156.10111,
'S':   87.03203,
'T':   101.04768,
'V':   99.06841,
'W':   186.07931,
'Y':   163.06333,
}

def main():
    """ The main() """
    parser = argparse.ArgumentParser(description='Description goes here.')
    parser.add_argument("sequence", type=str)
    args = parser.parse_args()
    weight = 0.0
    for aa in args.sequence.upper():
        weight += MONOISOTOPIC_MASS[aa]
    print('{0:0.3f}'.format(round(weight,3)))

if __name__ == '__main__':
    main()

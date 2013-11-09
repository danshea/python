#!/usr/bin/env python
'''
Author: djs
Date: 2012-07-07
Description: IEEE Floating Point (IEEE-754) for 32 bit numbers
(Single Precision Floats)

32 bit floating point number is broken down as follows

1-bit sign bit: 0 is positive number, 1 is a negative number

8-bit exponent: 2**8 possible values, there is a bias of 2**(8-1) - 1
that effectively partitions the value range.  i.e. - e = 1 is stored as 1 + 127
So 128 is stored in the 8-bit field '1000 0000'

23-bit unsigned mantissa: holds the mantissa m where m is the unsigned binary
value used to construct the mantissa 1.m
i.e. - 1.5 is stored as '100 0000 0000 0000 0000 0000'
Starting from the left values go from 2**-1 ... 2**-23

Valid 32-bit float value range for normalized are as follows:
+/- 2**-126 <= x <= (2 - 2**-23) * 2*127

Special Values:
(Note: mantissa is also referred to as fraction)
Zero: exponent is 0, mantissa is 0
Denormalized: exponent is all 0, mantissa is non-zero
Infinity: exponent of all 1's, fraction of all 0's
Not a Number (NaN): exponent of all 1's and a non-zero fraction

***********************************************************************
* WARNING: This code does not catch edge cases, underflow or overflow *
*          Use at your own risk!!!!!                                  *
***********************************************************************
'''

def exponent32(N):
    '''
    exponent32 returns the value of exponent without the bias.
    '''
    N = abs(N)
    e = 128
    while 2**e > N:
        e -= 1
    return e

def mantissa32(N):
    '''
    mantissa32 returns the string representing the binary value stored in the
    fraction portion of the floating point spec.
    '''
    N = abs(N)
    e = exponent32(N)
    remainder = N / 2.0**e
    remainder -= 1
    mantissa = [0 for i in xrange(23)]
    for b in xrange(1, 24):
        if 2**-b <= remainder:
            mantissa[b-1] = 1
            remainder = remainder - 2**-b
        else:
            mantissa[b-1] = 0
    return ''.join([str(i) for i in mantissa])

def float32(N):
    '''
    float32 returns a tuple of strings (sign, exponent, mantissa) representing
    the binary values to be stored in a 32-bit single precision float.
    '''
    bias = 127 # 2**(8-1) - 1
    if N > 0:
        sign = '0'
    else:
        sign = '1'
    exponent = bin(exponent32(N) + bias)[2:].zfill(8)
    mantissa = mantissa32(N)
    return (sign, exponent, mantissa)

def num32(sign, exponent, mantissa):
    '''
    num32 computes the floating point value of a tuple of strings
    (sign, exponent, mantissa) that represents the binary values of a 32-bit
    single precision float.
    '''
    bias = 127 # 2**(8-1) - 1
    mvalue = 1.0 + sum([int(i)*2**-(e+1) for e, i in enumerate(list(mantissa))])
    if int(sign):
        return -mvalue * 2**(int(exponent, 2) - bias)
    else:
        return  mvalue * 2**(int(exponent, 2) - bias)
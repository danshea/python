#!/usr/bin/env python3

'''
Author: djs
Date: 2012-07-08
Description: Python3 implementation of IEEE floating point

For a single precision floating point number we have the following:

1:8:23 sign:exponent:fraction

For a double precision floating point number we have the following:

1:11:52 sign:exponent:fraction
'''

def float64(num):
    '''
    float64 return a tuple (sign, exponent, fraction) that represents the ieee
    double precision floating point specification.
    '''
    if type(num) is int or type(num) is float:
        num = float(num)
    else:
        raise ValueError('argument must be float or int')

    # split the number into its integer and fraction components
    integer, fraction = str(num).split('.')

    # determine the sign bit
    integer = int(integer)

    if integer < 0:
        sign = '1'
    else:
        sign = '0'


    # compute the exponent by a series of bit shifts right to determine the
    # exponents value

    integer = abs(integer)

    exponent = 0
    for e in range(0, 1026):
        if integer >> e == 1:
            exponent = e
            break

    # compute the mantissa (fraction) first for the integer portion
    int_mantissa = integer / 2**exponent

    # determine the offset (if any) for the fraction and then compute the
    # mantissa (fraction) for the fraction portion
    offset = 1
    f = list(fraction)
    while f[0] == '0':
        f.pop(0)
        offset += 1
    fraction = int(fraction)
    f_mantissa = fraction / 2**exponent * 10**-offset

    mantissa = int_mantissa + f_mantissa - 1
    # print('mantissa = {0:0.52f}'.format(mantissa))

    # convert the mantissa (fraction) to binary
    b_mantissa = ''
    for i in range(1,53):
        if mantissa >= 2**-i:
            mantissa -= 2**-i
            b_mantissa += '1'
        else:
            b_mantissa += '0'

    # apply the bias
    bias = 1023
    exponent = exponent + bias

    return (sign, bin(exponent)[2:].zfill(11), b_mantissa)

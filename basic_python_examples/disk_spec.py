#!/usr/bin/env

'''
Quick and dirty way to check how long it takes to access the disk on your machine.

>>> execfile('disk_spec.py')
>>> disk_spec()
Results
Read Speed	Write Speed	Seek Speed
6.277238	10.954050	4.388958
'''

import os
import random
import tempfile
import time

def disk_spec(bufsize=1024, runs=1000000):
    filename = tempfile.mkstemp()[1]
    writespeed = 0.0
    readspeed = 0.0
    seekspeed = 0.0
    
    now = time.time()
    with open(filename, 'w') as fh:
        
        for i in xrange(runs):
            fh.write('x'*bufsize)   # write bufsize bytes to the buffer
            fh.flush()              # flush the buffer
        writespeed = time.time() - now
        fh.close()
    
    now = time.time()
    with open(filename, 'r') as fh:    
        fh.seek(0, 0)               # seek to the beginning of the file
        for i in xrange(runs):
            fh.read(bufsize)        # read in bufsize bytes
        readspeed = time.time() - now
        fh.close()
        
    now = time.time()
    with open(filename, 'r') as fh:
        fh.seek(0, 0)               # seek to the beginning of the file
        for i in xrange(runs):
            fh.seek(random.randint(0,bufsize*runs), 0) # seek to a random position in the file offset from the bgeginning of the file
        seekspeed = time.time() - now
        fh.close()                  # close the file
        
    os.unlink(filename)         # unlink the file as it is no longer needed
    
    print 'Results'
    print 'Read Speed\tWrite Speed\tSeek Speed'
    print '{0:f}\t{1:f}\t{2:f}'.format(readspeed, writespeed, seekspeed)
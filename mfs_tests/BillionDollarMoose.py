#!/usr/bin/env python
'''
Creates file system objects for moose testing
'''
import argparse
import collections
import os
import sys
import time

def create_directories(total_num, per_directory, basedir):
    '''Create a directory tree with per_directory directories per directory and a total number of total_num directories under basedir'''
    if not os.path.isdir(basedir):
        sys.stderr.write('{0:s} does not exist.\n'.format(basedir))
        sys.exit(1)
    # Keep a stack to traverse the directory structure we will create
    dirstack = collections.deque()
    # Push our top level directory onto the stack
    dirstack.append(basedir)
    # The current count of subdirectories created
    count = 0
    # While the stack is not empty, pop off a directory and make that the cwd
    while dirstack and count < total_num:
        cwd = dirstack.popleft()
        #print 'cd {0:s}'.format(cwd)
        os.chdir(cwd)
        # Create per_directory number of files in this directory and add those directories to the end of the stack
        if count + per_directory > total_num:
            end = total_num
        else:
            end = count + per_directory
        for i in xrange(count, end):
            count += 1
            subdir = os.path.join(cwd, str(i))
            #print 'mkdir {0:s}'.format(subdir)
            os.mkdir(subdir)
            dirstack.append(subdir)
    return count

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('total', type=int, help='The total number of file system directories to create')
    parser.add_argument('num_per_directory', type=int, help='The maximum number of top-level subdirectories each directory should have')
    parser.add_argument('basedir', type=str, help='The base directory to start creating directories under')
    args = parser.parse_args()
    mydir = os.getcwd()
    sys.stderr.write('Starting: {0:s}'.format(time.asctime()))
    start = time.time()
    count = create_directories(args.total, args.num_per_directory, args.basedir)
    now = time.time()
    os.chdir(mydir)
    sys.stderr.write('Finished: {0:s}\n'.format(time.asctime()))
    print 'Total number of directories created: {0:d}'.format(count)
    print 'Total elapsed time: {0:f}s'.format(now-start)
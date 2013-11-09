#!/usr/bin/env python

import datetime
import logging
import os
import os.path
import shlex
import shutil
import subprocess
import sys
import tempfile
import time

# Set up logging
logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s', level=logging.DEBUG)

def determine_fs_size(filesystem):
    # run df -B1 on the filesystem
    df_output = subprocess.check_output(['df','-B1',filesystem], shell=False)
    # creates a dictionary based on the df output for the given filesystem
    df_dict = dict(zip(df_output.split('\n')[0].split(),df_output.split('\n')[1].split()))
    # Sample output
    # {'Available': '1795633008640', 'Use%': '4%', 'Used': '56345862144', '1B-blocks': '1951088439296', 'Filesystem': '/dev/sdc1', 'Mounted': '/'}
    return(int(df_dict['Available']))

def choose_file_size(fs_size):
    '''choose_file_size(fs_size): choose file size based on total filesystem size'''
    i=10
    while (fs_size / 2**i > 2**15):
        i += 1
    if i > 10:
        return(2**(i - 1))
    else:
        return(2**10)

def fill_fs(available, fsize):
    nfiles = available / fsize
    for i in xrange(0, nfiles):
        try:
            dd_cmd = 'dd if={0:s} of={1:s} bs={2:d} count=1'.format('/dev/zero', os.path.join(os.getcwd(),str(i)+'.dat'), fsize)
            errno = subprocess.call(shlex.split(dd_cmd), shell=False)
        except subprocess.CalledProcessError:
            logging.fatal('\'{0:s}\' has failed, exit code {1:d}'.format(dd_cmd, errno))
            sys.exit(errno)
    return nfiles

def clean_fs(nfiles):
    for i in xrange(0, nfiles):
        try:
            rm_cmd = 'rm {0:s}'.format(os.path.join(os.getcwd(),str(i)+'.dat'))
            errno = subprocess.call(shlex.split(rm_cmd), shell=False)
        except subprocess.CalledProcessError:
            logging.fatal('\'{0:s}\' has failed, exit code {1:d}'.format(rm_cmd, errno))
            sys.exit(errno)
    return nfiles

def usage():
    print 'usage: {0:s} filesystem'.format(sys.argv[0])
    print 'Where \'filesystem\' is a valid mount point.'
    print '(Check output of \'df\' if unsure.)'
    sys.exit(1)

def timestamp():
    return time.strftime('%Y%m%d-%H:%M:%S',time.localtime())

def cleanup_failure(tmpdir):
    logging.fatal('Cleanup appears to have failed, please examine {0:s}'.format(tmpdir))

def main():
    if len(sys.argv) != 2:
        usage()
    else:
        # Track how many times we fill and clean the filesystem
        num_runs = 0
        # Determine where to create our temp directory on the filesystem
        # Is sys.argv[1] a valid mount point?
        df_output = subprocess.check_output(['df'], shell=False)
        valid_mounts = [df_output.split('\n')[i].split()[-1] for i in xrange(1,len(df_output.split('\n'))-1)]
        mount = sys.argv[1]
        if mount in valid_mounts:
            try:
                # Create a temporary working directory under the root of the mount
                tmpdir = tempfile.mkdtemp(dir=sys.argv[1])
                # Set this newly created as the cwd
                os.chdir(tmpdir)
                # Determine the amount of freely available space
                available = determine_fs_size(mount)
                # Determine optimal file size based on available space
                fsize = choose_file_size(available)
                logging.info('Beginning stress tests...')
                logging.info('Available space in bytes: {0:d}'.format(available))
                logging.info('Optimal file size in bytes: {0:d}'.format(fsize))
                logging.info('Files to be created/removed per run: {0:d}'.format(available/fsize))
                logging.info('{0:s}'.format('-'*60))
                while True:
                    # Fill it
                    nfiles = fill_fs(available, fsize)
                    # Clean it
                    rfiles = clean_fs(nfiles)
                    # Increase the run count and print some info about the last run
                    num_runs += 1
                    logging.info('Run({0:d}): Created {1:d} files, Removed {2:d} files'.format(num_runs, nfiles, rfiles))
                    logging.info('Run({0:d}): {1:d} of {2:d} ({3:0.2f}%)'.format(num_runs, (fsize*nfiles), available, (fsize*nfiles)/float(available)))
            except Exception as e:
                logging.info('{0:d} runs, {1:d} bytes/run'.format(num_runs, fsize*nfiles))
            finally:
                logging.info('Cleaning up now...')
                # Set the cwd to root level of the mount
                os.chdir(mount)
                # Cleanup the temporary working directory
                shutil.rmtree(tmpdir,onerror=cleanup_failure(tmpdir))
                logging.info('Cleanup complete!')
        else:
            print '{0:s} is not a valid mount point.'.format(sys.argv[1])
            print 'Valid mount points are {0:s}'.format(valid_mounts)
            usage()

if __name__ == '__main__':
    main()

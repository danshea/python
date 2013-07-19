#!/usr/bin/env python
#
# Name: io.py
# Date: 2013-07-18
# Description: writes files to an mfs mount point so that we may conduct failover
# testing under some i/o
#

import ConfigParser
import logging
import os
import random
import re
import signal
import sys
import time

def signal_handler(signal, frame):
    logging.debug('Cleanup and exit.')
    logging.info('Removing files.')
    for i in xrange(num_files):
        if os.path.isfile(os.path.join(mountpoint,str(i)+'.dat')):
            try:
                os.unlink(os.path.join(mountpoint,str(i)+'.dat'))
                logging.info('{0:s} removed.'.format(os.path.join(mountpoint,str(i)+'.dat')))
            except Exception as e:
                logging.warn(e)
        else:
            logging.warn('{0:s} does not exist.'.format(os.path.join(mountpoint,str(i)+'.dat'))) 
    sys.exit(0) 

# Setup the signal handler
signal.signal(signal.SIGINT, signal_handler)

# Open the configfile and parse it
try:
    configfile = 'io.cfg'
    if not os.path.isfile(configfile):
        raise IOError
    try:
        open(configfile, 'r')
    except Exception as e:
        raise e
    config = ConfigParser.ConfigParser()
    config.read(configfile)
    CONFIG='Configuration'
    num_files = config.get(CONFIG, 'num_files')
    file_size = config.get(CONFIG, 'file_size')
    buff_size = config.get(CONFIG, 'buff_size')
    mountpoint = config.get(CONFIG, 'mountpoint')
    logfile_prefix = config.get(CONFIG, 'logfile_prefix')
    logdir = config.get(CONFIG, 'logdir')
except Exception as e:
    sys.stderr.write(e)
    sys.exit(1)

# Initialize the logger
if os.path.isdir(logdir):
    logfile = logfile_prefix + time.strftime('%Y%m%d%H%M%S', time.localtime()) + '.log'
    logging.basicConfig(filename=os.path.join(logdir,logfile), level=logging.DEBUG, format='%(asctime)s %(message)s')
    logging.debug('Logger Initialized.')
else:
    sys.stderr.write('{0:s} does not exist or is not a directory!'.format(logdir))
    sys.exit(1)

logging.info('num_files = {0:s}'.format(num_files))
logging.info('file_size = {0:s}'.format(file_size))
logging.info('buff_size = {0:s}'.format(buff_size))
logging.info('mountpoint = {0:s}'.format(mountpoint))
logging.info('logfile_prefix = {0:s}'.format(logfile_prefix))
logging.info('logdir = {0:s}'.format(logdir))
logging.debug('Validating configuration values.')
# Validate the configuration values
try:
    num_files = int(num_files)
    if num_files <= 0:
        raise ValueError('num_files must be a positive integer')
except Exception as e:
    logging.fatal(e)
    sys.exit(1)

try:
    buff_size = int(buff_size)
    if buff_size <= 0:
        raise ValueError('buff_size must be a positive integer')
except Exception as e:
    logging.fatal(e)
    sys.exit(1)

try:
    # See if the file size is valid
    intmatch = re.match('[0-9]+[kMG]', file_size)
    floatmatch = re.match('[0-9]+\.[0-9]+[kMG]', file_size)
    if intmatch:
        unit = file_size[len(file_size)-1]
        file_size = int(file_size[0:len(file_size)-1])
    elif floatmatch:
        unit = file_size[len(file_size)-1]
        file_size = float(file_size[0:len(file_size)-1])
    else:
        raise ValueError('{0:s} is not valid, please check the configuration file.'.format(file_size))
    
    # Account for unit
    if unit == 'k':
        file_size *= 2**10
    elif unit == 'M':
        file_size *= 2**20
    elif unit == 'G':
        file_size *= 2**30
    else:
        raise ValueError('{0:s} is not a valid unit, please check the configuration file.'.format(unit))
except Exception as e:
    logging.fatal(e)
    sys.exit(1)

try:
    if not os.path.isdir(mountpoint):
        raise IOError('{0:s} not found.'.format(mountpoint))
except Exception as e:
    logging.fatal(e)

logging.debug('Configuration values validated.')

logging.debug('Entering file creation loop.')
for i in xrange(num_files):
    try:
        fh = open(os.path.join(mountpoint,str(i)+'.dat'), 'w')
        buffcount = int(file_size / buff_size)
        while buffcount > 0:
            buff = buff_size * chr(random.randint(65,90))
            fh.write(buff)
            buffcount -= 1
        logging.info('{0:s} written.'.format(os.path.join(mountpoint,str(i)+'.dat')))
    except Exception as e:
        logging.fatal(e)

logging.debug('Cleanup and exit.')
logging.info('Removing files.')
for i in xrange(num_files):
    if os.path.isfile(os.path.join(mountpoint,str(i)+'.dat')):
        try:
            os.unlink(os.path.join(mountpoint,str(i)+'.dat'))
            logging.info('{0:s} removed.'.format(os.path.join(mountpoint,str(i)+'.dat')))
        except Exception as e:
            logging.warn(e)
    else:
        logging.warn('{0:s} does not exist.'.format(os.path.join(mountpoint,str(i)+'.dat'))) 
sys.exit(0)

#!/usr/bin/env python
#
# Name: mfs_install
# Date: 2013-07-14
# Author: djs
# Version: 0.1
# Description: Install moosefs from a source tar file.  Configurable options are located in mfs_install.cfg
# Usage: install_mfs <version> <master|chunkserver>
# version should be in the form of the version appended to the source code tar file, i.e. - 1.7.2
# Specifiy master if installing a master
# Specify chunkserver if installing a chunkserver
#
# Exit status codes:
# 0 - success, program completed successfully
# 1 - configuration file not found or not readable
# 2 - source file does not exist
# 3 - command line argument parsing failed
# 4 - problem extracting files
# 5 - tarfile extraction directory does not exist
# 6 - configure returned non-zero exit status
# 7 - make returned non-zero exit status
# 8 - make install returned non-zero exit status
# 911 - mfs_install root directory does not exist
# 912 - nodetype was set to invalid value

# import modules necessary for this program
import argparse
import ConfigParser
import logging
import os
import shlex
import subprocess
import sys
import tarfile
import time

def main():
    if __name__ == '__main__':
        # NOTE: This program should be executed from within the root directory or it will not find the configuration file
        # Parse the configuration and assign values
        try:
            CONFIG = 'mfs_install.cfg'
            configfile = open(CONFIG, 'r')
        except IOError as detail:
            sys.stderr.write('{0:s} does not exist or is not readable.'.format(CONFIG))
            sys.stderr.write(detail)
            sys.exit(1)
        else:
            with configfile:
                SECTION = 'Configuration'
                config = ConfigParser.ConfigParser()
                config.read(configfile)
                root = config.get(SECTION, 'root')
                sourcedirectory = os.path.join(root, config.get(SECTION, 'sourcedirectory'))
                sourcebasename = config.get(SECTION, 'sourcebasename')
                sourcefile_extension = config.get(SECTION, 'sourcefile_extension')
                tarmode = config.get(SECTION, 'tarmode')
                mfsuser = config.get(SECTION, 'mfsuser')
                mfsgroup = config.get(SECTION, 'mfsgroup')
                loglevel = config.get(SECTION, 'loglevel')
                logdir = os.path.join(root, config.get(SECTION, 'logdir'))
                logfile = config.get(SECTION, 'logfile')
                master_flags = config.get(SECTION, 'master_configure_flags')
                chunk_flags = config.get(SECTION, 'chunkserver_configure_flags')
        
        # Make sure the root exists, if it doesn't we have some big problems!
        if not os.path.isdir(root):
            sys.stderr.write('{0:s} does not exist, check the confguration file and the installation path.')
            sys.exit(911)
        
        # Ensure that logdir, sourcedirectory exist, create them if necessary
        if not os.path.isdir(os.path.join(root, logdir)):
            os.mkdir(os.path.join(root, logdir))
            sys.stderr.write('{0:s} did not exist, creating.'.format(os.path.join(root, logdir)))
        if not os.path.isdir(os.path.join(root, sourcedirectory)):
            os.mkdir(os.path.join(root, sourcedirectory))
            sys.stderr.write('{0:s} did not exist, creating, this run will fail due to source file not being present.'.format(os.path.join(root, sourcedirectory)))
            sys.stderr.write('Please place the source file into {0:s} and re-run this command.'.format(os.path.join(root, sourcedirectory)))
            exit(2)
        
        # Initialize the logger
        logfile += '-' + time.strftime('%Y%m%d%H%M%S', time.localtime()) + '.log'
        logging.basicConfig(filename=os.path.join(logdir,logfile), level=logging.DEBUG)
        logging.basicConfig(format='%(asctime)s %(message)s')
        logging.debug('Logger Initialized.')
        
        # Parse the command line and ensure we get the arguments we expect and that they are valid
        try:
            parser = argparse.ArgumentParser()
            parser.add_argument('version', help='version appended to the source code tar file, i.e. - 1.7.2')
            parser.add_argument('nodetype', help='specify \'master\' or \'chunkserver\' installation type', choices=['master','chunkserver'])
            args = parser.parse_args()
            version = args.version
            nodetype = args.nodetype
            sourcefilename = sourcebasename + version + sourcefile_extension
        except Exception as detail:
            logging.fatal(detail)
            sys.exit(3)
        
        logging.info('Command line arguments parsed.')
        
        # Ensure the sourcefile exists
        if not os.path.isfile(os.path.join(sourcedirectory, sourcefilename)):
            logging.fatal('{0:s} does not exist.'.format(os.path.join(sourcedirectory, sourcefilename)))
            sys.stderr.write('{0:s} does not exist.'.format(os.path.join(sourcedirectory, sourcefilename)))
            sys.exit(2)
        
        # untar the source file
        try:
            sourcetar = tarfile.open(os.path.join(sourcedirectory, sourcefilename), mode=tarmode)
            sourcetar.extractall(sourcedirectory)
        except tarfile.TarError as detail:
            logging.fatal('Error retrieving files from {0:s}'.format(os.path.join(sourcedirectory, sourcefilename)))
            sys.stderr.write('Error retrieving files from {0:s}'.format(os.path.join(sourcedirectory, sourcefilename)))
            logging.fatal(detail)
            sys.exit(4)
        
        logging.info('untar of source file complete.')
        
        # One last check to make sure the subdirectory we expect has been created, if it has, set the CWD to that directory
        if not os.path.isdir(os.path.join(sourcedirectory,sourcebasename+version)):
            logging.fatal('{0:s} is not a valid directory.  Please check the source tarfile was properly extracted.'.format(os.path.join(sourcedirectory,sourcebasename+version)))
            sys.stderr.write('{0:s} is not a valid directory.  Please check the source tarfile was properly extracted.'.format(os.path.join(sourcedirectory,sourcebasename+version)))
            sys.exit(5)
        else:
            os.chdir(os.path.join(sourcedirectory,sourcebasename+version))
            logging.info('CWD is now {0:s}'.format(os.path.join(sourcedirectory,sourcebasename+version)))
        
        # Run configure
        try:
            if nodetype == 'master':
                subprocess.check_call(shlex.split(master_flags))
            if nodetype == 'chunkserver':
                subprocess.check_call(shlex.split(chunk_flags))
            else:
                logging.fatal('nodetype = {0:s} this should never happen!  nodetype should be one of master or chunkserver.'.format(nodetype))
                sys.stderr.write('nodetype = {0:s} this should never happen!  nodetype should be one of master or chunkserver.'.format(nodetype))
                sys.exit(912)
        except subprocess.CalledProcessError:
            logging.fatal('configure returned non-zero exit status')
            sys.stderr('configure returned non-zero exit status')
            logging.fatal(subprocess.CalledProcessError.output)
            sys.exit(6)
        
        logging.info('Configuration complete.')
        
        # compile the source
        try:
            subprocess.check_call(shlex.split('make'))
        except subprocess.CalledProcessError:
            logging.fatal('make returned non-zero exit status')
            sys.stderr.write('make returned non-zero exit status')
            logging.fatal(subprocess.CalledProcessError.output)
            sys.exit(7)
        
        logging.info('Compilation complete.')
        
        # upon successful compilation of the source code install the software
        try:
            subprocess.check_call(shlex.split('make install'))
        except subprocess.CalledProcessError:
            logging.fatal('make install returned non-zero exit status')
            sys.stderr.write('make install returned non-zero exit status')
            logging.fatal(subprocess.CalledProcessError.output)
            sys.exit(8)
        
        logging.info('Installation complete.')
        
        # cleanup and exit
        os.chdir(root)
        logging.info('CWD is now {0:s}'.format(root))
        logging.info('Done.')
        sys.exit(0)
    else:
        pass
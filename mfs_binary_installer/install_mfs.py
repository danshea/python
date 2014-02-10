#!/usr/bin/env python
'''
Author     : djs
Name       : install_mfs.py
Date       : 2013-12-21
Description: Moose Binary Installer
             This program installs the MooseFS file system using vendor
             supplied binaries and returns appropriate json formatted
             success or failure codes to salt stack for automation of
             the installation.
'''

import gzip
import tarfile
import os

class mfsInstaller(object):
    def __init__(self, version):
        '''
        Initialize the installer with the appropriate version number as
        an argument to the __init__ method
        '''
        self.localsoftware = os.path.join('root','LocalSoftware')
        self.version = version
        self.filename = 'mfs-{0:s}-debs.tgz'.format(version)
    
    
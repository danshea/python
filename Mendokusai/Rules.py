#!/usr/bin/env python3
'''
Author: dshea
Date: 2013-11-17
Description: Rules.py

Rules are simple python class structures that contain a test method and a run method

test(filename, stats): must return True or False
run(filename): perform some action on the filename given

TODO:
* Since this was a proof of concept, look at ways to structure it better and lock
  down unintended behaviour.

* Add exception handling and modify the main body to be run
  from the command line.

* Implement logging module to log actions and create a noop flag that merely walks
  through the rules and logs what action would have been taken on a given file.
'''
import os
import pwd
import grp

class Rule0:
    def test(filename, stats):
        (st_mode,st_ino,st_dev,st_nlink,st_uid,st_gid,st_size,st_atime,st_mtime,st_ctime) = stats
        if st_uid == 1000:
            print('{0:s} is owned by user {1:s}'.format(filename,pwd.getpwuid(st_uid)[0]))
            return True
        else:
            return False
    def run(filename):
        print('{0:s} method called.'.format(__name__))
        pass

class Rule1:
    def test(filename, stats):
        (st_mode,st_ino,st_dev,st_nlink,st_uid,st_gid,st_size,st_atime,st_mtime,st_ctime) = stats
        if st_gid == 1000:
            print('{0:s} is owned by group {1:s}'.format(filename,grp.getgrgid(st_gid)[0]))
            return True
        else:
            return False
    def run(filename):
        print('{0:s} method called.'.format(__name__))
        pass
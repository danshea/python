#!/usr/bin/env python3
'''
Author: dshea
Date: 2013-11-17
Description: Mendokusai.py

Mendokusai walks a directory structure and performs actions on the
files it encounters based on rules defined in Rules.py

Rules are simple python class structures that contain a test method and a run method
If the test method evaluates to true, then the run method is called and an action is
taken on that file.

TODO:
* Since this was a proof of concept, look at ways to structure it better and lock
  down unintended behaviour.

* Add exception handling and modify the main body to be run
  from the command line.

* Implement logging module to log actions and create a noop flag that merely walks
  through the rules and logs what action would have been taken on a given file.
'''
import inspect
import os
import pwd
import sys
import collections
import Rules

def run(start=os.getcwd()):
    '''run(start=os.getcwd()): Walk directory tree starting from start and apply rules'''
    # Walk the directory tree
    for (dirpath,dirnames,filenames) in os.walk(start):
        # Create an OrderedDict of files encountered, calling os.stat on each file
        file_stats = collections.OrderedDict()
        for filename in filenames:
            absolute_filename = os.path.join(dirpath,filename)
            file_stats[absolute_filename] = os.stat(absolute_filename)
            # Now we have an OrderedDict of every file at this directory level, process
            # the files based on some pre-defined rules written as python classes and
            # imported into this file as a list of those classes, if the test() method
            # evaluates to True, the action() method is then called against that file
            for (rule_name, rule) in inspect.getmembers(Rules, inspect.isclass):
                if rule.test(absolute_filename, file_stats[absolute_filename]) == True:
                    rule.run(absolute_filename)
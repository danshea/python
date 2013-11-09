#!/usr/bin/env python3

'''
Author: djs
Date: 2012-07-11
Description: Simple example of a Read Eval Print Loop
'''

import sys

def REPL():
    print('Simple shell, \\ is the line continuation character.')
    print('Ctrl-C exits')
    try:
        lines = ''
        print('$ ', end='')
        while True:
            line = sys.stdin.readline().strip()
            if len(line) and line[-1] == '\\':
                lines += line
                print('> ', end='')
                continue
            else:
                lines += line
                lines = lines.replace('\\','\n')
                print('received: {0}'.format(lines))
                lines = ''
                print('$ ', end='')  
    except KeyboardInterrupt:
        print('Goodbye.')

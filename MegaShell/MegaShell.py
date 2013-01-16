#!/usr/bin/env python
#
# Name: MegaShell
# Description: Shell wrapper for MegaCli that tries to be slightly more user
# friendly.
# Date: 2012-01-15
# Author: Dan Shea, Harvard Medical School
# Email: daniel_shea2@hms.harvard.edu
#
import argparse
import cmd
import ConfigParser
import datetime
import logging
import os
import os.path
import re
import shlex
import subprocess
import sys
import time

class MegaShell(cmd.Cmd):
    '''MegaShell command line interpreter'''
    def __init__(self, cfgfile='conf/MegaShell.conf',
                 logfile='log/MegaShell_{}.log'.format(time.strftime('%Y%m%d%H%M%S',
                                                                  time.localtime())),
                 noexec=False):
        cmd.Cmd.__init__(self)
        # This flag determines if the actual commands are executed
        self._noexec = noexec
        # creates self.logger for logging
        self._setup_logging(logfile)
        # creates dict self.configuration and loads cfgfile
        self._read_config(cfgfile)
        # Get the hostname
        self._hostname = subprocess.check_output('hostname').strip()
        # Set the prompt
        self.prompt = self.configuration['ps1'].format(hostname=self._hostname)
        # Set the intro banner
        if self._noexec:
            self.intro = self.configuration['banner'].format(current_time=time.asctime()) + \
            ' (_noexec)'
        else:
            self.intro = self.configuration['banner'].format(current_time=time.asctime())
        self._cmd_prefix = os.path.join(self.configuration['megaclipath'],self.configuration['megaclicmd'])


    def _setup_logging(self, filename):
        # Setup the logger to log the session activity
        # create logger with 'MegaShell'
        self.logger = logging.getLogger('MegaShell')
        self.logger.setLevel(logging.DEBUG)
        # create file handler which logs even debug messages
        fh = logging.FileHandler(filename)
        fh.setLevel(logging.DEBUG)
        # create console handler with a higher log level
        ch = logging.StreamHandler()
        ch.setLevel(logging.ERROR)
        # create formatter and add it to the handlers
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        # add the handlers to the logger
        self.logger.addHandler(fh)
        self.logger.addHandler(ch)

    def _read_config(self, filename):
        self.configuration = dict()
        try:
            config = ConfigParser.ConfigParser()
            config.read(filename)
            sections = config.sections()
            for section in sections:
                self.logger.debug('Loading configuration section {}:'.format(section))
                for option in config.options(section):
                    self.configuration[option] = config.get(section, option)
                    self.logger.debug('Loaded configuration {} = {} from {}'.format(option,
                                                                               self.configuration[option],
                                                                               filename))
                self.logger.debug('End loading configuration section {}'.format(section))
        except Exception as e:
            self.logger.fatal('{}'.format(e))
    # Helper function calls MegaCli commands
    def _run(self, cmd):
        self.logger.debug('Executing: {}'.format(cmd))
        try:
            output = subprocess.check_output(shlex.split(cmd))
            return output
        except subprocess.CalledProcessError as e:
            self.logger.error('{}'.format(e))
            return None

    # Commands available to the user in MegaShell
    def do_AlarmOff(self, adapter):
        '''AlarmOff [adapter]
        Turn the alarm off for the adapter, if adapter is not specified,
        turn off the alarm for all adapters'''
        self.logger.debug('AlarmOff called.')
        if adapter:
            cmd = self.configuration['alarmoff'].format(adapter=adapter)
        else:
            cmd = self.configuration['alarmoff'].format(adapter='ALL')
        cmd = self._cmd_prefix + ' ' + cmd
        if not self._noexec:
            output = self._run(cmd)
            if output:
                sys.stdout.write('{}'.format(output))
                self.logger.debug('Output: {}'.format(output))
            else:
                self.logger.debug('Output:')
        else:
            self.logger.debug('Command: {}'.format(cmd))

    def do_ListDrives(self, adapter):
        '''ListDrives [adapter]
        List the physical drives for the adapter, along with a brief description,
        if adapter is not specified, list all drives for all adapters
        '''
        self.logger.debug('ListDrives called.')
        if adapter:
            cmd = self.configuration['listdrives'].format(adapter=adapter)
        else:
            cmd = self.configuration['listdrives'].format(adapter='ALL')
        cmd = self._cmd_prefix + ' ' + cmd
        if not self._noexec:
            output = self._run(cmd)
            if output:
                sys.stdout.write('{}'.format(output))
                self.logger.debug('Output: {}'.format(output))
            else:
                self.logger.debug('Output:')
        else:
            self.logger.debug('Command: {}'.format(cmd))

    def do_ListAdapters(self, line):
        '''ListAdapters
        List the adapter id for all adapters, along with a brief description'''
        self.logger.debug('ListAdapters called.')
        cmd = self.configuration['listadapters']
        cmd = self._cmd_prefix + ' ' + cmd
        if not self._noexec:
            output = self._run(cmd)
            if output:
                adapters = re.findall('Adapter #\d+', output)
                descriptions = re.findall('Product Name.+', output)
                for adapter, description in zip(adapters, descriptions):
                    adapter = adapter.strip()
                    description = description.strip()
                    sys.stdout.write('{}:{}\n'.format(adapter,description))
                    self.logger.debug('Output: {}:{}'.format(adapter,description))
            else:
                self.logger.debug('Output:')
        else:
            self.logger.debug('Command: {}'.format(cmd))

    def do_LED(self, args):
        '''LED on|off adapter enclosure slot
        Turn the drive activity LED on or off for a physical drive'''
        self.logger.debug('LED called')
        args = args.split()
        if len(args) != 4:
            sys.stdout.write('usage: LED on|off adapter enclosure slot\n')
            self.logger.debug('LED called with improper arguments')
            pass
        else:
            led, adapter, enclosure, slot = args
            if led == 'on' or led == 'off':
                cmd = self.configuration['usediskactivityforlocate']
                cmd = self._cmd_prefix + ' ' + cmd
                if not self._noexec:
                    output = self._run(cmd)
                else:
                    self.logger.debug('Command: {}'.format(cmd))
                if led == 'on':
                    cmd = self.configuration['ledon'].format(adapter=adapter, enclosure=enclosure, slot=slot)
                    cmd = self._cmd_prefix + ' ' + cmd
                    if not self._noexec:
                        output = self._run(cmd)
                        sys.stdout.write('Adapter: {} Enclosure: {} Drive Slot: {} LED: on'.format(adapter,
                                                                                                   enclosure,
                                                                                                   slot))
                        self.logger.debug('Output: Adapter: {} Enclosure: {} Drive Slot: {} LED: on'.format(adapter,
                                                                                                            enclosure,
                                                                                                            slot))
                    else:
                        self.logger.debug('Command: {}'.format(cmd))
                else:
                    cmd = self.configuration['ledoff'].format(adapter=adapter, enclosure=enclosure, slot=slot)
                    cmd = self._cmd_prefix + ' ' + cmd
                    if not self._noexec:
                        output = self._run(cmd)
                        sys.stdout.write('Adapter: {} Enclosure: {} Drive Slot: {} LED: off'.format(adapter,
                                                                                                    enclosure,
                                                                                                    slot))
                        self.logger.debug('Output: Adapter: {} Enclosure: {} Drive Slot: {} LED: off'.format(adapter,
                                                                                                             enclosure,
                                                                                                             slot))
                    else:
                        self.logger.debug('Command: {}'.format(cmd))
            else:
                sys.stdout.write('usage: LED on|off adapter enclosure slot\n')
                self.logger.debug('LED called without on or off argument')
                pass

    def do_exit(self, line):
        '''Exits Megashell'''
        self.logger.debug('do_exit called, exiting.')
        return True

    def do_EOF(self, line):
        '''EOF drops us out of the Megashell'''
        self.logger.debug('do_EOF called, exiting.')
        return self.do_exit(line)

    # These definitions override cmd.Cmd methods
    def emptyline(self):
        '''Do nothing on an empty input'''
        pass

def main():
    # We need to parse any commands supplied on the command line when MegaShell is
    # invoked.
    # Create the argument parser
    parser = argparse.ArgumentParser()
    # Add optional arguments
    parser.add_argument('-c', '--conf', help='the configuration file to use')
    parser.add_argument('-l', '--log', help='the log file to use')
    parser.add_argument('-n', '--noexec', help='do not execute commands', action="store_true")
    args = parser.parse_args()
    # Instantiate a command interpreter and run it based on these options
    if args.noexec:
        noexec = True
    else:
        noexec = False
    if args.log and args.conf:
        megashell = MegaShell(cfgfile=args.conf, logfile=args.log, noexec=noexec)
    elif args.conf:
        megashell = MegaShell(cfgfile=args.conf, noexec=noexec)
    elif args.log:
        megashell = MegaShell(logfile=args.log, noexec=noexec)
    else:
        megashell = MegaShell(noexec=noexec)
    exit_code = megashell.cmdloop()
    sys.exit(exit_code)

if __name__ == '__main__':
    main()

#!/usr/bin/env python
import ConfigParser
import logging
import os
import paramiko
import random
import re
import shlex
import signal
import subprocess
import sys
import time

class ChaosMoose(object):
    def __init__(self, mastername='moose1-master.med.harvard.edu'):
        """Initialize the Chaos Moose! Optionally pass it the shared name of master nodes."""
        # Initialize Logging
        logfile = '.'.join(['-'.join(['ChaosMoose', time.strftime('%Y%m%d%H%M%S')]), 'log'])
        logging.basicConfig(filename=logfile, level=logging.DEBUG, format='%(asctime)s %(levelname)s: %(message)s')
        self.mastername = mastername
        self.method = None
        self.method_args = None
        logging.debug('mastername = {0:s}'.format(mastername))
        # Instantiate a paramiko ssh object
        self.ssh = paramiko.SSHClient()
        logging.debug('instantiated paramiko SSHClient object')
        # Tell the paramiko object to accept unknown keys
        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        logging.debug('set paramiko host key policy to AutoAddPolicy')

    def ssh_command(self, host, command):
        logging.info('ssh_command({0:s}, {1:s})'.format(host, command))
        self.ssh.connect(host)
        logging.debug('Connected to {0:s}'.format(host))
        stdin, stdout, stderr = self.ssh.exec_command(command)
        output = stdout.readlines()
        error = stdout.readlines()
        logging.debug('stdout: {0:s}'.format(''.join(output)))
        logging.debug('stderr: {0:s}'.format(''.join(error)))
        self.ssh.close()
        logging.debug('Session to {0:s} closed'.format(host))
        return (output, error)

    def current_leader(self):
        """Return host and port tuple of the current elected leader of the master nodes or (None, None) if there is no current leader."""
        executable = '/usr/sbin/mfshastarter'       # To be moved to a configuration file
        options = ' '.join(['-h', self.mastername]) # To be moved to a configuration file
        command = ' '.join([executable, options])
        # ssh to one of the masters to find the leader
        output, error = self.ssh_command(self.mastername, command)
        leader = None
        port = None
        for line in output:
            # If we encounter a DEAD master node, return (None, None) this should avoid inadvertently killing
            # more than one master node in the event that a master fails to recover, log the dead master's ip
            # in the log file
            deadmatch = re.search('DEAD', line)
            if deadmatch:
                leader,port = deadmatch.string.split(' : ')[0].split()[1].strip('(').strip(')').split(':')
                logging.debug('master node {0:s}:{1:s} is DEAD'.format(leader, port))
                return (None, None)
            match = re.search('LEADER', line)
            if match:
                # strip off extraneous bits, leaving us with the ip and the port of the master node
                leader,port = match.string.split(' : ')[0].split()[1].strip('(').strip(')').split(':')
        return (leader, port)

    '''Methods of execution go here'''
    
    def _clean_kill(self, leader, port):
        """Kill the leader via a clean shutdown, wait a bit and then restore the service."""
        executable = '/usr/sbin/mfsmaster' # To be moved to a configuration file
        try:
            options = 'stop'                   # To be moved to a configuration file
            command = ' '.join([executable, options])
            output, error = self.ssh_command(leader, command)
            self.random_sleep()
        except Exception as e:
            sys.stderr.write(e.message)
        finally:
                options = 'start'
                command = ' '.join([executable, options])
                output, error = self.ssh_command(leader, command)

    #TODO: Fix this so that the ssh interface is not the interface being brought down
    def _network_outage(self, leader, port):
        """Kill the leader by simulating a network outage, wait a bit and then restore the network."""
        try:
            executable = '/sbin/ifdown' # To be moved to a configuration file
            options = 'eth0'            # To be moved to a configuration file
            command = ' '.join([executable, options])
            output, error = self.ssh_command(leader, command)
            self.random_sleep()
        except Exception as e:
            sys.stderr.write(e.message)
        finally:
            # The interface should always be brought back up, no matter what
            executable = '/sbin/ifup' # To be moved to a configuration file
            command = ' '.join([executable, options])
            output, error = self.ssh_command(leader, command)
    
    '''End of methods of execution'''
    
    def kill_leader(self, leader, port):
        """Kill the current leader using a randomly chosen method of execution."""
        method_of_execution = [self._clean_kill] # To be moved to a configuration file
        self.method = random.choice(method_of_execution)
        self.method_args=(leader, port)
        self.method(*self.method_args)

    def monitor_election(self):
        """Monitor the status of the election process and return the new leader and port tuple when a new leader has been elected."""
        leader = None
        while leader == None:
            leader, port = self.current_leader()
            self.random_sleep(sleep_min=1, sleep_max=5)
        return (leader, port)
    
    def random_sleep(self, sleep_min=60, sleep_max=300):
        """Sleep for a random period of time defined by an inclusive sleep range, given in seconds."""
        sleepval = random.randint(sleep_min, sleep_max)
        try:
            logging.info('Going to sleep for {0:d}s'.format(sleepval))
            time.sleep(sleepval)
        except Exception as e:
            sys.stderr.write(e.message)
        
    def signal_handler(self, signal, frame):
        if self.method != None:
            self.method(*self.method_args)
        sys.stdout.write('\nChaos Moose has been contained, but for how long?!?!?\n')
        sys.stdout.write('Exiting.\n')
        sys.exit(0)
    
    def run(self):
        """Commit random acts of chaos against the master servers until an interrupt signal is received."""
        chaos_moose='''
            WWWWWW||WWWWWW
             W W W||W W W
                  ||
                ( OO )__________
                 /  |           \\
                /o o|W           \\
                \___/||_||__||_|| *
                     || ||  || ||
                    _||_|| _||_||
                   (__|__|(__|__|
        '''
        sys.stdout.write(chaos_moose)
        sys.stdout.write('Chaos Moose has been unleashed!!!!!\n')
        sys.stdout.write('      To end his reign of terror, press Ctrl+C\n')
        signal.signal(signal.SIGINT, self.signal_handler)
        while True:
            # Monitor the election until a new leader has been chosen
            leader, port = self.monitor_election()
            # Kill the current leader
            self.kill_leader(leader, port)
            # Sleep for a random interval of time
            self.random_sleep()
            # Continue to do this until Chaos Moose is sent an Interrupt Signal sending him back to the murky depths from whence he came
    
if __name__ == '__main__':
    cm = ChaosMoose()
    cm.run()
    sys.exit(0)
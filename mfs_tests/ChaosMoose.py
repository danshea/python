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
        logging.debug('mastername = {0:s}'.format(mastername))
        # Instantiate a paramiko ssh object
        self.ssh = paramiko.SSHClient()
        logging.debug('instantiated paramiko SSHClient object')
        # Tell the paramiko object to accept unknown keys
        self.ssh.set_mising_host_key_policy(paramiko.AutoAddPolicy())
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
        return output

    def current_leader(self):
        """Return host and port tuple of the current elected leader of the master nodes or (None, None) if there is no current leader."""
        executable = '/usr/sbin/mfshastarter'       # To be moved to a configuration file
        options = ' '.join(['-h', self.mastername]) # To be moved to a configuration file
        command = ' '.join([executable, options])
        # ssh to one of the masters to find the leader
        output = self.ssh_command(self.mastername, command)
        leader = None
        port = None
        for line in output:
            match = re.search('.LEADER.', output)
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
            output = self.ssh_command(leader, command)
            self.random_sleep()
        except Exception as e:
            sys.stderr.write(e)
        finally:
            # Check if the mfsmaster process is down, if it is down, restart it
            try:
                options = 'test'
                command = ' '.join([executable, options])
                output = self.ssh_command(leader, command)
                for line in output:
                    match = re.search('mfsmaster is not running')
                    if match:
                        options = 'start'
                        command = ' '.join([executable, options])
                        output = self.ssh_command(leader, command)
            except Exception as e:
                sys.stderr.write(e)

    def _network_outage(self, leader, port):
        """Kill the leader by simulating a network outage, wait a bit and then restore the network."""
        try:
            executable = '/sbin/ifdown' # To be moved to a configuration file
            options = 'eth0'            # To be moved to a configuration file
            command = ' '.join([executable, options])
            output = self.ssh_command(leader, command)
            self.random_sleep()
        except Exception as e:
            sys.stderr.write(e)
        finally:
            # The interface should always be brought back up, no matter what
            try:
                executable = '/sbin/ifup' # To be moved to a configuration file
                command = ' '.join([executable, options])
                output = self.ssh_command(leader, command)
            except Exception as e:
                sys.stderr.write(e)
    
    '''End of methods of execution'''
    
    def kill_leader(self, leader, port):
        """Kill the current leader using a randomly chosen method of execution."""
        method_of_execution = [self._clean_kill, self._network_outage] # To be moved to a configuration file
        random.choice(method_of_execution)(leader, port)

    def monitor_election(self):
        """Monitor the status of the election process and return the new leader and port tuple when a new leader has been elected."""
        while leader == None:
            leader, port = self.current_leader()
        return (leader, port)
    
    def random_sleep(self, sleep_min=60, sleep_max=3600):
        """Sleep for a random period of time defined by an inclusive sleep range, given in seconds."""
        time.sleep(random.randint(sleep_min, sleep_max))

    def signal_handler(signal, frame):
        sys.stdout.write('Chaos Moose has been contained, but for how long?!?!?\n')
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
        sys.stdout.write('To end his reign of terror, press Ctrl+C\n')
        signal.signal(signal.SIGINT, signal_handler)
        while True:
            # Monitor the election until a new leader has been chosen
            leader, port = self.monitor_election()
            # Kill the current leader
            self.kill_leader(leader, port)
            # Sleep for a random interval of time
            self.random.sleep()
            # Continue to do this until Chaos Moose is sent an Interrupt Signal sending him back to the murky depths from whence he came
    
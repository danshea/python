#!/usr/bin/env python

'''
Author: dshea
Date: 2012-06-20
Name: ssh_keymanager.py
Description: Application that will create ssh keys for use in authorized_keys file and store them in a sqlite db file, lookup, key_comparison and basic key_management functions will be provided via a web interface
'''

import sqlite3
import os
import ConfigParser
import logging
import subprocess
import sys
import tempfile
import time
import datetime
import md5
import sha
import base64

class KeyManager():

    def __init__(self, cfgfile='ssh_keymanager.conf', dbfile='authorized_keys.db', logfile='ssh_keymanager.log'):
        
        # The configuration file and associated variables for section names
        self.cfgfile          = cfgfile
        self._dbsection       = 'Database'
        self._dbfilekey       = 'dbfile'
        self._privatetablekey = 'private_key_table'
        self._publictablekey  = 'public_key_table'
        self._hasherkey       = 'hash_function'
        
        self._logsection   = 'Logging'
        self._logfilekey   = 'logfile'
        self._loglevelkey  = 'loglevel'
        
        self._sqlsection    = 'SQl Statements'
        self._createSQLkey  = 'createSQL'
        self._insertSQLkey  = 'insertSQL'
        self._updateSQLkey  = 'updateSQL'
        self._deleteSQLkey  = 'deleteSQL'
        
        self._sshsection   = 'SSH'
        self._keygenkey    = 'ssh_keygen_cmd'
        
        # WARNING: If a configuration exists, it will OVERRIDE these arguments!!!
        #          If no configuration exists, these will become the defaults in the newly create configuration file!!!
        self.dbfile       = dbfile
        self.hasher       = 'sha'
        self.logfile      = logfile
        self.loglevel     = 'DEBUG'
        self.privatetable = 'private_keys'
        self.publictable  = 'public_keys'
        
        try:
            self.keygencmd = subprocess.check_output('which ssh-keygen', shell=True).strip()
        except subprocess.CalledProcessError:
            sys.stderr.write('you must have ssh-keygen installed to use this program')
            sys.exit(1)
        
        # Read in the configuration file if it exists, otherwise use the defaults and create a configuration
        # file with the defaults
        self.config = ConfigParser.ConfigParser()
        if os.path.isfile(self.cfgfile) == False:
            '''
            It took me a while to figure this one out.  One cannot paramterize database identifiers, so we parameterize them using
            python's string formatting and then use the recommended parameterization technique for literal values.
            '''
            self.createSQL = 'CREATE TABLE {0:s} (username text, key text, digest text, timestamp text)'
            self.insertSQL = 'INSERT INTO {0:s} VALUES (?, ?, ?, ?)'
            self.updateSQL = 'UPDATE {0:s} SET {1:s} = ? WHERE {2:s} = ?'
            self.deleteSQL = 'DELETE FROM {0:s} WHERE {1:s} = ?'
            self.selectSQL = 'SELECT * FROM {0:s} WHERE {1:s} = ?'
            
            self.config.add_section(self._dbsection)
            self.config.set(self._dbsection, self._dbfilekey, self.dbfile)
            self.config.set(self._dbsection, self._privatetablekey, self.privatetable)
            self.config.set(self._dbsection, self._publictablekey, self.publictable)
            self.config.set(self._dbsection, self._hasherkey, self.hasher)
            
            self.config.add_section(self._logsection)
            self.config.set(self._logsection, self._logfilekey, self.logfile)
            self.config.set(self._logsection, self._loglevelkey, self.loglevel)
            
            self.config.add_section(self._sqlsection)
            self.config.set(self._sqlsection, self._createSQLkey, self.createSQL)
            self.config.set(self._sqlsection, self._insertSQLkey, self.insertSQL)
            self.config.set(self._sqlsection, self._updateSQLkey, self.updateSQL)
            self.config.set(self._sqlsection, self._deleteSQLkey, self.deleteSQL)
            
            self.config.add_section(self._sshsection)
            self.config.set(self._sshsection, self._keygenkey, self.keygencmd)
            
            with open(self.cfgfile, 'w') as fh:
                self.config.write(fh)
        else:
            self.config.read(cfgfile)
            self.dbfile       = self.config.get(self._dbsection,  self._dbfilekey)
            self.privatetable = self.config.get(self._dbsection,  self._privatetablekey)
            self.publictable  = self.config.get(self._dbsection,  self._publictablekey)
            self.hasher       = self.config.get(self._dbsection,  self._hasherkey)
            self.logfile      = self.config.get(self._logsection, self._logfilekey)
            self.loglevel     = self.config.get(self._logsection, self._loglevelkey)
            self.createSQL    = self.config.get(self._sqlsection, self._createSQLkey)
            self.insertSQL    = self.config.get(self._sqlsection, self._insertSQLkey)
            self.updateSQL    = self.config.get(self._sqlsection, self._updateSQLkey)
            self.deleteSQL    = self.config.get(self._sqlsection, self._deleteSQLkey)
            self.keygencmd    = self.config.get(self._sshsection, self._keygenkey)
        
        # Set up the logger
        logging.basicConfig(filename=self.logfile, filemode='w', format='%(asctime)s %(message)s', level=getattr(logging, self.loglevel.upper()))
        logging.info('Logging Started')
        
        # if dbfile doesn't exist, create it at initialization
        if os.path.isfile(self.dbfile) == False:
            with open(self.dbfile, 'w') as fh:
                fh.write('')
                fh.close()
                # Log that we have created the file
                logging.info('created dbfile: ' + self.dbfile)
            # Create the private and public tables to store the keys in and log these events
            self.conn = sqlite3.connect(self.dbfile)
            logging.debug('db connection opened')
            self.curs = self.conn.cursor()
            logging.debug('db cursor created')
            self.curs.execute(self.createSQL.format(self.privatetable))
            logging.info('created db table {0:s}, statement executed: {1:s}'.format(self.privatetable, self.createSQL.format(self.privatetable)))
            self.curs.execute(self.createSQL.format(self.publictable))
            logging.info('created db table {0:s}, statement executed: {1:s}'.format(self.publictable, self.createSQL.format(self.publictable)))
            self.conn.commit()
            logging.debug('db changes committed')
            self.conn.close()
            logging.debug('db connection closed')
        logging.info('KeyManager object successfully initialized')

    def connectdb(self):
        '''connectdb(): connect to a sqlite db for storing keys and create a cursor'''
        self.conn = sqlite3.connect(self.dbfile)
        logging.debug('db cursor created')
        self.curs = self.conn.cursor()
        logging.debug('db cursor created')

    def disconnectdb(self,commit=True):
        '''disconnectdb(): close a connection to a sqlite db'''
        if commit == True:
            self.conn.commit()
            logging.debug('db changes committed')
        self.conn.close()
        logging.debug('db connection closed')

    def _getkey(self, filename):
        with open(filename, 'r') as fh:
            #return reduce(lambda x,y: x+y, fh.readlines())
            return ''.join(fh.readlines())

    def _mkdigest(self, key):
        '''_mkdigest(key): return base16 encoding of SHA-1 digest of the private key'''
        s = eval(self.hasher).new(key) # This will be used as a unique serial for each entry to join both public and private tables
        return base64.b16encode(s.digest())

    def createkeys(self, username=''):
        '''createkeys(username): create an ssh key pair and store them in the db'''
        tmpprivate = tempfile.mkstemp()[1]
        os.unlink(tmpprivate) # unlink the temp file or ssh-keygen will prompt to overwrite and program will hang waiting for input
        tmppublic  = tmpprivate+'.pub'
        cmd = self.keygencmd+' -q -b 2048 -N \'\' -C \'ssh key for '+username+'\' -t rsa -f '+tmpprivate
        subprocess.check_output(cmd, shell=True)
        private_key = self._getkey(tmpprivate)
        public_key  = self._getkey(tmppublic)
        digest = self._mkdigest(private_key) # digest of private key
        timestamp = str(time.time()) # POSIX timestamp
        self.connectdb()
        self.curs.execute(self.insertSQL.format(self.privatetable), (username, private_key, digest, timestamp))
        logging.info('inserted key {0:s} into {1:s}'.format(private_key, self.privatetable))
        self.curs.execute(self.insertSQL.format(self.publictable), (username, public_key, digest, timestamp))
        logging.info('inserted key {0:s} into {1:s}'.format(public_key, self.publictable))
        self.disconnectdb()
        # remove the temporary files mkstemp created
        os.unlink(tmpprivate)
        os.unlink(tmppublic)
    
    def getkeys(self, username='', private=False):
        '''
        getkey(username): return the key for username or None if not present
        default behaviour is to return the public key, setting private=True will return the private key
        '''
        self.connectdb()
        results = []
        if private == True:
            self.curs.execute(self.selectSQL.format(self.privatetable, 'username'), (username,))
            results = [row for row in self.curs]
        else:
            self.curs.execute(self.selectSQL.format(self.publictable, 'username'), (username,))
            results = [row for row in self.curs]
        self.disconnectdb()
        return results
    
    def removekeypair(self, digest):
        '''
        removekeypair(digest): remove all instances of the key matching the digest from the db
        WARNING: this function removes the key pair, both the private and the public keys are deleted!!!
        '''
        self.connectdb()
        private_rows = self.curs.execute(self.deleteSQL.format(self.privatetable, 'digest'), (digest,)).rowcount
        public_rows  = self.curs.execute(self.deleteSQL.format(self.publictable, 'digest'), (digest,)).rowcount
        logging.info('deleted {0:d} rows from {1:s} with digest={2:s}'.format(private_rows, self.privatetable, digest))
        logging.info('deleted {0:d} rows from {1:s} with digest={2:s}'.format(public_rows, self.publictable, digest))
        self.disconnectdb()
        return (private_rows, public_rows)
    
    def checkkey(self, key, private=False):
        '''checkkey(key): check if the key exists in the db'''
        self.connectdb()
        if private==True:
            self.curs.execute(self.selectSQL.format(self.privatetable, 'key'), (key,))
        else:
            self.curs.execute(self.selectSQL.format(self.publictable, 'key'), (key,))
        results = [row for row in self.curs]
        self.disconnectdb()
        if len(results) != 0:
            return True
        else:
            return False
        
if __name__ == '__main__':
    user = 'dshea'
    km = KeyManager()
    km.createkeys(user)
    public = km.getkeys(user)[0][1]
    private, digest = km.getkeys(user, private=True)[0][1:3]
    print public
    print private
    print km.checkkey(public)
    print km.checkkey(private, private=True)
    print km.removekeypair(digest)
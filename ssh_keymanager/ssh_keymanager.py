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
            self.createSQL = 'CREATE TABLE ? (username text, key text, digest text, timestamp text)'
            self.insertSQL = 'INSERT INTO ? VALUES (?, ?, ?, ?)'
            self.updateSQL = 'UPDATE ? SET ? = ? WHERE ? = ?'
            self.deleteSQL = 'DELETE FROM ? WHERE ? = ?'
            self.selectSQL = 'SELECT * FROM ? WHERE ? = ?'
            
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
            self.privatetable = self.config.get(self._dbsection,  self._dbprivatetablekey)
            self.publictable  = self.config.get(self._dbsection,  self._dbpublictablekey)
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
            self.curs.execute(self.createSQL,(self.privatetable))
            logging.info('created db table '+self.privatetable+', statement executed: ' + self.createSQL)
            self.curs.execute(self.createSQL,(self.publictable))
            logging.info('created db table '+self.publictable+', statement executed: ' + self.createSQL)
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
        tmppublic  = tmpprivate+'.pub'
        cmd = [self.keygencmd, ' -q ', ' -b ', '2048', ' -C ', '\'ssh key for '+username+'\'',
               ' -t ', 'rsa', ' -f ', tmpfile]
        subprocess.check_output(cmd, shell=False)
        private_key = self._getkey(tmpprivate)
        public_key  = self._getkey(tmppublic)
        digest = self._mkdigest(private_key) # digest of private key
        timestamp = time.time() # POSIX timestamp
        self.connectdb()
        self.curs.execute(self.insertSQL, (self.privatetable, username, key, digest, timestamp))
        self.curs.execute(self.insertSQL, (self.publictable, username, key, digest, timestamp))
        self.disconnectdb()
        # remove the temporary files mkstemp created
        os.unlink(tmpfile)
        os.unlink(tmpfile+'.pub')
    
    def getkeys(self, username='', private=False):
        '''
        getkey(username): return the key for username or None if not present
        default behaviour is to return the public key, setting private=True will return the private key
        '''
        self.connectdb()
        if private == True:
            self.curs.execute(self.selectSQL, (self.privatetable, 'username', username))
        else:
            self.curs.execute(self.selectSQL, (self.publictable, 'username', username))
        results = [row for row in self.curs]
        self.disconnectdb()
        return results
    
    def removekeypair(self, digest):
        '''
        removekeypair(digest): remove all instances of the key matching the digest from the db
        WARNING: this function removes the key pair, both the private and the public keys are deleted!!!
        '''
        self.connectdb()
        private_rows = self.curs.execute(self.deleteSQL, (self.privatetable, 'digest', digest)).rowcount
        public_rows  = self.curs.execute(self.deleteSQL, (self.publictable, 'digest', digest)).rowcount
        logging.info('deleted '+private_rows+' rows from '+self.privatetable+'with digest='+digest)
        logging.info('deleted '+public_rows+' rows from'+self.publictable+'with digest='+digest)
        self.disconnectdb()
        return (private_rows, public_rows)
    
    def checkkey(self, key, private=False):
        '''checkkey(key): check if the key exists in the db'''
        self.connectdb()
        if private==True:
            self.cursor.execute(self.selectSQL, (self.privatetable, 'key', key))
        else:
            self.curs.execute(self.selectSQL, (self.publictable, 'key', key))
        results = [row for row in self.curs]
        self.disconnectdb()
        if len(results) != 0:
            return True
        else:
            return False
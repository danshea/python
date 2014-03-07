#!/usr/bin/env python
'''
Name: GEO2sqlite.py
Date: 2014-03-07
Author: Dan Shea
Description:
    BioPython provides a very basic parser for Gene Omnibus Expression files.
    I have attempted to extend some additional functionality to it in order to
    read in GEO files and load them into a sqlite database.
    
    The GEO soft file format is as follows:
    Lines that begin with ^ identify the entity_type
    Valid types are:
    DATABASE
    PLATFORM
    SERIES
    SAMPLE
    
    entity_attribute lines begin with !
    
    For SAMPLE entity_types, we want to construct a TABLE for each and insert
    the table_rows.  col_defs.keys() or table_rows[0] will give is the column
    names for the table.
'''

from Bio import Geo
import sys
import gzip
import os
import os.path
import re
import sqlite3
import getopt

def parseGeo(filename):
    '''parseGeo(filename): open a GEO soft file and return the records it contains.'''
    try:
            fh = gzip.open(filename)
            records = Geo.parse(fh)
            return(records)
    except:
        try:
                fh = open(filename)
                records = Geo.parse(fh)
                return(records)
        except:
            print 'Could not open filename'
            sys.exit(3)

def connect_db(filename=':memory:'):
    '''connect_db: connect to sqlite db, default is in-memory'''
    try:
        connection = sqlite3.connect(filename)
        return(connection)
    except:
        print 'Could not connect to database'
        sys.exit(4)

def execute_query(connection, query, many=False, values=None):
    if many is False:
        try:
            with connection:
                connection.execute(query)
        except:
            print 'Caught exception during execution of:'
            print '{0:s}'.format(query)
            return 5
    else:
        try:
            with connection:
                connection.executemany(query, values)
        except:
            print 'Caught exception during execution of:'
            print '{0:s}'.format(query)
            return 5

def createTables(outfile, records):
    '''createTables(records): create sqlite db tables based on the records contained in the GEO file'''
    for record in records:
        if re.match('DATABASE', record.entity_type, re.IGNORECASE) != None:
            pass
        elif re.match('PLATFORM', record.entity_type, re.IGNORECASE) != None:
            pass
        elif re.match('SERIES', record.entity_type, re.IGNORECASE) != None:
            pass
        elif re.match('SAMPLE', record.entity_type, re.IGNORECASE) != None:
            # SAMPLE entity_types contain data to be inserted
            #
            # **********NOTE**********
            #
            # Usually your SQL operations will need to use values from Python variables.
            # You shouldn't assemble your query using Python's string operations because
            # doing so is insecure; it makes your program vulnerable to an SQL injection
            # attack (see http://xkcd.com/327/ for humorous example of what can go wrong).
            #
            # Instead, use the DB-API's parameter substitution. Put ? as a placeholder
            # wherever you want to use a value, and then provide a tuple of values as the
            # second argument to the cursor's execute() method.
            
            # create a parameter string based on the columns
            column_names = record.col_defs.keys()
            column_string = ','.join(['{{{0:d}:s}}'.format(i) for i in range(len(column_names))])
            # create a TABLE with the name entity_type
            # Since we are not allowed to parameterize table names, we must use string substitution here!
            create_table = 'CREATE TABLE IF NOT EXISTS {0:s} ({1:s})'.format(record.entity_id, column_string).format(*column_names)
            # commits if successful, rollback if failure
            execute_query(connect_db(outfile), create_table)
            print 'Created table: {0:s}'.format(record.entity_id)
            # create a parameterized string based on the number of values to insert
            value_string = ','.join('?'*len(record.col_defs.keys()))
            # Bulk INSERT values into the newly created table
            insert_many = 'INSERT INTO {0:s} ({1:s}) VALUES ({2:s})'.format(record.entity_id, column_string, value_string).format(*column_names)
            execute_query(connect_db(outfile), insert_many, many=True, values=record.table_rows[1:])
            print 'Inserted {0:d} rows into {1:s}'.format(len(record.table_rows[1:]), record.entity_id)            
    return 0

def usage(argv):
    print 'USAGE: {0:s}'.format(argv[0])
    print 'Options:'
    print '-h, --help\t\tDisplay this message.'
    print '-i, --infile\t\tThe input GEO soft file (can be gzip compressed.)'
    print '-o, --outfile\t\tThe output sqlite database file.'
    return 1

def main(argv=None):
    infile = None
    outfile = None
    if argv is None:
        argv = sys.argv
    try:
        optlist, args = getopt.getopt(argv[1:], 'hi:o:', ['help','infile=','outfile='])
    except:
        return(usage(argv))
    for o, a in optlist:
        if o in ('-h', '--help'):
            usage(argv)
        elif o in ('-i', '--infile'):
            infile = a
            if not os.path.isfile(a):
                print '{0:s} is not a file.'.format(a)
                return 3
        elif o in ('-o', '--outfile'):
            if os.path.isfile(a):
                print '{0:s} exists.'.format(a)
                print 'Exiting without modifying file.'
                return 4
            else:
                outfile = a
    if infile is None or outfile is None:
        usage(argv)        
    else:
        records = parseGeo(infile)
        return(createTables(outfile, records))

if __name__ == '__main__':
    sys.exit(main())
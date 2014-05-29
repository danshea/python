#!/usr/bin/env python

import ldap
import sys

def list_users(group, uri='ldap://wqldap-mdc1.med.harvard.edu', base='dc=med,dc=harvard,dc=edu',):
    try:
        lconn = ldap.initialize(uri)
        results = lconn.search_s(base, ldap.SCOPE_SUBTREE, '(&(objectclass=posixGroup)(cn={0:s}))'.format(group.lower()))
        for dn, result in results:
            print 'Members of {0:s}'.format(result['cn'][0])
            for member in result['memberUid']:
                print member
        return(0)
    except Exception as e:
        print e.mesg()
        return(1)

def main():
    if len(sys.argv) != 2:
        print 'usage: {0:s} ldapGroup'.format(sys.argv[0])
        return(1)
    else:
        return(list_users(sys.argv[1]))

if __name__ == '__main__':
    sys.exit(main())
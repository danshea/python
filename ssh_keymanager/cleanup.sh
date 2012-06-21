#!/bin/bash

files="ssh_keymanager.conf ssh_keymanager.log authorized_keys.db"
for file in $files; do
    [ -f ${file} ] && rm ${file} && echo "removed ${file}"
done
exit 0

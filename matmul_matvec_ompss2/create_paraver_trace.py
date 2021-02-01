#! /usr/bin/env python

# The TRACE.mpits file is often broken
# Need to list the files in set-0 in a TRACE.mpits file

import os

files = sorted(os.listdir('set-0'))

pwd = os.getcwd()

fp = open('TRACE.mpits', 'w')
for fname in files:
    if fname.endswith('.mpit'):
        print('%s/set-0/%s named' % (pwd, fname), file=fp)
fp.close()

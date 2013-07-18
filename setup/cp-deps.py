#!/usr/bin/env python
# copy dependencies

import sys, os, shutil

srcdir = sys.argv[1]
destdir = sys.argv[2]
files = sys.argv[3:]

# copy subset of python modules for packaging with arghmm

for fn in files:
    dirname = os.path.dirname(destdir + "/" + fn)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    fn1 = srcdir + "/" + fn
    fn2 = destdir + "/" + fn
    print "%s --> %s" % (fn1, fn2)
    shutil.copy(fn1, fn2)

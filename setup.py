#!/usr/bin/env python
#
# setup for the ArgWeaver library package
#
# use the following to install:
#   python setup.py install
#

from distutils.core import setup, Extension
import os

import arghmm
VERSION = arghmm.PROGRAM_VERSION_TEXT

scripts = [os.path.join('bin', x) for x in os.listdir('bin')]

setup(
    name='arghmm',
    version=VERSION,
    description='Ancestral recombination graph sampling method',
    long_description = """
            """,
    author='Matt Rasmussen',
    author_email='matt.rasmus@gmail.edu',

    packages=['arghmm', 'arghmm.deps.rasmus', 'arghmm.deps.compbio'],
    scripts=scripts,
)

#!/usr/bin/env python
#
# setup for the ARGweaver library package
#
# use the following to install:
#   python setup.py install
#

from distutils.core import setup, Extension
import os

import argweaver
VERSION = argweaver.PROGRAM_VERSION_TEXT

scripts = [os.path.join('bin', x) for x in os.listdir('bin')]

setup(
    name='argweaver',
    version=VERSION,
    description='Ancestral recombination graph sampling method',
    long_description = """
            """,
    author='Matt Rasmussen',
    author_email='matt.rasmus@gmail.edu',

    packages=['argweaver', 'argweaver.deps.rasmus', 'argweaver.deps.compbio'],
    scripts=scripts,
)

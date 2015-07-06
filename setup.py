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


def get_files(path, ext=''):
    """
    Get all files in a directory with a extension.
    """
    files = []
    for filename in os.listdir(path):
        if filename.endswith(ext):
            files.append(os.path.join(path, filename))
    return files


scripts = get_files('bin')
lib_src = get_files('src/argweaver', '.cpp')


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

    ext_modules=[
        Extension(
            'libargweaver',
            lib_src,
        )
    ],
)

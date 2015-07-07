#!/usr/bin/env python
#
# setup for the ARGweaver library package
#
# use the following to install:
#   python setup.py install
#

from distutils import log
from distutils.core import Command
from distutils.core import Extension
from distutils.core import setup
from distutils.command.build import build as _build
from distutils.command.build_clib import build_clib as _build_clib
from distutils.dep_util import newer_group
from distutils.sysconfig import customize_compiler
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


# Default to a C++ compiler, if possible.
if os.system('which g++') == 0:
    os.environ.setdefault('CC', 'g++')
    os.environ.setdefault('CXX', 'g++')


class build(_build):
    """
    Add more subcommands to the build command.
    """
    sub_commands = _build.sub_commands + [('build_prog', lambda cmd: True)]


class build_prog(_build_clib):
    """
    Build subcommand for compiling programs.

    This is implemented as subclass of build_clib since compiling a
    C/C++ program is very similar to compiling C/C++ library.
    """

    description = 'build compiled programs'

    # Config. Ideally, we could add this to setup().
    progs = [
        (
            'arg-sample',
            {
                'sources': lib_src + ['src/arg-sample.cpp'],
            }
        ),
        (
            'arg-summarize',
            {
                'sources': lib_src + ['src/arg-summarize.cpp'],
            }
        ),
        (
            'smc2bed',
            {
                'sources': lib_src + ['src/smc2bed.cpp'],
            }
        ),
    ]

    def initialize_options(self):
        # Call super class.
        _build_clib.initialize_options(self)

        # Directory for built binaries.
        self.build_bin = None

    def finalize_options(self):
        # Call super class.
        _build_clib.finalize_options(self)

        self.set_undefined_options('build', ('build_scripts', 'build_bin'))

        self.libraries = self.progs
        if self.libraries:
            self.check_library_list(self.libraries)

    def build_libraries(self, libraries):
        if not os.path.exists(self.build_bin):
            os.makedirs(self.build_bin)

        for (prog_name, build_info) in libraries:
            sources = build_info.get('sources')
            if sources is None or not isinstance(sources, (list, tuple)):
                raise DistutilsSetupError, \
                      ("in 'libraries' option ('%s'), " +
                       "'sources' must be present and must be " +
                       "a list of source filenames") % prog_name
            sources = list(sources)

            # Skip build, if program already built.
            prog_path = os.path.join(self.build_bin, prog_name)
            if not (self.force or newer_group(sources, prog_path, 'newer')):
                log.debug("skipping '%s' program (up-to-date)", prog_name)
                return

            log.info("building '%s' program", prog_name)

            macros = build_info.get('macros')
            include_dirs = build_info.get('include_dirs')
            objects = self.compiler.compile(sources,
                                            output_dir=self.build_temp,
                                            macros=macros,
                                            include_dirs=include_dirs,
                                            debug=self.debug)
            self.compiler.link_executable(objects, prog_name,
                                          output_dir=self.build_bin,
                                          debug=self.debug)


setup(
    name='argweaver',
    version=VERSION,
    description='Ancestral recombination graph sampling method',
    long_description = """
            """,
    author='Matt Rasmussen',
    author_email='matt.rasmus@gmail.edu',
    cmdclass={
        'build': build,
        'build_prog': build_prog,
    },
    packages=[
        'argweaver',
        'argweaver.deps.rasmus',
        'argweaver.deps.compbio',
    ],
    scripts=scripts,
    ext_modules=[
        Extension(
            'libargweaver',
            lib_src,
        )
    ],
)

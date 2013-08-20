"""

Tests for installing ArgWeaver.

"""

import os

import arghmm

from rasmus.testing import make_clean_dir


def run_cmd(cmd):
    assert os.system(cmd) == 0


def test_install():
    """
    Test installing ArgWeaver.
    """

    make_clean_dir("test/data/install")
    run_cmd("make install prefix=test/data/install > /dev/null")
    run_cmd("PYTHONPATH=test/data/install python -c 'import arghmm'")
    assert os.path.exists("test/data/install/bin/arg-sample")

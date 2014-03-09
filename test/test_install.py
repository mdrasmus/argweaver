"""

Tests for installing ARGweaver.

"""

import os

import argweaver

from rasmus.testing import make_clean_dir


def run_cmd(cmd):
    assert os.system(cmd) == 0


def test_install():
    """
    Test installing ARGweaver.
    """

    make_clean_dir("test/tmp/install")
    run_cmd("python setup.py clean > /dev/null")
    run_cmd("make install prefix=test/tmp/install > /dev/null")
    run_cmd("PYTHONPATH=test/tmp/install python -c 'import argweaver'")
    assert os.path.exists("test/tmp/install/bin/arg-sample")

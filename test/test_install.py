"""

Tests for installing ARGweaver.

"""

import os

import argweaver

from rasmus.testing import make_clean_dir


def run_cmd(cmd, ret_code=0):
    """
    Run shell command and assert success exit code.
    """
    assert os.system(cmd) == ret_code


def test_install_prog():
    """
    Test installing ARGweaver program.

    Use Makefile to install.
    """

    make_clean_dir("test/tmp/install")
    run_cmd("python setup.py clean > /dev/null")
    run_cmd("make install prefix=test/tmp/install > /dev/null")
    run_cmd("PYTHONPATH=test/tmp/install python -c 'import argweaver; "
            "assert argweaver.argweaverc.argweaverclib'")
    assert os.path.exists("test/tmp/install/bin/arg-sample")


def test_install_lib():
    """
    Test installing ARGweaver python lib.

    Use setup.py to install and check that module can be imported.
    Also ensure that c library is installed.
    """

    make_clean_dir("test/tmp/install")
    run_cmd("python setup.py clean > /dev/null")
    run_cmd(
        "python setup.py install --prefix=test/tmp/install "
        "--install-lib=test/tmp/install/lib/python/site-packages > /dev/null")
    run_cmd("cd test; PYTHONPATH=tmp/install/lib/python/site-packages "
            "python -c 'import argweaver; "
            "assert argweaver.argweaverc.argweaverclib'")

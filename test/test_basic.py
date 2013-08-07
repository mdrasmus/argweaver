"""

Tests for various HMM calculations

"""

import os

import arghmm

from rasmus.testing import make_clean_dir


def test_read_sites():
    """
    Test reading site files.
    """

    make_clean_dir("test/data/test_read_sites")
    os.system("""bin/arg-sim \
        -k 10 -L 10000 \
        -N 1e4 -r 1.5e-8 -m 2.5e-8 \
        --ntimes 10 --maxtime 400e3  \
        -o test/data/test_read_sites/0 > /dev/null""")

    sites = arghmm.arghmm_read_sites(
        'test/data/test_read_sites/0.sites', -1, -1)
    arghmm.arghmm_delete_sites(sites)

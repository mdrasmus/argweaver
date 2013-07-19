
import os

import arghmm

from rasmus.testing import make_clean_dir


def test_prog_small():

    make_clean_dir("test/data/test_prog_small")
    os.system("""bin/arg-sim \
        -k 4 -L 100000 \
        -N 1e4 -r 1.5e-8 -m 2.5e-8 \
        --ntimes 10 --maxtime 400e3  \
        -o test/data/test_prog_small/0 > /dev/null""")

    make_clean_dir("test/data/test_prog_small/0.sample")
    os.system("""bin/arg-sample -q \
-s test/data/test_prog_small/0.sites \
-x 1 -N 1e4 -r 1.5e-8 -m 2.5e-8 \
--ntimes 10 --maxtime 400e3 -c 20 \
-n 100 \
-o test/data/test_prog_small/0.sample/out""")

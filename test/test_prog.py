
import os
import subprocess

import argweaver

from compbio import arglib
from compbio import phylo
from rasmus import treelib
from rasmus.testing import make_clean_dir


EXIT_ERROR = 1


def run_cmd(cmd, retcode=0, set_pypath=True):
    """Run a command check its return code"""
    if set_pypath:
        cmd = "PYTHONPATH=. " + cmd
    print cmd
    result = subprocess.call(cmd, shell=True)
    assert result == retcode, (result, retcode)


def sites_split(names, col):
    part1 = []
    part2 = []

    c = col[0]
    for i in range(len(col)):
        if col[i] == c:
            part1.append(names[i])
        else:
            part2.append(names[i])
    return min([part1, part2], key=len)


def test_prog_small():
    """
    Test arg-sample on a small simulated dataset
    """

    make_clean_dir("test/tmp/test_prog_small")
    run_cmd("""bin/arg-sim \
        -k 4 -L 100000 \
        -N 1e4 -r 1.5e-8 -m 2.5e-8 \
        --ntimes 10 --maxtime 400e3  \
        -o test/tmp/test_prog_small/0 > /dev/null""")

    make_clean_dir("test/tmp/test_prog_small/0.sample")
    run_cmd("""bin/arg-sample -q \
        -s test/tmp/test_prog_small/0.sites \
        -x 1 -N 1e4 -r 1.5e-8 -m 2.5e-8 \
        --ntimes 10 --maxtime 400e3 -c 20 \
        -n 10 --sample-step 1 \
        -o test/tmp/test_prog_small/0.sample/out""")


def test_prog_popsizes():
    """
    Ensure arg-sample can use multiple population sizes from a file.
    """

    make_clean_dir("test/tmp/test_prog_popsizes")

    # Setup popsize config file.
    popsizes = range(10000, 0, -1000)
    with open('test/tmp/test_prog_popsizes/popsizes', 'w') as out:
        out.write('\n'.join(map(str, popsizes)))

    run_cmd("""bin/arg-sim \
        -k 4 -L 100000 \
        -N 1e4 -r 1.5e-8 -m 2.5e-8 \
        --ntimes 10 --maxtime 400e3  \
        -o test/tmp/test_prog_popsizes/0 > /dev/null""")

    make_clean_dir("test/tmp/test_prog_popsizes/0.sample")
    run_cmd("""bin/arg-sample -q \
        -s test/tmp/test_prog_popsizes/0.sites \
        -x 1 -N 1e4 -r 1.5e-8 -m 2.5e-8 \
        --ntimes 10 --maxtime 400e3 -c 20 \
        --popsize-file test/tmp/test_prog_popsizes/popsizes \
        -n 10 --sample-step 1 \
        -o test/tmp/test_prog_popsizes/0.sample/out""")

    # Write too few population sizes to config file.
    popsizes_bad = range(10000, 5000, -1000)
    with open('test/tmp/test_prog_popsizes/popsizes-bad', 'w') as out:
        out.write('\n'.join(map(str, popsizes_bad)))

    # Expect an error.
    run_cmd("""bin/arg-sample -q \
        -s test/tmp/test_prog_popsizes/0.sites \
        -x 1 -N 1e4 -r 1.5e-8 -m 2.5e-8 \
        --ntimes 10 --maxtime 400e3 -c 20 \
        --popsize-file test/tmp/test_prog_popsizes/popsizes-bad \
        -n 10 --sample-step 1 \
        -o test/tmp/test_prog_popsizes/0.sample/out-bad""",
            EXIT_ERROR)


def test_prog_resume():
    """
    Ensure that a run can be resumed
    """

    test_prog_small()

    run_cmd("""bin/arg-sample -q \
        -s test/tmp/test_prog_small/0.sites \
        -x 1 -N 1e4 -r 1.5e-8 -m 2.5e-8 \
        --ntimes 20 --maxtime 400e3 -c 20 \
        -n 20 --sample-step 1 --resume \
        -o test/tmp/test_prog_small/0.sample/out""")


def test_tmrca():

    if not os.path.exists("test/tmp/test_prog_small/0.sample"):
        test_prog_small()

    run_cmd("""bin/arg-extract-tmrca \
        test/tmp/test_prog_small/0.sample/out.%d.smc.gz \
        > test/tmp/test_prog_small/0.tmrca.txt""")


def test_popsize():

    if not os.path.exists("test/tmp/test_prog_small/0.sample"):
        test_prog_small()

    run_cmd("""bin/arg-extract-popsize \
        test/tmp/test_prog_small/0.sample/out.%d.smc.gz \
        > test/tmp/test_prog_small/0.popsize.txt""")


def test_breaks():

    if not os.path.exists("test/tmp/test_prog_small/0.sample"):
        test_prog_small()

    run_cmd("""bin/arg-extract-breaks \
        test/tmp/test_prog_small/0.sample/out.%d.smc.gz \
        > test/tmp/test_prog_small/0.breaks.txt""")


def test_recomb():

    if not os.path.exists("test/tmp/test_prog_small/0.sample"):
        test_prog_small()

    run_cmd("""bin/arg-extract-recomb \
        test/tmp/test_prog_small/0.sample/out.%d.smc.gz \
        > test/tmp/test_prog_small/0.recomb.txt""")


def test_treelen():

    if not os.path.exists("test/tmp/test_prog_small/0.sample"):
        test_prog_small()

    run_cmd("""bin/arg-extract-treelen \
        test/tmp/test_prog_small/0.sample/out.%d.smc.gz \
        > test/tmp/test_prog_small/0.treelen.txt""")


def test_ages():

    if not os.path.exists("test/tmp/test_prog_small/0.sample"):
        test_prog_small()

    run_cmd("""bin/arg-extract-ages \
        test/tmp/test_prog_small/0.sample/out.%d.smc.gz \
        test/tmp/test_prog_small/0.sites \
        > test/tmp/test_prog_small/0.ages.txt""")


def test_arg_smc():

    if not os.path.exists("test/tmp/test_prog_small/0.sample"):
        test_prog_small()

    run_cmd("""bin/arg2smc \
        test/tmp/test_prog_small/0.arg \
        test/tmp/test_prog_small/0.smc""")
    run_cmd("""bin/smc2arg \
        test/tmp/test_prog_small/0.smc \
        test/tmp/test_prog_small/0.smc.arg""")


def test_db():

    if not os.path.exists("test/tmp/test_prog_small/0.sample"):
        test_prog_small()

    run_cmd("""bin/arg-db \
        --db test/tmp/test_prog_small/arg.db \
        --add-smc test/tmp/test_prog_small/0.sample/out.0.smc.gz
        """)
    run_cmd("""bin/sites-db \
        --db test/tmp/test_prog_small/sites.db \
        --add-sites test/tmp/test_prog_small/0.sites
        """)


def _test_prog_infsites():

    make_clean_dir("test/tmp/test_prog_infsites")

    run_cmd("""bin/arg-sim \
        -k 40 -L 200000 \
        -N 1e4 -r 1.5e-8 -m 2.5e-8 --infsites \
        --ntimes 20 --maxtime 400e3 \
        -o test/tmp/test_prog_infsites/0""")

    make_clean_dir("test/tmp/test_prog_infsites/0.sample")
    run_cmd("""bin/arg-sample \
        -s test/tmp/test_prog_infsites/0.sites \
        -N 1e4 -r 1.5e-8 -m 2.5e-8 \
        --ntimes 5 --maxtime 100e3 -c 1 \
        --climb 0 -n 20 --infsites \
        -x 1 \
        -o test/tmp/test_prog_infsites/0.sample/out""")

    arg = argweaver.read_arg(
        "test/tmp/test_prog_infsites/0.sample/out.0.smc.gz")
    sites = argweaver.read_sites("test/tmp/test_prog_infsites/0.sites")
    print "names", sites.names
    print

    noncompats = []
    for block, tree in arglib.iter_local_trees(arg):
        tree = tree.get_tree()
        treelib.remove_single_children(tree)
        phylo.hash_order_tree(tree)
        for pos, col in sites.iter_region(block[0]+1, block[1]+1):
            assert block[0]+1 <= pos <= block[1]+1, (block, pos)
            split = sites_split(sites.names, col)
            node = arglib.split_to_tree_branch(tree, split)
            if node is None:
                noncompats.append(pos)
                print "noncompat", block, pos, col
                print phylo.hash_tree(tree)
                print tree.leaf_names()
                print "".join(col[sites.names.index(name)]
                              for name in tree.leaf_names())
                print split
                print
    print "num noncompats", len(noncompats)

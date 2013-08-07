
import os

import arghmm

from compbio import arglib
from compbio import phylo
from rasmus import treelib
from rasmus.testing import make_clean_dir


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
        -n 10 --sample-step 1 \
        -o test/data/test_prog_small/0.sample/out""")


def test_prog_resume():
    """
    Ensure that a run can be resumed
    """

    test_prog_small()

    os.system("""bin/arg-sample -q \
        -s test/data/test_prog_small/0.sites \
        -x 1 -N 1e4 -r 1.5e-8 -m 2.5e-8 \
        --ntimes 20 --maxtime 400e3 -c 20 \
        -n 20 --sample-step 1 --resume \
        -o test/data/test_prog_small/0.sample/out""")


def _test_prog_infsites():

    make_clean_dir("test/data/test_prog_infsites")

    os.system("""arg-sim \
        -k 40 -L 200000 \
        -N 1e4 -r 1.5e-8 -m 2.5e-8 --infsites \
        --ntimes 20 --maxtime 400e3 \
        -o test/data/test_prog_infsites/0""")

    make_clean_dir("test/data/test_prog_infsites/0.sample")
    os.system("""arg-sample \
        -s test/data/test_prog_infsites/0.sites \
        -N 1e4 -r 1.5e-8 -m 2.5e-8 \
        --ntimes 5 --maxtime 100e3 -c 1 \
        --climb 0 -n 20 --infsites \
        -x 1 \
        -o test/data/test_prog_infsites/0.sample/out""")

    arg = arghmm.read_arg("test/data/test_prog_infsites/0.sample/out.0.smc.gz")
    sites = arghmm.read_sites("test/data/test_prog_infsites/0.sites")
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

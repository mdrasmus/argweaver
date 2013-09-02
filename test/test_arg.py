
from itertools import izip

import arghmm
from arghmm import arghmmc

from compbio import arglib
from compbio import phylo


def arg_equal(arg, arg2):

    # test recomb points
    recombs = sorted(x.pos for x in arg if x.event == "recomb")
    recombs2 = sorted(x.pos for x in arg2 if x.event == "recomb")
    assert recombs == recombs2

    # check local tree topologies
    for (start, end), tree in arglib.iter_tree_tracks(arg):
        pos = (start + end) / 2.0

        arglib.remove_single_lineages(tree)
        tree1 = tree.get_tree()

        tree2 = arg2.get_marginal_tree(pos)
        arglib.remove_single_lineages(tree2)
        tree2 = tree2.get_tree()

        hash1 = phylo.hash_tree(tree1)
        hash2 = phylo.hash_tree(tree2)
        print
        print pos
        print hash1
        print hash2
        assert hash1 == hash2

    # check sprs
    sprs1 = arglib.iter_arg_sprs(arg, use_leaves=True)
    sprs2 = arglib.iter_arg_sprs(arg2, use_leaves=True)

    for (pos1, recomb1, coal1), (pos2, recomb2, coal2) in izip(sprs1, sprs2):
        recomb1 = (sorted(recomb1[0]), recomb1[1])
        recomb2 = (sorted(recomb2[0]), recomb2[1])
        coal1 = (sorted(coal1[0]), coal1[1])
        coal2 = (sorted(coal2[0]), coal2[1])

        print
        print (pos1, recomb1, coal1)
        print (pos2, recomb2, coal2)

        # check pos, leaves, time
        assert pos1 == pos2
        assert recomb1 == recomb2
        assert coal1 == coal2


def test_arg_convert():
    """
    Test conversion for python to C args
    """

    k = 10
    n = 1e4
    rho = 1.5e-8 * 20
    length = 10000
    times = arghmm.get_time_points(ntimes=20, maxtime=200000)

    arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                 times=times)

    # convert to C++ and back
    trees, names = arghmmc.arg2ctrees(arg, times)
    arg2 = arghmmc.ctrees2arg(trees, names, times)

    arg_equal(arg, arg2)

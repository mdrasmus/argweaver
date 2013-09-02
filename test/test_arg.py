
from itertools import izip
import StringIO

import nose

import arghmm
from arghmm import arghmmc

from compbio import arglib
from compbio import phylo


def arg_equal(arg, arg2):

    # test recomb points
    recombs = sorted(x.pos for x in arg if x.event == "recomb")
    recombs2 = sorted(x.pos for x in arg2 if x.event == "recomb")
    nose.tools.assert_equal(recombs, recombs2)

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
        nose.tools.assert_equal(hash1, hash2)

    # check sprs
    sprs1 = arglib.iter_arg_sprs(arg, use_leaves=True)
    sprs2 = arglib.iter_arg_sprs(arg2, use_leaves=True)

    for (pos1, recomb1, coal1), (pos2, recomb2, coal2) in izip(sprs1, sprs2):
        recomb1 = (sorted(recomb1[0]), recomb1[1])
        recomb2 = (sorted(recomb2[0]), recomb2[1])
        coal1 = (sorted(coal1[0]), coal1[1])
        coal2 = (sorted(coal2[0]), coal2[1])

        # check pos, leaves, time
        nose.tools.assert_equal(pos1, pos2)
        nose.tools.assert_equal(recomb1, recomb2)
        nose.tools.assert_equal(coal1, coal2)


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


def test_node_numbering():
    """
    Test node numbering across ARG.

    A node should keep the same numbering until it is broken by
    recombination. The new recoal node should take the index of the broken
    node.
    """

    k = 10
    n = 1e4
    rho = 1.5e-8 * 20
    length = 10000
    times = arghmm.get_time_points(ntimes=20, maxtime=200000)

    arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                 times=times)

    (ptrees, ages, sprs, blocks), all_nodes = (
        arghmmc.get_treeset(arg, times))

    # check nodes list
    nnodes = len(all_nodes[0])
    last_nodes = None
    for i, nodes in enumerate(all_nodes):
        if last_nodes:
            recombj = sprs[i][0]
            brokenj = ptrees[i][recombj]
            for j in range(nnodes):
                if j != brokenj:
                    nose.tools.assert_equal(last_nodes[j], nodes[j])

        last_nodes = nodes


_smc_example = "NAMES	n0	n2	n3	n1\nREGION	chr	1	100000\nTREE	1	60200	((2:531.657854[&&NHX:age=0.000000],(3:151.328043[&&NHX:age=0.000000],1:151.328043[&&NHX:age=0.000000])6:380.329810[&&NHX:age=151.328043])4:0.000000[&&NHX:age=531.657854],0:531.657854[&&NHX:age=0.000000])5[&&NHX:age=531.657854];\nSPR	60200	2	151.328043	5	3889.916442\nTREE	60201	68400	(2:3889.916442[&&NHX:age=0.000000],((3:151.328043[&&NHX:age=0.000000],1:151.328043[&&NHX:age=0.000000])6:380.329810[&&NHX:age=151.328043],0:531.657854[&&NHX:age=0.000000])5:3358.258588[&&NHX:age=531.657854])4[&&NHX:age=3889.916442];\nSPR	68400	2	1487.533324	5	3889.916442\nTREE	68401	81840	(2:3889.916442[&&NHX:age=0.000000],((3:151.328043[&&NHX:age=0.000000],1:151.328043[&&NHX:age=0.000000])6:380.329810[&&NHX:age=151.328043],0:531.657854[&&NHX:age=0.000000])5:3358.258588[&&NHX:age=531.657854])4[&&NHX:age=3889.916442];\nSPR	81840	2	1487.533324	5	3889.916442\nTREE	81841	88580	(2:3889.916442[&&NHX:age=0.000000],((3:151.328043[&&NHX:age=0.000000],1:151.328043[&&NHX:age=0.000000])6:380.329810[&&NHX:age=151.328043],0:531.657854[&&NHX:age=0.000000])5:3358.258588[&&NHX:age=531.657854])4[&&NHX:age=3889.916442];\nSPR	88580	3	151.328043	4	9927.778924\nTREE	88581	100000	(3:9927.778924[&&NHX:age=0.000000],(2:3889.916442[&&NHX:age=0.000000],(1:531.657854[&&NHX:age=0.000000],0:531.657854[&&NHX:age=0.000000])5:3358.258588[&&NHX:age=531.657854])4:6037.862482[&&NHX:age=3889.916442])6[&&NHX:age=9927.778924];"  # nopep8


def test_arg2smc():
    """
    Test that an ARG be converted to SMC and back.
    """
    infile = StringIO.StringIO(_smc_example)
    smc = arghmm.read_smc(infile)
    arg = arghmm.smc2arg(smc)

    smc2 = list(arghmm.arg2smc(arg, names=smc[0]["names"]))
    arg2 = arghmm.smc2arg(smc2)

    arg_equal(arg, arg2)

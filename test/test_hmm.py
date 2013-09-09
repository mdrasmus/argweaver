"""

Tests for various HMM calculations

"""

from itertools import izip
from math import exp
from math import log

import argweaver
import argweaver.popsize
from argweaver import argweaverc

from compbio import arglib
from rasmus import util
from rasmus.testing import fequal


def test_trans_two():
    """
    Calculate transition probabilities for k=2

    Only calculate a single matrix
    """

    k = 2
    n = 1e4
    rho = 1.5e-8 * 20
    length = 1000
    times = argweaver.get_time_points(ntimes=5, maxtime=200000)
    time_steps = [times[i] - times[i-1]
                  for i in range(1, len(times))]
    time_steps.append(200000*10000.0)
    popsizes = [n] * len(times)

    arg = arglib.sample_arg(k, 2*n, rho, start=0, end=length)

    argweaver.discretize_arg(arg, times)
    print "recomb", arglib.get_recomb_pos(arg)

    arg = argweaver.make_trunk_arg(0, length, "n0")

    pos = 10
    tree = arg.get_marginal_tree(pos)
    nlineages = argweaver.get_nlineages_recomb_coal(tree, times)
    states = list(argweaver.iter_coal_states(tree, times))
    mat = argweaver.calc_transition_probs(
        tree, states, nlineages,
        times, time_steps, popsizes, rho)

    nstates = len(states)

    def coal(j):
        return 1.0 - exp(-time_steps[j]/(2.0 * n))

    def recoal2(k, j):
        p = coal(j)
        for m in range(k, j):
            p *= 1.0 - coal(m)
        return p

    def recoal(k, j):
        if j == nstates-1:
            return exp(- sum(time_steps[m] / (2.0 * n)
                             for m in range(k, j)))
        else:
            return ((1.0 - exp(-time_steps[j]/(2.0 * n))) *
                    exp(- sum(time_steps[m] / (2.0 * n)
                              for m in range(k, j))))

    def isrecomb(i):
        return 1.0 - exp(-max(rho * 2.0 * times[i], rho))

    def recomb(i, k):
        treelen = 2*times[i] + time_steps[i]
        if k < i:
            return 2.0 * time_steps[k] / treelen / 2.0
        else:
            return time_steps[k] / treelen / 2.0

    def trans(i, j):
        a = states[i][1]
        b = states[j][1]

        p = sum(recoal(k, b) * recomb(a, k)
                for k in range(0, min(a, b)+1))
        p += sum(recoal(k, b) * recomb(a, k)
                 for k in range(0, min(a, b)+1))
        p *= isrecomb(a)
        if i == j:
            p += 1.0 - isrecomb(a)
        return p

    for i in range(len(states)):
        for j in range(len(states)):
            print isrecomb(states[i][1])
            print states[i], states[j], mat[i][j], log(trans(i, j))
            fequal(mat[i][j], log(trans(i, j)))

        # recombs add up to 1
        fequal(sum(recomb(i, k) for k in range(i+1)), 0.5)

        # recoal add up to 1
        fequal(sum(recoal(i, j) for j in range(i, nstates)), 1.0)

        # recomb * recoal add up to .5
        fequal(sum(sum(recoal(k, j) * recomb(i, k)
                       for k in range(0, min(i, j)+1))
                   for j in range(0, nstates)), 0.5)

        fequal(sum(trans(i, j) for j in range(len(states))), 1.0)


def test_trans():
    """
    Calculate transition probabilities
    """

    k = 4
    n = 1e4
    rho = 1.5e-8 * 20
    length = 1000
    times = argweaver.get_time_points(ntimes=4, maxtime=200000)
    popsizes = [n] * len(times)

    arg = arglib.sample_arg(k, 2*n, rho, start=0, end=length)
    argweaver.discretize_arg(arg, times)

    pos = 10
    tree = arg.get_marginal_tree(pos)

    assert argweaverc.assert_transition_probs(tree, times, popsizes, rho)


def test_trans_switch():
    """
    Calculate transition probabilities for switch matrix

    Only calculate a single matrix
    """

    k = 12
    n = 1e4
    rho = 1.5e-8 * 20
    length = 1000
    times = argweaver.get_time_points(ntimes=20, maxtime=200000)
    popsizes = [n] * len(times)

    recombs = []

    while len(recombs) == 0:
        arg = argweaver.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                        times=times)
        recombs = [x.pos for x in arg if x.event == "recomb"]

    pos = recombs[0]
    tree = arg.get_marginal_tree(pos-.5)
    rpos, r, c = arglib.iter_arg_sprs(arg, start=pos-.5).next()
    spr = (r, c)

    assert argweaverc.assert_transition_switch_probs(
        tree, spr, times, popsizes, rho)


def test_trans_internal():
    """
    Calculate transition probabilities for internal branch re-sampling

    Only calculate a single matrix
    """

    k = 5
    n = 1e4
    rho = 1.5e-8 * 20
    length = 1000
    times = argweaver.get_time_points(ntimes=5, maxtime=200000)
    popsizes = [n] * len(times)

    arg = arglib.sample_arg(k, 2*n, rho, start=0, end=length)
    argweaver.discretize_arg(arg, times)

    pos = 10
    tree = arg.get_marginal_tree(pos)

    assert argweaverc.assert_transition_probs_internal(
        tree, times, popsizes, rho)


def test_trans_switch_internal():
    """
    Calculate transition probabilities for switch matrix and internal branches

    Only calculate a single matrix
    """

    k = 10
    n = 1e4
    rho = 1.5e-8 * 20
    length = int(100e3) / 20
    times = argweaver.get_time_points(ntimes=20, maxtime=200000)
    popsizes = [n] * len(times)

    arg = argweaver.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                    times=times)
    trees, names = argweaverc.arg2ctrees(arg, times)

    assert argweaverc.assert_transition_probs_switch_internal(
        trees, times, popsizes, rho)


def test_emit():
    """
    Calculate emission probabilities
    """

    k = 10
    n = 1e4
    rho = 1.5e-8 * 20
    mu = 2.5e-8 * 20
    length = int(1e3) / 20
    times = argweaver.get_time_points(ntimes=20, maxtime=200000)

    arg = argweaver.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                    times=times)

    muts = argweaver.sample_arg_mutations(arg, mu, times)
    seqs = argweaver.make_alignment(arg, muts)

    new_name = "n%d" % (k-1)
    arg = argweaver.remove_arg_thread(arg, new_name)

    trees, names = argweaverc.arg2ctrees(arg, times)
    seqs2, nseqs, seqlen = argweaverc.seqs2cseqs(seqs, names + [new_name])

    assert argweaverc.argweaver_assert_emit(trees, len(times), times, mu,
                                            seqs2, nseqs, seqlen)


def test_emit_internal():
    """
    Calculate emission probabilities for internal branches
    """

    k = 10
    n = 1e4
    rho = 1.5e-8 * 20
    mu = 2.5e-8 * 20
    length = int(10e3) / 20
    times = argweaver.get_time_points(ntimes=20, maxtime=200000)

    arg = argweaver.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                    times=times)

    muts = argweaver.sample_arg_mutations(arg, mu, times)
    seqs = argweaver.make_alignment(arg, muts)

    trees, names = argweaverc.arg2ctrees(arg, times)
    seqs2, nseqs, seqlen = argweaverc.seqs2cseqs(seqs, names)

    assert argweaverc.argweaver_assert_emit_internal(
        trees, len(times), times, mu, seqs2, nseqs, seqlen)


def test_prior_counts():
    """
    Calculate initial tree prior
    """

    a = 10
    n = 1e4
    t = 1e3

    for b in range(1, a):
        x = argweaverc.prob_coal_counts_matrix(a, b, t, 2*n)
        y = argweaver.popsize.prob_coal_counts(a, b, t, 2*n)
        fequal(x, y)


def test_forward():

    k = 4
    n = 1e4
    rho = 1.5e-8 * 20
    mu = 2.5e-8 * 20
    length = int(100e3 / 20)
    times = argweaver.get_time_points(ntimes=100)

    arg = arglib.sample_arg_smc(k, 2*n, rho, start=0, end=length)
    muts = arglib.sample_arg_mutations(arg, mu)
    seqs = arglib.make_alignment(arg, muts)

    print "muts", len(muts)
    print "recomb", len(arglib.get_recomb_pos(arg))

    argweaver.discretize_arg(arg, times)

    # remove chrom
    new_name = "n%d" % (k - 1)
    arg = argweaver.remove_arg_thread(arg, new_name)

    carg = argweaverc.arg2ctrees(arg, times)

    util.tic("C fast")
    probs1 = argweaverc.argweaver_forward_algorithm(carg, seqs, times=times)
    util.toc()

    util.tic("C slow")
    probs2 = argweaverc.argweaver_forward_algorithm(carg, seqs, times=times,
                                                    slow=True)
    util.toc()

    for i, (col1, col2) in enumerate(izip(probs1, probs2)):
        for a, b in izip(col1, col2):
            fequal(a, b, rel=.0001)


def test_arg_joint():
    """
    Compute joint probability of an ARG
    """

    k = 2
    n = 1e4
    rho = 1.5e-8 * 20
    mu = 2.5e-8 * 20
    length = 10000
    times = argweaver.get_time_points(ntimes=20, maxtime=200000)

    arg = argweaver.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                    times=times)
    muts = argweaver.sample_arg_mutations(arg, mu, times=times)
    seqs = arglib.make_alignment(arg, muts)

    lk = argweaver.calc_joint_prob(arg, seqs, mu=mu, rho=rho, times=times)
    print lk

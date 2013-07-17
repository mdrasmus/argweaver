"""

Tests for various HMM calculations

"""

import unittest, random

from rasmus.common import *
from rasmus.testing import *
from compbio import arglib

import arghmm
import arghmm.popsize


class TestHMM (unittest.TestCase):

    def test_trans_two(self):
        """
        Calculate transition probabilities for k=2

        Only calculate a single matrix
        """

        k = 2
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 1000
        times = arghmm.get_time_points(ntimes=5, maxtime=200000)

        arg = arglib.sample_arg(k, 2*n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        arghmm.discretize_arg(arg, times)
        print "recomb", arglib.get_recomb_pos(arg)

        new_name = "n%d" % (k-1)
        arg = arghmm.make_trunk_arg(0, length, "n0")
        model = arghmm.ArgHmm(arg, seqs, new_name=new_name,
                              popsize=n, rho=rho, mu=mu,
                              times=times)

        pos = 10
        tree = arg.get_marginal_tree(pos)
        model.check_local_tree(pos, force=True)
        mat = arghmm.calc_transition_probs(
            tree, model.states[pos], model.nlineages,
            model.times, model.time_steps, model.popsizes, rho)

        states = model.states[pos]
        nstates = len(states)

        def coal(j):
            return 1.0 - exp(-model.time_steps[j]/(2.0 * n))

        def recoal2(k, j):
            p = coal(j)
            for m in range(k, j):
                p *= 1.0 - coal(m)
            return p

        def recoal(k, j):
            if j == nstates-1:
                return exp(- sum(model.time_steps[m] / (2.0 * n)
                              for m in range(k, j)))
            else:
                return ((1.0 - exp(-model.time_steps[j]/(2.0 * n))) * 
                        exp(- sum(model.time_steps[m] / (2.0 * n)
                                  for m in range(k, j))))

        def isrecomb(i):
            return 1.0 - exp(-max(rho * 2.0 * model.times[i], rho))

        def recomb(i, k):
            treelen = 2*model.times[i] + model.time_steps[i]
            if k < i:
                return 2.0 * model.time_steps[k] / treelen / 2.0
            else:
                return model.time_steps[k] / treelen / 2.0

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


    def test_trans(self):
        """
        Calculate transition probabilities
        """

        k = 4
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 1000
        times = arghmm.get_time_points(ntimes=4, maxtime=200000)
        popsizes = [n] * len(times)

        arg = arglib.sample_arg(k, 2*n, rho, start=0, end=length)
        arghmm.discretize_arg(arg, times)

        pos = 10
        tree = arg.get_marginal_tree(pos)

        assert arghmm.assert_transition_probs(tree, times, popsizes, rho)


    def test_trans_switch(self):
        """
        Calculate transition probabilities for k=2

        Only calculate a single matrix
        """

        k = 12
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 1000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        popsizes = [n] * len(times)

        recombs = []

        while len(recombs) == 0:
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            recombs = [x.pos for x in arg if x.event == "recomb"]

        pos = recombs[0]
        tree = arg.get_marginal_tree(pos-.5)
        rpos, r, c = arglib.iter_arg_sprs(arg, start=pos-.5).next()
        spr = (r, c)

        assert arghmm.assert_transition_switch_probs(tree, spr,
                                                     times, popsizes, rho)


    def test_trans_internal(self):
        """
        Calculate transition probabilities for k=2

        Only calculate a single matrix
        """

        k = 5
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 1000
        times = arghmm.get_time_points(ntimes=5, maxtime=200000)
        popsizes = [n] * len(times)

        arg = arglib.sample_arg(k, 2*n, rho, start=0, end=length)
        arghmm.discretize_arg(arg, times)

        pos = 10
        tree = arg.get_marginal_tree(pos)

        assert arghmm.assert_transition_probs_internal(
            tree, times, popsizes, rho)


    def test_trans_switch_internal(self):
        """
        Calculate transition probabilities for k=2

        Only calculate a single matrix
        """

        k = 10
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(100e3) / 20
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        popsizes = [n] * len(times)

        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        trees, names = arghmm.arg2ctrees(arg, times)

        assert arghmm.assert_transition_probs_switch_internal(
            trees, times, popsizes, rho)


    def test_emit(self):
        """
        Calculate emission probabilities
        """

        k = 10
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(1e3) / 20
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)

        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)

        muts = arghmm.sample_arg_mutations(arg, mu, times)
        seqs = arghmm.make_alignment(arg, muts)

        new_name = "n%d" % (k-1)
        arg = arghmm.remove_arg_thread(arg, new_name)

        trees, names = arghmm.arg2ctrees(arg, times)
        seqs2, nseqs, seqlen = arghmm.seqs2cseqs(seqs, names + [new_name])

        assert arghmm.arghmm_assert_emit(trees, len(times), times, mu,
                                         seqs2, nseqs, seqlen)

    def test_emit_internal(self):
        """
        Calculate emission probabilities
        """

        k = 10
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(10e3) / 20
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)

        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)

        muts = arghmm.sample_arg_mutations(arg, mu, times)
        seqs = arghmm.make_alignment(arg, muts)

        trees, names = arghmm.arg2ctrees(arg, times)
        seqs2, nseqs, seqlen = arghmm.seqs2cseqs(seqs, names)

        assert arghmm.arghmm_assert_emit_internal(trees, len(times), times, mu,
                                                  seqs2, nseqs, seqlen)



    def test_prior_counts(self):
        """
        Calculate initial tree prior
        """

        a = 10
        n = 1e4
        t = 1e3

        for b in range(1, a):
            x = arghmm.prob_coal_counts_matrix(a, b, t, 2*n)
            y = arghmm.popsize.prob_coal_counts(a, b, t, 2*n)
            fequal(x, y)



    def test_forward_c(self):

        k = 4
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(100e3 / 20)
        times = arghmm.get_time_points(ntimes=100)

        arg = arglib.sample_arg_smc(k, 2*n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        print "muts", len(muts)
        print "recomb", len(arglib.get_recomb_pos(arg))

        arghmm.discretize_arg(arg, times)

        # remove chrom
        new_name = "n%d" % (k - 1)
        arg = arghmm.remove_arg_thread(arg, new_name)

        carg = arghmm.arg2ctrees(arg, times)

        util.tic("C fast")
        probs1 = arghmm.arghmm_forward_algorithm(carg, seqs, times=times)
        util.toc()

        util.tic("C slow")
        probs2 = arghmm.arghmm_forward_algorithm(carg, seqs, times=times,
                                                 slow=True)
        util.toc()


        for i, (col1, col2) in enumerate(izip(probs1, probs2)):
            for a, b in izip(col1, col2):
                try:
                    fequal(a, b, rel=.0001)
                except:
                    print model.states[i]
                    print i, col1
                    print i, col2
                    raise



    def test_arg_joint(self):
        """
        Compute joint probability of an ARG
        """

        k = 2
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        refine = 0

        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg, mu, times=times)
        seqs = arglib.make_alignment(arg, muts)

        lk = arghmm.calc_joint_prob(arg, seqs, mu=mu, rho=rho, times=times)
        print lk



#=============================================================================
if __name__ == "__main__":

    test_main()



import unittest, random

import arghmm

from rasmus.common import *
from rasmus import stats, hmm
from rasmus.testing import *

from compbio import coal, arglib, fasta


def test_arg_equal(arg, arg2):

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



class Basic (unittest.TestCase):

    def test_remove_thread(self):
        """
        Remove a leaf from an ARG
        """

        k = 3
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8
        length = 1000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        
        times = arghmm.get_time_points(10)
        arghmm.discretize_arg(arg, times)
        chrom = "n%d" % (k-1)
        arg = arghmm.remove_arg_thread(arg, chrom)
        
        recomb = arglib.get_recomb_pos(arg)
        print "recomb", recomb

        tree = arg.get_marginal_tree(-.5)

        draw_tree_names(tree.get_tree(), minlen=5, maxlen=5)        
        print sorted([(x.pos, x.event) for x in tree])

        


    def test_recomb(self):
        """
        Investigate the fact that some recombinations are not visible
        """

        k = 3
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8
        length = 1000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        
        times = arghmm.get_time_points(10)
        arghmm.discretize_arg(arg, times)
        arg.set_ancestral()
        arg.prune()

        recombs = arglib.get_recomb_pos(arg)

        # find recombs by walking
        recombs2 = []
        i = 0
        while True:
            tree = arg.get_marginal_tree(i-.5)
            recomb = arghmm.find_tree_next_recomb(tree, i+1, tree=True)
            if recomb:
                recombs2.append(recomb.pos)
                i = recomb.pos
            else:
                break

        # these are suppose to differ because some recombination occur
        # in the hole of ancestral sequence intervals
        print recombs
        print recombs2

        arglib.write_arg("tmp/b.arg", arg)


    def test_nlineages(self):
        """
        Test lineage counting
        """
        k = 4
        n = 1e4
        rho = 1.5e-8 * 1
        length = 1000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)

        times = arghmm.get_time_points(ntimes=6)
        arghmm.discretize_arg(arg, times)
        tree = arg.get_marginal_tree(0)

        nlineages, nrecombs, ncoals = arghmm.get_nlineages_recomb_coal(
            tree, times)
        
        treelib.draw_tree_names(tree.get_tree(), scale=4e-3)

        print list(arghmm.iter_coal_states(tree, times))
        print nlineages
        self.assert_(nlineages == sorted(nlineages, reverse=True))

        print nlineages
        print nrecombs
        print ncoals


    def test_states(self):
        """
        Test state enumeration
        """
        k = 2
        n = 1e4
        rho = 1.5e-8 * 100
        length = 1000

        for i in xrange(20):
            arg = arglib.sample_arg(k, n, rho, start=0, end=length)

            times = arghmm.get_time_points(10)
            arghmm.discretize_arg(arg, times)
            tree = arg.get_marginal_tree(0)

            states = list(arghmm.iter_coal_states(tree, times))

            treelib.draw_tree_names(tree.get_tree(), scale=4e-4,
                                    minlen=6, maxlen=6)
            print states


    def test_pars_seq(self):
        """
        Test parsimony ancestral sequence inference
        """

        k = 10
        n = 1e4
        rho = 1.5e-8
        mu = 2.5e-8 * 100
        length = 1000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        pos = int(muts[0][2])
        tree = arg.get_marginal_tree(pos)

        print "pos =", pos
        treelib.draw_tree_names(tree.get_tree(), scale=4e-4, minlen=5)

        arglib.remove_single_lineages(tree)
        ancestral = arghmm.emit.parsimony_ancestral_seq(tree, seqs, pos)
        util.print_dict(ancestral)


    def test_thread(self):
        """
        Test thread retrieval
        """
        
        k = 10
        n = 1e4
        rho = 1.5e-8 * 10
        mu = 2.5e-8 * 100
        length = 1000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        for (block, tree), threadi in izip(
            arglib.iter_tree_tracks(arg),
            arghmm.iter_chrom_thread(arg, arg["n9"], by_block=True)):
            print block
            print threadi
            treelib.draw_tree_names(tree.get_tree(), minlen=5, scale=4e-4)


    def test_plot_thread(self):
        """
        Test thread retrieval
        """
        
        k = 60
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(1000e3) / 20
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)

        node = arg.leaves().next()
        x = range(length)
        y = cget(arghmm.iter_chrom_thread(arg, node, by_block=False), 1)
        
        p = plot(x, y, style='lines')

        pause()



    def test_prior(self):
        """
        Calculate state priors
        """

        k = 10
        n = 1e4
        rho = 1.5e-8
        mu = 2.5e-8
        length = 1000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        times = arghmm.get_time_points()
        arghmm.discretize_arg(arg, times)
        new_name = "n%d" % (k-1)
        arg = arghmm.remove_arg_thread(arg, new_name)

        model = arghmm.ArgHmm(arg, seqs, new_name=new_name, times=times)

        prior = [model.prob_prior(0, j)
                 for j in xrange(model.get_num_states(0))]
        print prior
        print sum(map(exp, prior))
        fequal(sum(map(exp, prior)), 1.0, rel=.01)


    def test_trans_single(self):
        """
        Calculate transition probabilities

        Only calculate a single matrix
        """

        k = 4
        n = 1e4
        rho = 1.5e-8
        mu = 2.5e-8
        length = 1000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)
        
        times = arghmm.get_time_points(10)
        arghmm.discretize_arg(arg, times)
        print "recomb", arglib.get_recomb_pos(arg)

        new_name = "n%d" % (k-1)
        arg = arghmm.remove_arg_thread(arg, new_name)
        model = arghmm.ArgHmm(arg, seqs, new_name=new_name, times=times)

        pos = 10
        tree = arg.get_marginal_tree(pos)
        mat = arghmm.calc_transition_probs(
            tree, model.states[pos], model.nlineages,
            model.times, model.time_steps, model.popsizes, rho)
        print model.states[pos]
        pc(mat)

        for row in mat:
            print sum(map(exp, row))

        

    def test_trans_switch_single(self):
        """
        Calculate transitions probabilities for switching between blocks

        Only calculate a single matrix
        """

        k = 5
        n = 1e4
        rho = 1.5e-8 * 100
        mu = 2.5e-8
        length = 1000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        #arglib.write_arg("tmp/a.arg", arg)
        #arg = arglib.read_arg("tmp/a.arg")
        #arg.set_ancestral()

        
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)
        
        times = arghmm.get_time_points(5)
        arghmm.discretize_arg(arg, times)

        new_name = "n%d" % (k-1)
        arg = arghmm.remove_arg_thread(arg, new_name)
        
        model = arghmm.ArgHmm(arg, seqs, new_name=new_name, times=times)

        # get recombs
        recombs = list(x.pos for x in arghmm.iter_visible_recombs(arg))
        print "recomb", recombs

        pos = recombs[0] + 1
        tree = arg.get_marginal_tree(pos-.5)
        last_tree = arg.get_marginal_tree(pos-1-.5)

        print "states1>>", model.states[pos-1]
        print "states2>>", model.states[pos]

        treelib.draw_tree_names(last_tree.get_tree(), minlen=5, maxlen=5)
        treelib.draw_tree_names(tree.get_tree(), minlen=5, maxlen=5)

        print "pos>>", pos
        recomb = [x for x in tree
                  if x.event == "recomb" and x.pos+1 == pos][0]
        mat = arghmm.calc_transition_probs_switch(
            tree, last_tree, recomb.name,
            model.states[pos-1], model.states[pos],
            model.nlineages, model.times,
            model.time_steps, model.popsizes, rho)
        pc(mat)

    '''
    def test_trans_all(self):
        """
        Calculate all transition probabilities
        """

        k = 10
        n = 1e4
        rho = 1.5e-8 * 100
        mu = 2.5e-8
        length = 1000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)
        
        times = arghmm.get_time_points(10)
        arghmm.discretize_arg(arg, times)
        new_name = "n%d" % (k-1)
        arg = arghmm.remove_arg_thread(arg, new_name)
        
        model = arghmm.ArgHmm(arg, seqs, new_name=new_name, times=times)

        # display recombs
        print "recomb", list(x.pos for x in arghmm.iter_visible_recombs(arg))
        print "recomb", model.recomb_pos[1:-1]
        
        
        tree = arg.get_marginal_tree(-.5)
        states = list(arghmm.iter_coal_states(tree, times))
        recomb = arghmm.find_tree_next_recomb(tree, 0)
        
        for pos in xrange(1, length):
            print "site", pos, (recomb.pos if recomb else -1)
            
            if recomb and pos == recomb.pos+1:
                print "switching local trees"
                # start next tree
                last_tree = tree
                tree = arg.get_marginal_tree(pos-.5)

                states = list(arghmm.iter_coal_states(tree, times))
                assert model.states[pos] == states
                
                # we are after local tree change
                mat = arghmm.calc_transition_probs_switch(
                    tree, last_tree, recomb.name,
                    model.states[pos-1], model.states[pos],
                    model.nlineages, model.times,
                    model.time_steps, model.popsizes, model.rho)

                recomb = arghmm.find_tree_next_recomb(arg, pos)

            else:
                assert model.states[pos] == states
                
                mat = arghmm.calc_transition_probs(
                    tree, model.states[pos], model.nlineages,
                    model.times, model.time_steps, model.popsizes, model.rho)
    '''

    def test_trans2(self):
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
        Calculate transition probabilities for k=2

        Only calculate a single matrix
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



    def test_emit_argmax(self):
        """
        Calculate emission probabilities
        """
        
        k = 10
        n = 1e4
        rho = 0.0
        mu = 2.5e-8 * 100
        length = 10000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)
        
        times = arghmm.get_time_points(10)
        arghmm.discretize_arg(arg, times)
        
        new_name = "n%d" % (k-1)
        thread = list(arghmm.iter_chrom_thread(arg, arg[new_name]))
        arg = arghmm.remove_arg_thread(arg, new_name)
        
        model = arghmm.ArgHmm(arg, seqs, new_name=new_name, times=times)
        
        nstates = model.get_num_states(1)
        probs = [0.0 for j in xrange(nstates)]
        for i in xrange(1, length):
            if i % 100 == 0:
                print i
            for j in xrange(nstates):
                probs[j] += model.prob_emission(i, j)
        print

        # is the maximum likelihood emission matching truth
        data = sorted(zip(probs, model.states[0]), reverse=True)
        pc(data[:20])
        
        state = (thread[0][0], times.index(thread[0][1]))
        
        print data[0][1], state
        assert data[0][1] == state



    def test_emit_parsimony(self):
        """
        Calculate emission probabilities with parsimony
        """
        
        k = 10
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(100e3) / 20
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)

        x = []; y = []
        for i in range(20):
            print i
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times)
            seqs = arghmm.make_alignment(arg, muts)

            x.append(arghmm.calc_likelihood(
                arg, seqs, mu=mu, times=times, delete_arg=False))
            y.append(arghmm.calc_likelihood_parsimony(
                arg, seqs, mu=mu, times=times, delete_arg=False))

        p = plot(x, y, xlab="true likelihood", ylab="parsimony likelihood")
        p.plot([min(x), max(x)], [min(x), max(x)], style="lines")
        pause()
            



    def test_trans_c(self):

        k = 10
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 1000
        times = arghmm.get_time_points(20)
        
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        arghmm.discretize_arg(arg, times)
        new_chrom = "n%d" % (k-1)
        arg = arghmm.remove_arg_thread(arg, new_chrom)
        
        print "recomb", arglib.get_recomb_pos(arg)

        model = arghmm.ArgHmm(arg, seqs, new_name="n%d" % (k-1), times=times,
                              mu=mu, rho=rho, popsize=n)

        i = 1
        nstates1 = model.get_num_states(i-1)
        for i in xrange(0, 1):
            print i, nstates1
            nstates2 = model.get_num_states(i)
            states = model.states[i]
            trans = model.transmat
            trans2 = arghmm.calc_transition_probs_c(
                model.local_tree, model.states[i], model.nlineages,
                model.times, model.time_steps, model.popsizes, model.rho)

            print times
            t = model.local_tree.get_tree()
            treelib.remove_single_children(t)
            t.write()

            for a in xrange(nstates1):
                for b in xrange(nstates1):
                    try:
                        print states[a], states[b], trans[a][b], trans2[a][b]
                        fequal(trans[a][b], trans2[a][b])
                    except:
                        print "!", states[a], states[b], trans[a][b], trans2[a][b]
                        #raise
                    
            nstates1 = nstates2


    def test_trans_switch_c(self):

        k = 10
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 5000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        times = arghmm.get_time_points(20)
        arghmm.discretize_arg(arg, times)
        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)
        arg = arglib.smcify_arg(arg)
        print "recomb", arglib.get_recomb_pos(arg)

        model = arghmm.ArgHmm(arg, seqs, new_name="n%d" % (k-1), times=times)

        for pos in model.recomb_pos[:-1][1:50]:
            nstates1 = model.get_num_states(pos)
            print pos, nstates1
            nstates2 = model.get_num_states(pos+1)

            last_tree = model.arg.get_marginal_tree(pos-.5)
            tree = model.arg.get_marginal_tree(pos+1-.5)

            nlineages = arghmm.get_nlineages_recomb_coal(last_tree, model.times)
            recomb = arghmm.find_tree_next_recomb(arg, pos)

            states1 = model.states[pos]
            states2 = model.states[pos+1]

            trans = arghmm.calc_transition_probs_switch(
                tree, last_tree, recomb.name,
                states1, states2,
                nlineages, model.times,
                model.time_steps, model.popsizes, model.rho)

            for a in xrange(nstates1):
                fequal(sum(map(exp, trans[a])), 1.0, rel=.01)

            trans2 = arghmm.calc_transition_probs_switch_c(
                tree, last_tree, recomb.name,
                states1, states2,
                nlineages, model.times,
                model.time_steps, model.popsizes, model.rho, raw=False)

            for a in xrange(nstates1):
                for b in xrange(nstates2):
                    #if trans[a][b] == 0.0:
                    #    print a, states1[a], states2[b], \
                    #          trans[a].index(0.0), trans2[a].index(0.0), \
                    #          trans[a][b] == trans2[a][b]
                    
                    if trans[a][b] in (0.0, -util.INF) or \
                       trans2[a][b] in (0.0, -util.INF):
                        assert trans[a][b] == trans2[a][b], (
                            trans[a][b], trans2[a][b])
                    fequal(trans[a][b], trans2[a][b])

    def test_emit_c(self):
        
        k = 10
        n = 1e4
        rho = 0.0
        mu = 2.5e-8 * 100
        length = 100
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)
        
        times = arghmm.get_time_points(10)
        arghmm.discretize_arg(arg, times)
        
        new_name = "n%d" % (k-1)
        thread = list(arghmm.iter_chrom_thread(arg, arg[new_name]))
        
        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)
        
        model = arghmm.ArgHmm(arg, seqs, new_name=new_name, times=times)
        
        nstates = model.get_num_states(1)
        emit = arghmm.iter_trans_emit_matrices(model, length).next()[4]
        
        for i in xrange(1, length):
            for j in xrange(nstates):
                try:
                    fequal(emit[i][j], model.prob_emission(i, j))
                except:
                    print i, j, emit[i][j], model.prob_emission(i, j)
                    print model.states[i][j]
                    raise


    #========================
    # HMM algorithms

    def test_backward(self):
        """
        Run backward algorithm
        """
        
        k = 3
        n = 1e4
        rho = 1.5e-8 * 100
        mu = 2.5e-8 * 100
        length = 10000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        times = arghmm.get_time_points(ntimes=10)
        arghmm.discretize_arg(arg, times)

        tree = arg.get_marginal_tree(0)
        print tree.root.age
        treelib.draw_tree_names(tree.get_tree(), minlen=5, scale=4e-4)

        # remove chrom
        new_name = "n%d" % (k-1)
        arg = arghmm.remove_arg_thread(arg, new_name)

        model = arghmm.ArgHmm(arg, seqs, new_name=new_name, times=times,
                              rho=rho, mu=mu)
        print "states", len(model.states[0])
        print "recomb", model.recomb_pos
        print "muts", len(muts)

        probs = hmm.backward_algorithm(model, length, verbose=True)

        for pcol in probs:
            p = sum(map(exp, pcol))
            print p, " ".join("%.3f" % f for f in map(exp, pcol))
            


    def test_post(self):

        k = 6
        n = 1e4
        rho = 1.5e-8 * 10
        mu = 2.5e-8 * 10
        length = 10000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)
        print "muts", len(muts)
        print "recombs", len(arglib.get_recomb_pos(arg))
        
        times = arghmm.get_time_points(ntimes=10)
        arghmm.discretize_arg(arg, times)

        tree = arg.get_marginal_tree(0)
        treelib.draw_tree_names(tree.get_tree(), minlen=5, scale=4e-4)

        # remove chrom
        new_name = "n%d" % (k - 1)        
        keep = set(arg.leaf_names()) - set([new_name])
        arglib.subarg_by_leaf_names(arg, keep)
        arg = arglib.smcify_arg(arg)

        model = arghmm.ArgHmm(arg, seqs, new_name=new_name,
                              times=times, rho=rho, mu=mu)
        print "states", len(model.states[0])

        probs = arghmm.get_posterior_probs(model, length, verbose=True)
        
        for pcol in probs:
            p = sum(map(exp, pcol))
            print p, " ".join("%.3f" % f for f in map(exp, pcol))
            fequal(p, 1.0, rel=1e-2)
                


    def test_post2(self):

        k = 2
        n = 1e4
        rho = 1.5e-8 * 10
        mu = 2.5e-8 * 10
        length = 10000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        print "muts", len(muts)

        times = arghmm.get_time_points()
        arghmm.discretize_arg(arg, times)

        thread = list(arghmm.iter_chrom_thread(arg, arg["n1"], by_block=False))
        tree = arg.get_marginal_tree(0)
        print tree.root.age
        treelib.draw_tree_names(tree.get_tree(), minlen=5, scale=4e-4)
        p = plot(cget(thread, 1), style="lines", ymin=0)

        #alignlib.print_align(seqs)

        # remove chrom
        keep = ["n0"]
        arglib.subarg_by_leaf_names(arg, keep)
        arg = arglib.smcify_arg(arg)

        model = arghmm.ArgHmm(arg, seqs, new_name="n1", times=times,
                              rho=rho, mu=mu)
        print "states", len(model.states[0])

        probs = arghmm.get_posterior_probs(model, length, verbose=True)
        
        high = list(arghmm.iter_posterior_times(model, probs, .95))
        low = list(arghmm.iter_posterior_times(model, probs, .05))
        p.plot(high, style="lines")
        p.plot(low, style="lines")

        pause()


    def test_post3(self):

        k = 3
        n = 1e4
        rho = 1.5e-8 * 3
        mu = 2.5e-8 * 100
        length = 10000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        arg.prune()
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        times = arghmm.get_time_points(ntimes=10)
        arghmm.discretize_arg(arg, times)


        tree = arg.get_marginal_tree(0)
        treelib.draw_tree_names(tree.get_tree(), minlen=5, scale=4e-4)

        thread = list(arghmm.iter_chrom_thread(arg, arg["n2"], by_block=False))
        p = plot(cget(thread, 1), style="lines", ymin=0)

        # remove chrom
        keep = ["n0", "n1"]
        arglib.subarg_by_leaf_names(arg, keep)
        arg.set_ancestral()
        arg.prune()

        
        model = arghmm.ArgHmm(arg, seqs, new_name="n2", times=times,
                              rho=rho, mu=mu)
        print "states", len(model.states[0])        
        print "muts", len(muts)
        print "recomb", len(model.recomb_pos) - 2, model.recomb_pos[1:-1]


        p.plot(model.recomb_pos, [1000] * len(model.recomb_pos),
               style="points")

        probs = arghmm.get_posterior_probs(model, length, verbose=True)
        
        high = list(arghmm.iter_posterior_times(model, probs, .95))
        low = list(arghmm.iter_posterior_times(model, probs, .05))
        p.plot(high, style="lines")
        p.plot(low, style="lines")

        pause()



    def test_post_plot(self):

        k = 6
        n = 1e4
        rho = 1.5e-8 * 50
        mu = 2.5e-8 * 50
        length = 10000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        times = arghmm.get_time_points(ntimes=30)
        arghmm.discretize_arg(arg, times)

        # save
        #arglib.write_arg("test/data/k4.arg", arg)
        #fasta.write_fasta("test/data/k4.fa", seqs)

        new_name = "n%d" % (k-1)
        thread = list(arghmm.iter_chrom_thread(arg, arg[new_name],
                                               by_block=False))    
        p = plot(cget(thread, 1), style="lines", ymin=times[1],
                 ylog=10)

        # remove chrom
        new_name = "n%d" % (k-1)
        arg = arghmm.remove_arg_thread(arg, new_name)

        model = arghmm.ArgHmm(arg, seqs, new_name=new_name, times=times,
                              rho=rho, mu=mu)
        print "states", len(model.states[0])
        print "muts", len(muts)
        print "recomb", len(model.recomb_pos) - 2, model.recomb_pos[1:-1]

        p.plot(model.recomb_pos, [10000] * len(model.recomb_pos),
               style="points")

        probs = arghmm.get_posterior_probs(model, length, verbose=True)
        print "done"
        
        high = list(arghmm.iter_posterior_times(model, probs, .95))
        low = list(arghmm.iter_posterior_times(model, probs, .05))
        p.gnuplot("set linestyle 2")
        p.plot(high, style="lines")
        p.gnuplot("set linestyle 2")
        p.plot(low, style="lines")


        write_list("test/data/post_real.txt", cget(thread, 1))
        write_list("test/data/post_high.txt", high)
        write_list("test/data/post_low.txt", low)

        pause()


    def test_norecomb_plot(self):

        k = 50
        n = 1e4
        rho = 1.5e-8 * .0001
        rho2 = 1.5e-8 * 10
        mu = 2.5e-8 * 100
        length = 1000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)


        times = arghmm.get_time_points(ntimes=20)
        arghmm.discretize_arg(arg, times)

        tree = arg.get_marginal_tree(0)
        treelib.draw_tree_names(tree.get_tree(), minlen=5, scale=4e-4)
        
        # get thread
        new_name = "n%d" % (k-1)
        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)
        arg.set_ancestral()
        arg.prune()

        model = arghmm.ArgHmm(arg, seqs, new_name=new_name, times=times,
                              rho=rho2, mu=mu)
        print "states", len(model.states[0])
        print "muts", len(muts)

        # simulate a new thread
        states = list(islice(hmm.sample_hmm_states(model), 0, arg.end))
        data = list(hmm.sample_hmm_data(model, states))
        
        seqs[new_name] = "".join(data)
        #alignlib.print_align(seqs)

        thread = [model.times[model.states[i][s][1]]
                  for i, s in enumerate(states)]
        p = plot(thread, style="lines")
        

        probs = arghmm.get_posterior_probs(model, length, verbose=True)
        print "done"
        
        high = list(arghmm.iter_posterior_times(model, probs, .75))
        low = list(arghmm.iter_posterior_times(model, probs, .25))
        p.plot(high, style="lines")
        p.plot(low, style="lines")

        pause()
        


    def test_post_real(self):

        k = 3
        n = 1e4
        rho = 1.5e-8
        mu = 2.5e-8
        length = 100000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)
        
        #arg = arglib.read_arg("test/data/real.arg")
        #seqs = fasta.read_fasta("test/data/real.fa")

        #arglib.write_arg("test/data/real.arg", arg)
        #fasta.write_fasta("test/data/real.fa", seqs)

        times = arghmm.get_time_points(maxtime=50000, ntimes=20)
        arghmm.discretize_arg(arg, times)

        new_name = "n%d" % (k - 1)
        thread = list(arghmm.iter_chrom_thread(arg, arg[new_name],
                                               by_block=False))
        tree = arg.get_marginal_tree(0)
        print tree.root.age
        treelib.draw_tree_names(tree.get_tree(), minlen=5, scale=4e-4)
        p = plot(cget(thread, 1), style="lines", ymin=10, ylog=10)

        #alignlib.print_align(seqs)

        # remove chrom
        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)
        arg = arglib.smcify_arg(arg)

        model = arghmm.ArgHmm(arg, seqs, new_name=new_name, times=times,
                              rho=rho, mu=mu)
        
        print "states", len(model.states[0])
        #print "muts", len(muts)
        print "recomb", len(model.recomb_pos) - 2, model.recomb_pos[1:-1]        

        probs = arghmm.get_posterior_probs(model, length, verbose=True)
        
        high = list(arghmm.iter_posterior_times(model, probs, .95))
        low = list(arghmm.iter_posterior_times(model, probs, .05))
        p.plot(high, style="lines")
        p.plot(low, style="lines")

        pause()



    def test_determ(self):

        k = 8
        n = 1e4
        rho = 1.5e-8
        mu = 2.5e-8
        length = 100000
        
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)
        times = arghmm.get_time_points(maxtime=50000, ntimes=20)
        arghmm.discretize_arg(arg, times)

        new_name = "n%d" % (k - 1)
        thread = list(arghmm.iter_chrom_thread(arg, arg[new_name],
                                               by_block=False))
        thread_clades = list(arghmm.iter_chrom_thread(
            arg, arg[new_name], by_block=True, use_clades=True))

        # remove chrom
        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)
        arg = arglib.smcify_arg(arg)

        model = arghmm.ArgHmm(arg, seqs, new_name=new_name, times=times,
                              rho=rho, mu=mu)


        for i, rpos in enumerate(model.recomb_pos[1:-1]):
            pos = rpos + 1

            model.check_local_tree(pos, force=True)

            #recomb = arghmm.find_tree_next_recomb(arg, pos - 1)
            tree = arg.get_marginal_tree(pos-.5)
            last_tree = arg.get_marginal_tree(pos-1-.5)
            states1 = model.states[pos-1]
            states2 = model.states[pos]
            (recomb_branch, recomb_time), (coal_branch, coal_time) = \
                arghmm.find_recomb_coal(tree, last_tree, pos=rpos)
            recomb_time = times.index(recomb_time)
            coal_time = times.index(coal_time)

            determ = arghmm.get_deterministic_transitions(
                states1, states2, times,
                tree, last_tree,
                recomb_branch, recomb_time,
                coal_branch, coal_time)


            leaves1, time1, block1 = thread_clades[i]
            leaves2, time2, block2 = thread_clades[i+1]
            if new_name in leaves1:
                leaves1.remove(new_name)
            if new_name in leaves2:
                leaves2.remove(new_name)
            node1 = arghmm.arg_lca(arg, leaves1, None, pos-1).name
            node2 = arghmm.arg_lca(arg, leaves2, None, pos).name

            state1 = (node1, times.index(time1))
            state2 = (node2, times.index(time2))
            print pos, state1, state2
            try:
                statei1 = states1.index(state1)
                statei2 = states2.index(state2)
            except:
                print "states1", states1
                print "states2", states2
                raise

            statei3 = determ[statei1]
            print "  ", statei1, statei2, statei3, states2[statei3]
            if statei2 != statei3 and statei3 != -1:
                tree = tree.get_tree()
                treelib.remove_single_children(tree)
                last_tree = last_tree.get_tree()
                treelib.remove_single_children(last_tree)

                print "block1", block1
                print "block2", block2
                print "r=", (recomb_branch, recomb_time)
                print "c=", (coal_branch, coal_time)

                print "tree"
                treelib.draw_tree_names(tree, minlen=8, maxlen=8)
            
                print "last_tree"
                treelib.draw_tree_names(last_tree, minlen=8, maxlen=8)
                assert False

    #=========================================================================


    def test_forward_c(self):

        k = 20
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(200e3 / 20)
        times = arghmm.get_time_points(ntimes=20)

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
            #sum2 = stats.logsum(col2)
            #col2 = [exp(x-sum2) for x in col2]
            
            for a, b in izip(col1, col2):
                try:
                    #print a, b
                    fequal(a, b, rel=.0001)
                except:
                    print model.states[i]
                    print i, col1
                    print i, col2
                    raise


    '''
    def test_forward_c(self):

        k = 10
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 5000
        arg = arglib.sample_arg(k, 2*n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)
        
        print "muts", len(muts)
        print "recomb", len(arglib.get_recomb_pos(arg))

        times = arghmm.get_time_points(ntimes=20)
        arghmm.discretize_arg(arg, times)

        # remove chrom
        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)
        arg = arglib.smcify_arg(arg)

        model = arghmm.ArgHmm(arg, seqs, new_name="n%d" % (k-1), times=times,
                              rho=rho, mu=mu)
        print "states", len(model.states[0])
        
        util.tic("C")
        probs1 = list(arghmm.forward_algorithm(model, length, verbose=True))
        util.toc()
        
        util.tic("python")
        probs2 = list(arghmm.py_forward_algorithm2(model, length, verbose=True))
        util.toc()

        for i, (col1, col2) in enumerate(izip(probs1, probs2)):
            sum2 = stats.logsum(col2)
            col2 = [exp(x-sum2) for x in col2]
            
            for a, b in izip(col1, col2):
                try:
                    fequal(a, b, rel=.01)
                except:
                    print model.states[i]
                    print i, col1
                    print i, col2
                    raise
                '''

    def test_backward_c(self):

        k = 3
        n = 1e4
        rho = 1.5e-8 * 30
        mu = 2.5e-8 * 100
        length = 100
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        arg.prune()
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        print arglib.get_recomb_pos(arg)
        print "muts", len(muts)
        print "recomb", len(arglib.get_recomb_pos(arg))

        times = arghmm.get_time_points(ntimes=10)
        arghmm.discretize_arg(arg, times)

        tree = arg.get_marginal_tree(0)
        print tree.root.age
        treelib.draw_tree_names(tree.get_tree(), minlen=5, scale=4e-4)

        # remove chrom
        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)

        model = arghmm.ArgHmm(arg, seqs, new_name="n%d" % (k-1), times=times,
                              rho=rho, mu=mu)
        print "states", len(model.states[0])

        util.tic("C")
        probs1 = list(arghmm.backward_algorithm(model, length, verbose=True))
        util.toc()     

        util.tic("python")
        probs2 = list(hmm.backward_algorithm(model, length, verbose=True))
        util.toc()

        print "probs1"
        pc(probs1)

        print "probs2"
        pc(probs2)
        

        for col1, col2 in izip(probs1, probs2):
            for a, b in izip(col1, col2):
                fequal(a, b)
        


    def test_post_c(self):

        k = 3
        n = 1e4
        rho = 1.5e-8 * 30
        mu = 2.5e-8 * 100
        length = 100
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        arg.prune()
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        print arglib.get_recomb_pos(arg)
        print "muts", len(muts)
        print "recomb", len(arglib.get_recomb_pos(arg))

        times = arghmm.get_time_points(ntimes=10)
        arghmm.discretize_arg(arg, times)

        tree = arg.get_marginal_tree(0)
        print tree.root.age
        treelib.draw_tree_names(tree.get_tree(), minlen=5, scale=4e-4)

        # remove chrom
        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)

        model = arghmm.ArgHmm(arg, seqs, new_name="n%d" % (k-1), times=times,
                              rho=rho, mu=mu)
        print "states", len(model.states[0])

        util.tic("C")
        probs1 = list(arghmm.get_posterior_probs(model, length, verbose=True))
        util.toc()

        util.tic("python")
        probs2 = list(hmm.get_posterior_probs(model, length, verbose=True))
        util.toc()

        print "probs1"
        pc(probs1)

        print "probs2"
        pc(probs2)
        

        for col1, col2 in izip(probs1, probs2):
            for a, b in izip(col1, col2):
                fequal(a, b)


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



    def test_arg_convert(self):
        """
        Test conversion for python to C args
        """

        k = 10
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)

        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)

        # convert to C++ and back
        trees, names = arghmm.arg2ctrees(arg, times)
        arg2 = arghmm.ctrees2arg(trees, names, times)

        test_arg_equal(arg, arg2)


    #======================================================================
    # misc

    def test_state_corr(self):

        k = 12
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(1e3) / 20
        times = arghmm.get_time_points(ntimes=20, maxtime=200e3)
        
        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg, mu, times)
        seqs = arglib.make_alignment(arg, muts)
                
        # remove chrom
        new_name = "n%d" % (k-1)
        arg = arghmm.remove_arg_thread(arg, new_name)
        
        model = arghmm.ArgHmm(arg, seqs, new_name=new_name, times=times,
                              rho=rho, mu=mu)
        print "states", len(model.states[0])


        nstates = len(model.states[0])
        prior = [-util.INF] * nstates
        prior[random.randint(0, nstates)] = 0.0

        probs1 = list(arghmm.forward_algorithm(model, length, verbose=True))
        probs2 = list(arghmm.forward_algorithm(model, length, prior=prior,
                                               verbose=True))

        model.rho *= 1e-9
        probs3 = list(arghmm.forward_algorithm(model, length, prior=prior,
                                               verbose=True))

        p = plot(vsubs(probs1[length-1], mean(probs1[length-1])))
        p.plot(vsubs(probs2[length-1], mean(probs2[length-1])))
        p.plot(vsubs(probs3[length-1], mean(probs3[length-1])))

        pause()


    def test_compress_align(self):
        """Test the compression of sequence alignments"""

        k = 12
        n = 1e4
        rho = 1.5e-8
        mu = 2.5e-8
        length = 200e3
        times = arghmm.get_time_points(ntimes=20, maxtime=200e3)
        compress = 20
        
        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg, mu, times)
        seqs = arglib.make_alignment(arg, muts)

        seqs2, cols = arghmm.compress_align(seqs, compress)
        print seqs2.alignlen(), length / compress
        delta = [cols[i] - cols[i-1] for i in range(1, len(cols))]
        
        plot(cols)
        plothist(delta, width=1)
        
        variant = [arghmm.is_variant(seqs, i) for i in range(seqs.alignlen())]
        print histtab(variant)
        print histtab(mget(variant, cols))
        
        pause()


    #------------------------------------
    # LD

    def test_ld(self):

        col1 = "AAAAAACCCCCCCCC"
        col2 = "AAAACCCCCCCCCCC"

        self.assertEqual(arghmm.find_high_freq_allele(col1), "C")
        self.assertEqual(arghmm.find_high_freq_allele(col2), "C")

        print_dict(arghmm.find_pair_allele_freqs(col1, col2))

        print arghmm.calc_ld_D(col1, col2)
        print arghmm.calc_ld_Dp(col1, col2)
        print arghmm.calc_ld_r2(col1, col2)


        # case 2
        col1 = "AAAAAACCCCCCCCC"
        col2 = "AAAACCACCCCCCCC"

        self.assertEqual(arghmm.find_high_freq_allele(col1), "C")
        self.assertEqual(arghmm.find_high_freq_allele(col2), "C")

        print_dict(arghmm.find_pair_allele_freqs(col1, col2))

        print arghmm.calc_ld_D(col1, col2)
        print arghmm.calc_ld_Dp(col1, col2)
        print arghmm.calc_ld_r2(col1, col2)



    def test_ld_block(self):

        k = 30
        n = 1e4
        rho = 1.5e-8
        mu = 2.5e-8
        length = 200e3
        times = arghmm.get_time_points(ntimes=20, maxtime=200e3)
        compress = 20
        
        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg, mu, times)
        seqs = arghmm.make_alignment(arg, muts)
        sites = arghmm.seqs2sites(seqs)

        #cols = transpose(seqs.values())[::10000]
        cols = mget(sites, sites.positions)
        cols = cols[:1000]

        ld = arghmm.calc_ld_matrix(cols, arghmm.calc_ld_Dp)
        
        heatmap(ld, width=2, height=2)
        pause()

        
        

#=============================================================================
if __name__ == "__main__":

    test_main()


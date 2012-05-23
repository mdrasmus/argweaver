
import unittest, random

import arghmm

from rasmus.common import *
from rasmus import stats, hmm
from rasmus.testing import *

from compbio import coal, arglib, fasta



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

        k = 10
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
        
        times = arghmm.get_time_points(10)
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
        
        
        for row in mat:
            fequal(sum(map(exp, row)), 1.0)


    def test_transprobs_all(self):
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


            for row in mat:
                fequal(sum(map(exp, row)), 1.0)

        
    def test_trans(self):
        """
        Calculate transition probabilities
        """

        k = 10
        n = 1e4
        rho = 1.5e-8
        mu = 2.5e-8
        length = 1000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        times = arghmm.get_time_points(20)
        arghmm.discretize_arg(arg, times)
        new_name = "n%d" % (k - 1)
        arg = arghmm.remove_arg_thread(arg, new_name)

        model = arghmm.ArgHmm(arg, seqs, new_name=new_name, times=times)
        print "recomb", model.recomb_pos[1:-1]

        i = 1
        nstates1 = model.get_num_states(i-1)
        nstates2 = model.get_num_states(i)
        treelen = sum(x.get_dist() for x in model.local_tree)
        for a in xrange(nstates1):
            trans = [model.prob_transition(i-1, a, i, b)
                     for b in xrange(nstates2)]
            #print trans
            print a, sum(map(exp, trans))
            print map(exp, trans)
            fequal(sum(map(exp, trans)), 1.0, rel=.01)

            # is non-transition greater than no-recomb prob
            node, timei = model.states[i][a]
            blen = model.times[timei]
            treelen2 = treelen + blen
            if node == model.local_tree.root.name:
                treelen2 += blen - model.local_tree.root.age
            norecomb = -model.rho * (treelen2 - treelen)
            
            print trans[a], norecomb, norecomb <= trans[a]
            

    def test_trans_switch(self):
        """
        Calculate transitions probabilities for switching between blocks
        """

        k = 10
        n = 1e4
        rho = 1.5e-8
        mu = 2.5e-8
        length = 1000
        arg = arglib.sample_arg(k, n, rho*100, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        times = arghmm.get_time_points(20)
        arghmm.discretize_arg(arg, times)
        new_name = "n%d" % (k - 1)
        arg = arghmm.remove_arg_thread(arg, new_name)

        model = arghmm.ArgHmm(arg, seqs, new_name=new_name, times=times)
        print "recomb", model.recomb_pos[1:-1]

        i = 1
        nstates1 = model.get_num_states(i-1)
        for i in xrange(1, length):
            print i, nstates1
            nstates2 = model.get_num_states(i)
            for a in xrange(nstates1):
                trans = [model.prob_transition(i-1, a, i, b)
                         for b in xrange(nstates2)]
                fequal(sum(map(exp, trans)), 1.0, rel=.01)
            nstates1 = nstates2


    def test_emit(self):
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



    def test_trans_c(self):

        k = 10
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 1000
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

        i = 1
        nstates1 = model.get_num_states(i-1)
        for i in xrange(0, 1):
            print i, nstates1
            nstates2 = model.get_num_states(i)
            for a in xrange(nstates1):
                trans = [model.prob_transition(i-1, a, i, b)
                         for b in xrange(nstates2)]
                fequal(sum(map(exp, trans)), 1.0, rel=.01)
            trans = model.transmat
            trans2 = arghmm.calc_transition_probs_c(
                model.local_tree, model.states[i], model.nlineages,
                model.times, model.time_steps, model.popsizes, model.rho)

            for a in xrange(nstates1):
                for b in xrange(nstates1):
                    fequal(trans[a][b], trans2[a][b])
                    
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


    def test_viterbi(self):
        """
        Run viterbi algorithm on HMM
        """

        k = 10
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)
        print "muts", len(muts)
        print "recombs", len(arglib.get_recomb_pos(arg))

        times = arghmm.get_time_points(ntimes=20)
        arghmm.discretize_arg(arg, times)

        new_name = "n%d" % (k-1)
        thread = list(arghmm.iter_chrom_thread(arg, arg[new_name], by_block=False))
        tree = arg.get_marginal_tree(0)
        print tree.root.age
        treelib.draw_tree_names(tree.get_tree(), minlen=5, scale=4e-4)
        p = plot(cget(thread, 1), style="lines", ymin=0)
        

        # remove chrom
        arg = arghmm.remove_arg_thread(arg, new_name)
        
        model = arghmm.ArgHmm(arg, seqs, new_name=new_name, times=times,
                              rho=rho, mu=mu)
        print "states", len(model.states[0])

        path = hmm.viterbi(model, length, verbose=True)
        thread2 = [model.states[i][j] for i, j in enumerate(path)]
        thread2 = [(node, model.times[t]) for node, t in thread2]
        p.plot(cget(thread2, 1), style="lines")

        pause()


    '''
    def test_viterbi_node(self):


        k = 10
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 1000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        times = arghmm.get_time_points(ntimes=10)
        arghmm.discretize_arg(arg, times)

        tree = arg.get_marginal_tree(0)
        treelib.draw_tree_names(tree.get_tree(), minlen=5, scale=4e-4)

        # TODO: need to recode correct thread for comparison

        # plot true thread chrom
        new_name = "n%d" % (k-1)
        thread = list((node, times.index(age))
                      for (node, age) in arghmm.iter_chrom_thread(
            arg, arg[new_name], by_block=False))

        # remove chrom
        arg = arghmm.remove_arg_thread(arg, new_name)
        
        model = arghmm.ArgHmm(arg, seqs, new_name=new_name, times=times,
                              rho=rho, mu=mu)
        print "states", len(model.states[0])
        print "muts", len(muts)
        print "recombs", len(model.recomb_pos) - 2

        path = hmm.viterbi(model, length, verbose=True)

        node2id = dict((x.name, i) for i, x in enumerate(arg))

        # plot inferred thread
        thread2 = [model.states[i][j] for i, j in enumerate(path)]

        nodes1 = [node2id[p[0]] for p in thread]
        nodes2 = [node2id[p[0]] for p in thread2]

        p = plot(nodes1, style="lines")
        p.plot(nodes2, style="lines")

        for i, (a, b) in enumerate(izip(thread, thread2)):
            if i % 100 == 0:
                print a, b

        pause()
    '''


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

        k = 7
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

        tree = arg.get_marginal_tree(0)
        treelib.draw_tree_names(tree.get_tree(), minlen=5, scale=4e-4)


        new_name = "n%d" % (k-1)
        thread = list(arghmm.iter_chrom_thread(arg, arg[new_name],
                                               by_block=False))    
        p = plot(cget(thread, 1), style="lines", ymin=30,
                 ylog=10)

        # remove chrom
        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)
        arg = arglib.smcify_arg(arg)

        model = arghmm.ArgHmm(arg, seqs, new_name="n%d" % (k-1), times=times,
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

        k = 10
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        print arglib.get_recomb_pos(arg)
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
            for a, b in izip(col1, col2):
                try:
                    fequal(a, b, rel=.01)
                except:
                    print model.states[i]
                    print i, col1
                    print i, col2
                    raise
        

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


    def test_forward2(self):

        k = 5
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)        
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)
        

        print arglib.get_recomb_pos(arg)
        print "muts", len(muts)
        print "recomb", len(arglib.get_recomb_pos(arg))

        times = arghmm.get_time_points(ntimes=30)
        arghmm.discretize_arg(arg, times)

        # remove chrom
        new_name = "n%d" % (k-1)
        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)
        arg = arglib.smcify_arg(arg)

        model = arghmm.ArgHmm(arg, seqs, new_name=new_name, times=times,
                              rho=rho, mu=mu)
        print "states", len(model.states[0])

        
        arghmm.forward_algorithm2(model, length, verbose=True)

        




#=============================================================================
if __name__ == "__main__":

    test_main()


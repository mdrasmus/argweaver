# test dlcoal.duploss

import unittest, random

import arghmm

from rasmus.common import *
from rasmus import stats, hmm
from rasmus.testing import *

from compbio import coal, arglib, fasta


class Basic (unittest.TestCase):

    def test_prune(self):

        k = 3
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8
        length = 1000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        #arglib.write_arg("tmp/a.arg", arg)
        #arg = arglib.read_arg("tmp/a.arg")
        
        times = arghmm.get_time_points(10)
        arghmm.discretize_arg(arg, times)

        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)
        arg.set_ancestral()
        arg.prune()
        
        recomb = arglib.get_recomb_pos(arg)
        print "recomb", recomb

        pos = recomb[0]
        tree = arg.get_marginal_tree(pos+.5)

        draw_tree_names(tree.get_tree(), minlen=5, maxlen=5)
        
        for row in sorted((x.pos, x.name, x.parents, x.data["ancestral"])
                          for x in arg if x.event == "recomb"):
            print row
        
        print sorted([(x.pos, x.event) for x in tree])
        nodes = [x for x in tree
                 if x.event == "recomb" and x.pos == pos]
        print nodes
        assert len(nodes) == 1


    def test_recomb(self):

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
            recomb = arghmm.find_tree_next_recomb(tree, i+1)
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
        k = 10
        n = 1e4
        rho = 1.5e-8 * 100
        length = 1000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)

        times = arghmm.get_time_points()
        arghmm.discretize_arg(arg, times)
        tree = arg.get_marginal_tree(0)
        nlineages = arghmm.get_nlineages(tree, times)

        treelib.draw_tree_names(tree.get_tree(), scale=4e-4)

        print list(arghmm.iter_coal_states(tree, times))
        print nlineages
        #self.assert_(nlineages == sorted(nlineages, reverse=True))


    def test_states(self):
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
            assert (tree.root.name, 10) in states
            print states


    def test_pars_seq(self):

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

        ancestral = arghmm.parsimony_ancestral_seq(tree, seqs, pos)
        util.print_dict(ancestral)


    def test_thread(self):
        
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

        


    def test_add_thread(self):
        
        k = 10
        n = 1e4
        rho = 1.5e-8 * 10
        mu = 2.5e-8 * 100
        length = 100
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)

        thread = arghmm.iter_chrom_thread(arg, arg["n9"], by_block=False)



    def test_prior(self):

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

        model = arghmm.ArgHmm(arg, seqs, new_name="n%d" % k)

        prior = [model.prob_prior(0, j)
                 for j in xrange(model.get_num_states(0))]
        print prior
        print sum(map(exp, prior))
        fequal(sum(map(exp, prior)), 1.0, rel=.01)


    def test_transprobs(self):

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

        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)        
        model = arghmm.ArgHmm(arg, seqs, new_name="n%d" % (k-1), times=times)

        pos = 10
        tree = arg.get_marginal_tree(pos)
        mat = arghmm.calc_transition_probs(
            tree, model.states[pos], model.nlineages,
            model.times, model.time_steps, model.popsizes, rho)
        print model.states[pos]
        pc(mat)

        for row in mat:
            print sum(map(exp, row))



    def test_recomb_coal(self):

        k = 10
        n = 1e4
        rho = 1.5e-8 * 100
        mu = 2.5e-8
        length = 1000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        arglib.write_arg("tmp/a.arg", arg)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)
        
        #arg = arglib.read_arg("tmp/a.arg")
        #arg.set_ancestral()
        #find_recomb_coal(tree, last_tree, recomb_name=None, pos=None)
        
        times = arghmm.get_time_points(10)
        arghmm.discretize_arg(arg, times)

        # get recombs
        recombs = list(x.pos for x in arghmm.iter_visible_recombs(arg))
        print "recomb", recombs

        pos = recombs[0] + 1
        tree = arg.get_marginal_tree(pos-.5)
        last_tree = arg.get_marginal_tree(pos-1-.5)

        treelib.draw_tree_names(last_tree.get_tree(), minlen=5, maxlen=5)
        treelib.draw_tree_names(tree.get_tree(), minlen=5, maxlen=5)

        model = arghmm.ArgHmm(arg, seqs, new_name="n%d" % (k-1), times=times)



        

    def test_transprobs_switch(self):

        k = 10
        n = 1e4
        rho = 1.5e-8 * 100
        mu = 2.5e-8
        length = 1000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        arglib.write_arg("tmp/a.arg", arg)
        #arg = arglib.read_arg("tmp/a.arg")
        #arg.set_ancestral()

        
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)
        
        times = arghmm.get_time_points(10)
        arghmm.discretize_arg(arg, times)

        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)
        arg.set_ancestral()
        arg.prune()
        
        model = arghmm.ArgHmm(arg, seqs, new_name="n%d" % (k-1), times=times)

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

        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)
        arg.set_ancestral()
        arg.prune()
        
        model = arghmm.ArgHmm(arg, seqs, new_name="n%d" % (k-1), times=times)

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
        print "recomb", arglib.get_recomb_pos(arg)

        model = arghmm.ArgHmm(arg, seqs, new_name="n%d" % k, times=times)

        i = 1
        nstates1 = model.get_num_states(i-1)
        nstates2 = model.get_num_states(i)
        for a in xrange(nstates1):
            trans = [model.prob_transition(i-1, a, i, b)
                     for b in xrange(nstates2)]
            print a
            print trans
            print sum(map(exp, trans))
            fequal(sum(map(exp, trans)), 1.0, rel=.01)
            

    def test_trans2(self):

        # test with recombinations (default_state code)

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
        print "recomb", arglib.get_recomb_pos(arg)

        model = arghmm.ArgHmm(arg, seqs, new_name="n%d" % k, times=times)

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
        
        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)
        
        model = arghmm.ArgHmm(arg, seqs, new_name=new_name, times=times)
        
        nstates = model.get_num_states(1)
        probs = [0.0 for j in xrange(nstates)]
        for i in xrange(1, length):
            for j in xrange(nstates):
                probs[j] += model.prob_emission(i, j)
        
        data = sorted(zip(probs, model.states[0]), reverse=True)
        pc(data[:20])
        
        state = (thread[0][0], times.index(thread[0][1]))
        
        print data[0][1], state
        assert data[0][1] == state



    def test_viterbi(self):

        k = 10
        n = 1e4
        rho = 1.5e-8 
        mu = 2.5e-8 * 10
        length = 100000
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
        keep = set(arg.leaf_names()) - set([new_name])
        arglib.subarg_by_leaf_names(arg, keep)

        model = arghmm.ArgHmm(arg, seqs, new_name=new_name, times=times,
                              rho=rho, mu=mu)
        print "states", len(model.states[0])

        path = hmm.viterbi(model, length, verbose=True)
        thread2 = [model.states[i][j] for i, j in enumerate(path)]
        thread2 = [(node, model.times[t]) for node, t in thread2]
        p.plot(cget(thread2, 1), style="lines")

        pause()



    def test_viterbi2(self):

        k = 2
        n = 1e4
        rho = 1.5e-8
        mu = 2.5e-8 * 10
        length = 10000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)
        print "muts", len(muts)
        print "recomb", len(arglib.get_recomb_pos(arg))

        times = arghmm.get_time_points()
        arghmm.discretize_arg(arg, times)

        thread = list(arghmm.iter_chrom_thread(arg, arg["n1"], by_block=False))
        tree = arg.get_marginal_tree(0)
        print tree.root.age
        treelib.draw_tree_names(tree.get_tree(), minlen=5, scale=4e-4)
        p = plot(cget(thread, 1), style="lines", ymin=0)


        # remove chrom
        keep = ["n0"]
        arglib.subarg_by_leaf_names(arg, keep)

        model = arghmm.ArgHmm(arg, seqs, new_name="n1", times=times,
                              rho=rho, mu=mu)
        print "states", len(model.states[0])

        path = hmm.viterbi(model, length, verbose=True)
        thread2 = [model.states[i][j] for i, j in enumerate(path)]
        thread2 = [(node, model.times[t]) for node, t in thread2]
        
        p.plot(cget(thread2, 1), style="lines")

        pause()



    def test_viterbi3(self):

        k = 3
        n = 1e4
        rho = 1.5e-8 * 3
        mu = 2.5e-8 * 100
        length = 10000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        print "muts", len(muts)

        times = arghmm.get_time_points(10)
        arghmm.discretize_arg(arg, times)

        thread = list(arghmm.iter_chrom_thread(arg, arg["n2"], by_block=False))
        tree = arg.get_marginal_tree(0)
        print tree.root.age
        treelib.draw_tree_names(tree.get_tree(), minlen=5, scale=4e-4)
        p = plot(cget(thread, 1), style="lines", ymin=0)

        #alignlib.print_align(seqs)

        # remove chrom
        keep = ["n0", "n1"]
        arglib.subarg_by_leaf_names(arg, keep)

        model = arghmm.ArgHmm(arg, seqs, new_name="n2", times=times,
                              rho=rho, mu=mu)
        
        print model.states[0]

        path = hmm.viterbi(model, length, verbose=True)
        thread2 = [model.states[i][j] for i, j in enumerate(path)]
        thread2 = [(node, model.times[t]) for node, t in thread2]
        
        p.plot(cget(thread2, 1), style="lines")

        pause()


    def test_viterbi_node(self):

        k = 20
        n = 1e4
        rho = 1.5e-8 
        mu = 2.5e-8
        length = 20000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        times = arghmm.get_time_points(ntimes=10)
        arghmm.discretize_arg(arg, times)

        tree = arg.get_marginal_tree(0)
        treelib.draw_tree_names(tree.get_tree(), minlen=5, scale=4e-4)

        # plot true thread chrom
        new_name = "n%d" % (k-1)
        thread = list((node, times.index(age))
                      for (node, age) in arghmm.iter_chrom_thread(
            arg, arg[new_name], by_block=False))

        # remove chrom
        keep = set(arg.leaf_names()) - set([new_name])
        arglib.subarg_by_leaf_names(arg, keep)
        arg.set_ancestral()
        arg.prune()

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



    def test_backward(self):

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
        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)

        arg.set_ancestral()
        arg.prune()

        model = arghmm.ArgHmm(arg, seqs, new_name="n%d" % (k-1), times=times,
                              rho=rho, mu=mu)
        print "states", len(model.states[0])
        print "recomb", model.recomb_pos
        print "muts", len(muts)

        probs = hmm.backward_algorithm(model, length, verbose=True)

        for pcol in probs:
            p = sum(map(exp, pcol))
            print p, " ".join("%.3f" % f for f in map(exp, pcol))
            


    def test_post(self):

        k = 3
        n = 1e4
        rho = 1.5e-8 * 3
        mu = 2.5e-8 * 100
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
        arg.set_ancestral()
        arg.prune()

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
        rho = 1.5e-8 * 3
        mu = 2.5e-8 * 100
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

        k = 5
        n = 1e4
        rho = 1.5e-8 * 10
        mu = 2.5e-8 * 100
        length = 10000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        times = arghmm.get_time_points(ntimes=20)
        arghmm.discretize_arg(arg, times)

        # save
        arglib.write_arg("test/data/k4.arg", arg)
        fasta.write_fasta("test/data/k4.fa", seqs)

        tree = arg.get_marginal_tree(0)
        treelib.draw_tree_names(tree.get_tree(), minlen=5, scale=4e-4)


        new_name = "n%d" % (k-1)
        thread = list(arghmm.iter_chrom_thread(arg, arg[new_name],
                                               by_block=False))    
        p = plot(cget(thread, 1), style="lines", ymin=100,
                 ylog=10)

        # remove chrom
        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)
        arg.set_ancestral()
        arg.prune()

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

        k = 4
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


    def test_sample(self):

        k = 5
        n = 1e4
        rho = 1.5e-8 * 10
        mu = 2.5e-8 * 100
        length = 10000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        times = arghmm.get_time_points(ntimes=20)
        arghmm.discretize_arg(arg, times)

        # save
        arglib.write_arg("test/data/sample.arg", arg)
        fasta.write_fasta("test/data/sample.fa", seqs)

        tree = arg.get_marginal_tree(0)
        treelib.draw_tree_names(tree.get_tree(), minlen=5, scale=4e-4)


        new_name = "n%d" % (k-1)
        thread = list(arghmm.iter_chrom_thread(arg, arg[new_name],
                                               by_block=False))    
        p = plot(cget(thread, 1), style="lines", ymin=100,
                 ylog=10)

        # remove chrom
        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)
        arg.set_ancestral()
        arg.prune()

        model = arghmm.ArgHmm(arg, seqs, new_name="n%d" % (k-1), times=times,
                              rho=rho, mu=mu)
        print "states", len(model.states[0])
        print "muts", len(muts)
        print "recomb", len(model.recomb_pos) - 2, model.recomb_pos[1:-1]

        p.plot(model.recomb_pos, [10000] * len(model.recomb_pos),
               style="points")

        fw = probs_forward = arghmm.forward_algorithm(model, n, verbose=True)
        probs = arghmm.get_posterior_probs(model, length, verbose=True,
                                           probs_forward=fw)
        
        high = list(arghmm.iter_posterior_times(model, probs, .95))
        low = list(arghmm.iter_posterior_times(model, probs, .05))
        p.plot(high, style="lines")
        p.plot(low, style="lines")


        for i in xrange(1):
            path = arghmm.sample_posterior(model, length, verbose=True,
                                           probs_forward=fw)
            thread2 = [times[model.states[pos][state][1]]
                       for pos, state in enumerate(path)]

            p.gnuplot("set linestyle 3")
            p.plot(thread2, style="lines")

        pause()


    def test_forward_c(self):

        k = 3
        n = 1e4
        rho = 1.5e-8 * 3
        mu = 2.5e-8 * 100
        length = 100
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        arg.prune()
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        print arglib.get_recomb_pos(arg)
        print "muts", len(muts)
        print "recomb", len(arglib.get_recomb_pos(arg))

        times = arghmm.get_time_points(ntimes=30)
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
        probs1 = list(arghmm.forward_algorithm(model, length, verbose=True))
        util.toc()

        util.tic("python")
        probs2 = list(hmm.forward_algorithm(model, length, verbose=True))
        util.toc()

        print "probs1"
        pc(probs1)

        print "probs2"
        pc(probs2)
        

        for col1, col2 in izip(probs1, probs2):
            for a, b in izip(col1, col2):
                fequal(a, b, rel=.01)

        

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



#=============================================================================
if __name__ == "__main__":

    test_main()


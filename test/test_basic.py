
import unittest, random

import arghmm

from rasmus.common import *
from rasmus import stats, hmm
from rasmus.testing import *

from compbio import coal, arglib, fasta



def add_node(arg, node, time, pos, event):

    node2 = arg.new_node(event=event, age=time, children=[node], pos=pos)
    if event == "coal":
        node2.pos = 0

    parent = arg.get_local_parent(node, pos-.5)
    if parent:
        node.parents[node.parents.index(parent)] = node2
        parent.children[parent.children.index(node)] = node2
        node2.parents.append(parent)
    else:
        node.parents.append(node2)

    return node2


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


    def test_remove_cycles(self):

        arg = arglib.ARG(0, 10)

        a = arg.new_node(name="a", event="gene", age=0)
        b = arg.new_node(name="b", event="gene", age=0)
        r = arg.new_node(name="r", event="recomb", age=10, pos=1)
        d = arg.new_node(name="d", event="recomb", age=15, pos=4)
        c = arg.new_node(name="c", event="coal", age=20)
        t = arg.new_node(name="t", event="coal", age=30)
        t2 = arg.new_node(name="t2", event="coal", age=40)

        a.parents = [r]
        r.children = [a]

        r.parents = [c, d]
        c.children = [r, d]
        d.children = [r]
        d.parents = [c, t2]
        c.parents = [t]

        b.parents = [t]
        t.children = [c, b]
        
        t.parents = [t2]
        t2.children = [t, d]

        treelib.draw_tree_names(arg.get_marginal_tree(0).get_tree(),
                                maxlen=8, minlen=8)
        treelib.draw_tree_names(arg.get_marginal_tree(2).get_tree(),
                                maxlen=8, minlen=8)

        arg.set_ancestral()
        arglib.remove_self_cycles(arg)

        print
        print
        treelib.draw_tree_names(arg.get_marginal_tree(0).get_tree(),
                                maxlen=8, minlen=8)
        treelib.draw_tree_names(arg.get_marginal_tree(2).get_tree(),
                                maxlen=8, minlen=8)

        


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
        #arglib.write_arg("tmp/a.arg", arg)
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
        arg = arglib.smcify_arg(arg)
        print "recomb", arglib.get_recomb_pos(arg)

        model = arghmm.ArgHmm(arg, seqs, new_name="n%d" % k, times=times)

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

        # test transition switch matrix

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
        arg = arglib.smcify_arg(arg)
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
        mu = 2.5e-8 * 3
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
        rho = 1.5e-8 * 40
        mu = 2.5e-8 * 40
        length = 10000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        times = arghmm.get_time_points(ntimes=20)
        arghmm.discretize_arg(arg, times)

        # save
        #arglib.write_arg("test/data/k4.arg", arg)
        #fasta.write_fasta("test/data/k4.fa", seqs)

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


    #=========================================================================
    
    def test_forward_c(self):

        k = 5
        n = 1e4
        rho = 1.5e-8 * 3
        mu = 2.5e-8 * 100
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
        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)
        arg = arglib.smcify_arg(arg)

        model = arghmm.ArgHmm(arg, seqs, new_name="n%d" % (k-1), times=times,
                              rho=rho, mu=mu)
        print "states", len(model.states[0])
        
        util.tic("C")
        probs1 = list(arghmm.forward_algorithm(model, length, verbose=True))
        util.toc()

        util.tic("C2")
        probs1 = list(arghmm.forward_algorithm2(model, length, verbose=True))
        util.toc()

        util.tic("python")
        probs2 = list(hmm.forward_algorithm(model, length, verbose=True))
        util.toc()

        #print "probs1"
        #pc(probs1)

        #print "probs2"
        #pc(probs2)
        

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


    def test_trans_c(self):

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
        arg = arglib.smcify_arg(arg)
        print "recomb", arglib.get_recomb_pos(arg)

        model = arghmm.ArgHmm(arg, seqs, new_name="n%d" % k, times=times)

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



    #=========================================================================

    def test_sample(self):

        k = 5
        n = 1e4
        rho = 1.5e-8 * 10
        mu = 2.5e-8 * 100
        length = 10000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        #arg = arglib.read_arg("test/data/sample.arg")
        #seqs = fasta.read_fasta("test/data/sample.fa")

        times = arghmm.get_time_points(ntimes=20)
        arghmm.discretize_arg(arg, times)

        # save
        arglib.write_arg("test/data/sample.arg", arg)
        fasta.write_fasta("test/data/sample.fa", seqs)

        new_name = "n%d" % (k-1)
        thread = list(arghmm.iter_chrom_thread(arg, arg[new_name],
                                               by_block=False))    
        p = plot(cget(thread, 1), style="lines", ymin=100,
                 ylog=10)

        # remove chrom
        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)
        arg = arglib.smcify_arg(arg)

        arglib.write_arg("test/data/sample-prune.arg", arg)

        model = arghmm.ArgHmm(arg, seqs, new_name="n%d" % (k-1), times=times,
                              rho=rho, mu=mu)
        print "states", len(model.states[0])
        #print "muts", len(muts)
        print "recomb", len(model.recomb_pos) - 2, model.recomb_pos[1:-1]

        p.plot(model.recomb_pos, [10000] * len(model.recomb_pos),
               style="points")

        fw = probs_forward = arghmm.forward_algorithm(model, length,
                                                      verbose=True)
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
            p.plot(thread2, style="lines")
            util.write_list("test/data/sample.thread", path)

        pause()


    def test_sample2(self):

        k = 5
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        times = arghmm.get_time_points(ntimes=20)
        arghmm.discretize_arg(arg, times)

        # save
        arglib.write_arg("test/data/sample_recomb.arg", arg)
        fasta.write_fasta("test/data/sample_recomb.fa", seqs)
        #arg = arglib.read_arg("test/data/sample_recomb.arg")
        #seqs = fasta.read_fasta("test/data/sample_recomb.fa")
        

        # get new chrom
        new_name = "n%d" % (k-1)
        thread = list(arghmm.iter_chrom_thread(arg, arg[new_name],
                                               by_block=False))
        
        # remove chrom
        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)
        arg = arglib.smcify_arg(arg)
        

        # setup model
        model = arghmm.ArgHmm(arg, seqs, new_name="n%d" % (k-1), times=times,
                              rho=rho, mu=mu)
        print "states", len(model.states[0])
        print "muts", len(muts)
        print "recomb", len(model.recomb_pos) - 2, model.recomb_pos[1:-1]
        
        # sample a chrom thread
        fw = probs_forward = arghmm.forward_algorithm(model, length,
                                                      verbose=True)
        util.tic("sample thread")
        path = arghmm.sample_posterior(model, length, verbose=True,
                                       probs_forward=fw)
        util.toc()

        thread2 = list(arghmm.iter_thread_from_path(model, path))

        out = open("test/data/sample_recomb.trees", "w")
        for pos, (node_name, coal_time) in enumerate(thread2):
            tree = arg.get_marginal_tree(pos-.5)
            arglib.remove_single_lineages(tree)
            node = tree[node_name]
            node2 = add_node(tree, node, coal_time, -1, "coal")
            if not node2.parents:
                tree.root = node2
            leaf = tree.new_node(name=new_name, event="gene", age=0)
            leaf.parents.append(node2)
            node2.children.append(leaf)
            
            tree = tree.get_tree()
            out.write(str(pos)+"\t"+str(pos+1)+"\t")
            tree.write(out, oneline=True)
            out.write("\n")
        out.close()


    def test_sample_recomb(self):

        k = 2
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000

        rx = []
        ry = []

        for i in range(20):
            arg = arglib.sample_arg(k, n, rho, start=0, end=length)
            seqs = dict((l, "A"*length) for l in arg.leaf_names())
            times = arghmm.get_time_points(ntimes=20)
            arghmm.discretize_arg(arg, times)
            
            # count initial recomb count
            nrecombs = ilen(arghmm.iter_visible_recombs(arg))

            # get new chrom thread
            new_name = "n%d" % (k-1)
            thread_clades = list(arghmm.iter_chrom_thread(
                arg, arg[new_name], by_block=True, use_clades=True))
            #thread_clades = [arghmm.get_clade_point(arg, x[0], x[1], pos)
            #                 for pos, x in enumerate(thread)]

            # remove chrom
            util.tic("setup G_{n-1}")
            new_name = "n%d" % (k-1)
            keep = ["n%d" % i for i in range(k-1)]
            arglib.subarg_by_leaf_names(arg, keep)
            arg = arglib.smcify_arg(arg)
            nrecombs2 = ilen(arghmm.iter_visible_recombs(arg))
            new_recombs = nrecombs - nrecombs2
            util.toc()

            # setup model
            model = arghmm.ArgHmm(arg, seqs, new_name=new_name, times=times,
                                  rho=rho, mu=mu)
            thread = []
            for leaves, age, block in thread_clades:
                node = arghmm.arg_lca(arg, leaves, None, block[0]).name
                for i in range(block[0], block[1]):
                    thread.append((node, age))
            

            util.tic("sample recomb")
            for i in range(4):
                recombs = list(arghmm.sample_recombinations_thread(
                    model, thread))
                rx.append(new_recombs)
                ry.append(len(recombs))
            util.toc()

        p = plot(dither(rx, .2), dither(ry, .2),
                 xlab="actual new recombs", ylab="sampled new recombs")
        p.plot([0, max(rx)], [0, max(rx)], style="lines")
        
        pause()


    def test_add_thread(self):

        k = 5
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 40000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        arghmm.discretize_arg_recomb(arg)
        arg = arglib.smcify_arg(arg)
        arg.set_ancestral()
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        times = arghmm.get_time_points(ntimes=20)
        arghmm.discretize_arg(arg, times)

        # save
        arglib.write_arg("test/data/add_thread.arg", arg)
        fasta.write_fasta("test/data/add_thread.fa", seqs)
        #arg = arglib.read_arg("test/data/sample_recomb.arg")
        #seqs = fasta.read_fasta("test/data/sample_recomb.fa")

        # TODO: I could recode the thread using clade points...
        
        # get new chrom
        new_name = "n%d" % (k-1)
        thread = list(arghmm.iter_chrom_thread(arg, arg[new_name],
                                               by_block=False))
        p = plot(cget(thread, 1), style="lines", ymin=10, ylog=10)

        # remove chrom
        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)
        arg.set_ancestral()
        arglib.remove_self_cycles(arg)
        assert not arglib.has_self_cycles(arg)

        #arg = arglib.smcify_arg(arg)
        #arg.set_ancestral()
        #arg.prune()
        #arg.set_ancestral()

        #s = list(arglib.iter_self_cycles(arg))
        #while s:
        #    print "cycles", [x.name for x in s]
        #    arglib.remove_self_cycles(arg)
        #    arg.set_ancestral()
        #    s = list(arglib.iter_self_cycles(arg))
        

        # setup model
        model = arghmm.ArgHmm(arg, seqs, new_name="n%d" % (k-1), times=times,
                              rho=rho, mu=mu)
        print "states", len(model.states[0])
        print "muts", len(muts)
        print "recomb", len(model.recomb_pos) - 2, model.recomb_pos[1:-1]
        r = list(arghmm.iter_visible_recombs(arg))
        if len(r) > 0:
            p.plot([x.pos for x in r], [max(x.age,10) for x in r],
                   style="points")

        util.tic("sample recombs")
        recombs = list(arghmm.sample_recombinations_thread(model, thread))
        util.toc()
        r = [x for x in recombs if x[1] == new_name]
        if len(r) > 0:
            p.plot(cget(r, 0), [max(x[2], 10) for x in r], style="points")
        r = [x for x in recombs if x[1] != new_name]
        if len(r) > 0:
            p.plot(cget(r, 0), [max(x[2], 10) for x in r], style="points")


        arg3 = arglib.read_arg("test/data/add_thread.arg")
        arg = arghmm.add_arg_thread2(arg, new_name, thread, recombs,
                                     arg3=arg3)
        arglib.assert_arg(arg)

        
        # check thread
        thread2 = list(arghmm.iter_chrom_thread(arg, arg[new_name],
                                                by_block=False))
        p.plot(cget(thread2, 1), style="lines")        

        q = plot(dither(map(lambda x: log(clamp(x, 10, None)),
                            cget(thread, 1)), .1),
                 dither(map(lambda x: log(clamp(x, 10, None)),
                            cget(thread2, 1)), .1))
        
        pause()


    def test_add_thread_sample(self):
        """
        Test adding a sampled thread to an ARG
        """

        k = 5
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        arghmm.discretize_arg_recomb(arg)
        arg.set_ancestral()
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        times = arghmm.get_time_points(ntimes=20)
        arghmm.discretize_arg(arg, times)

        # save
        arglib.write_arg("test/data/sample_recomb.arg", arg)
        fasta.write_fasta("test/data/sample_recomb.fa", seqs)

        # get new chrom
        new_name = "n%d" % (k-1)
        thread = list(arghmm.iter_chrom_thread(arg, arg[new_name],
                                               by_block=False))    
        p = plot(cget(thread, 1), style="lines", ymin=8, ylog=10)
        
        # remove chrom
        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)
        arg = arglib.smcify_arg(arg)
        
        # setup model
        model = arghmm.ArgHmm(arg, seqs, new_name=new_name, times=times,
                              rho=rho, mu=mu)
        print "states", len(model.states[0])
        print "muts", len(muts)
        print "recomb", len(model.recomb_pos) - 2, model.recomb_pos[1:-1]
        r = list(arghmm.iter_visible_recombs(arg))
        if len(r) > 0:
            p.plot([x.pos for x in r], [max(x.age,10) for x in r],
                   style="points")

        # sample a chrom thread
        util.tic("forward algorithm")
        fw = probs_forward = arghmm.forward_algorithm(
            model, length, verbose=True)
        util.toc()
        util.tic("sample thread")
        path = arghmm.sample_posterior(model, length, verbose=True,
                                       probs_forward=fw)
        util.toc()

        thread2 = list(arghmm.iter_thread_from_path(model, path))
        p.plot(cget(thread2, 1), style="lines")


        util.tic("sample recombs")
        recombs = list(arghmm.sample_recombinations_thread(model, thread2))
        util.toc()
        r = [x for x in recombs if x[1] == new_name]
        if len(r) > 0:
            p.plot(cget(r, 0), [max(x[2], 10) for x in r], style="points")
        r = [x for x in recombs if x[1] != new_name]
        if len(r) > 0:
            p.plot(cget(r, 0), [max(x[2], 10) for x in r], style="points")


        #arg3 = arglib.read_arg("test/data/sample_recomb.arg")
        #arg = arghmm.add_arg_thread2(arg, new_name, thread2, recombs,
        #                             arg3=arg3)

        util.tic("add thread")
        arg = arghmm.add_arg_thread(arg, new_name, thread2, recombs)
        util.toc()
        
        arglib.write_arg("test/data/sample_recomb2.arg", arg)
        arglib.assert_arg(arg)

        # check thread
        thread3 = list(arghmm.iter_chrom_thread(arg, arg[new_name],
                                                by_block=False))
        p.plot(cget(thread3, 1), style="lines")        

        q = plot(dither(map(lambda x: log(clamp(x, 10, None)),
                            cget(thread2, 1)), .1),
                 dither(map(lambda x: log(clamp(x, 10, None)),
                            cget(thread3, 1)), .1))
        
        pause()


    def test_sample_arg(self):
        """
        Fully sample an ARG from stratch
        """

        k = 8
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        
        # save
        arglib.write_arg("test/data/sample_arg.arg", arg)
        fasta.write_fasta("test/data/sample_arg.fa", seqs)
        #arg = arglib.read_arg("test/data/sample_arg.arg")
        #seqs = fasta.read_fasta("test/data/sample_arg.fa")
        

        # get new chrom
        new_name = "n%d" % (k-1)
        thread = list(arghmm.iter_chrom_thread(arg, arg[new_name],
                                               by_block=False))    
        p = plot(cget(thread, 1), style="lines", ymin=8, ylog=10)


        def add_chrom(arg, new_name):
            util.tic("adding %s..." % new_name)

            times = arghmm.get_time_points(ntimes=20)
            model = arghmm.ArgHmm(arg, seqs, new_name=new_name,
                                  times=times, rho=rho, mu=mu)
            util.logger("states", len(model.states[0]))

            util.tic("forward algorithm")
            fw = probs_forward = arghmm.forward_algorithm2(
                model, length, verbose=True)
            util.toc()
            
            util.tic("sample thread")
            path = arghmm.sample_posterior(
                model, length, verbose=True, probs_forward=fw)
            util.toc()

            util.tic("sample recombs")
            thread2 = list(arghmm.iter_thread_from_path(model, path))
            recombs = list(arghmm.sample_recombinations_thread(model, thread2))
            util.toc()

            util.tic("add thread")
            arg = arghmm.add_arg_thread(arg, new_name, thread2, recombs)
            util.toc()

            util.toc()
            return arg, thread2

        util.tic("sample ARG")
        arg = arghmm.make_trunk_arg(arg.start, arg.end, name="n0")
        for j in xrange(1, k):
            new_name = "n%d" % j
            arg, thread2 = add_chrom(arg, new_name)
        util.toc()
        
        arglib.write_arg("test/data/sample_arg2.arg", arg)

        # check thread
        thread3 = list(arghmm.iter_chrom_thread(arg, arg[new_name],
                                                by_block=False))
        p.plot(cget(thread2, 1), style="lines")
        p.plot(cget(thread3, 1), style="lines")

        q = plot(dither(map(lambda x: log(clamp(x, 10, None)),
                            cget(thread2, 1)), .1),
                 dither(map(lambda x: log(clamp(x, 10, None)),
                            cget(thread3, 1)), .1))
        
        pause()



    def test_sample_arg2(self):
        """
        Fully sample an ARG from stratch
        """

        k = 5
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 1000
        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        
        # save
        arglib.write_arg("test/data/sample_arg.arg", arg)
        fasta.write_fasta("test/data/sample_arg.fa", seqs)
        #arg = arglib.read_arg("test/data/sample_arg.arg")
        #seqs = fasta.read_fasta("test/data/sample_arg.fa")
        

        # get new chrom
        new_name = "n%d" % (k-1)
        thread = list(arghmm.iter_chrom_thread(arg, arg[new_name],
                                               by_block=False))    
        p = plot(cget(thread, 1), style="lines", ymin=8, ylog=10)


        def add_chrom(arg, new_name):
            util.tic("adding %s..." % new_name)

            times = arghmm.get_time_points(ntimes=20)
            model = arghmm.ArgHmm(arg, seqs, new_name=new_name,
                                  times=times, rho=rho, mu=mu)
            util.logger("states", len(model.states[0]))

            util.tic("forward algorithm")
            fw = probs_forward = arghmm.forward_algorithm(
                model, length, verbose=True)
            util.toc()
            
            util.tic("sample thread")
            path = arghmm.sample_posterior(
                model, length, verbose=True, probs_forward=fw)
            util.toc()

            util.tic("sample recombs")
            thread2 = list(arghmm.iter_thread_from_path(model, path))
            recombs = list(arghmm.sample_recombinations_thread(model, thread2))
            util.toc()

            util.tic("add thread")
            arg = arghmm.add_arg_thread(arg, new_name, thread2, recombs)
            util.toc()

            util.toc()
            return arg, thread2

        for i in xrange(5):
            util.tic("sample ARG")
            arg = arghmm.make_trunk_arg(arg.start, arg.end, name="n0")
            for j in xrange(1, k):
                new_name = "n%d" % j
                arg, thread2 = add_chrom(arg, new_name)
            util.toc()
            p.plot(cget(thread2, 1), style="lines")
        
        pause()




#=============================================================================
if __name__ == "__main__":

    test_main()


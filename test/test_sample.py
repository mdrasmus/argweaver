

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


class Sample (unittest.TestCase):


    def test_sample_thread(self):

        k = 10
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
        arglib.write_arg("test/data/sample.arg", arg)
        fasta.write_fasta("test/data/sample.fa", seqs)

        new_name = "n%d" % (k-1)
        thread = list(arghmm.iter_chrom_thread(arg, arg[new_name],
                                               by_block=False))    
        p = plot(cget(thread, 1), style="lines", ymin=10,
                 ylog=10)

        # remove chrom
        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)
        arg = arglib.smcify_arg(arg)

        arglib.write_arg("test/data/sample-prune.arg", arg)

        model = arghmm.ArgHmm(arg, seqs, new_name="n%d" % (k-1), times=times,
                              rho=rho, mu=mu)
        print "states", len(model.states[0])
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
            #util.write_list("test/data/sample.thread", path)

        pause()


    def test_sample_thread2(self):

        k = 2
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000


        x = []
        y = []

        for i in range(1):
            arg = arglib.sample_arg(k, n, rho, start=0, end=length)
            muts = arglib.sample_arg_mutations(arg, mu)
            seqs = arglib.make_alignment(arg, muts)
            times = arghmm.get_time_points(ntimes=20)
            arghmm.discretize_arg(arg, times)

            new_name = "n%d" % (k-1)
            thread = list(arghmm.iter_chrom_thread(arg, arg[new_name],
                                                   by_block=False))

            # remove chrom
            keep = ["n%d" % i for i in range(k-1)]
            arglib.subarg_by_leaf_names(arg, keep)
            arg = arglib.smcify_arg(arg)
            model = arghmm.ArgHmm(arg, seqs, new_name="n%d" % (k-1),
                                  times=times, rho=rho, mu=mu)

            matrices = list(arghmm.iter_trans_emit_matrices(model, length))
            fw = probs_forward = arghmm.forward_algorithm(
                model, length, matrices=matrices, verbose=True)
            for i in xrange(20):
                path = arghmm.sample_posterior(model, length, verbose=True,
                                               matrices=matrices,
                                               probs_forward=fw)
                thread2 = list(arghmm.iter_thread_from_path(model, path))
                x.extend(cget(thread, 1)[::100])
                y.extend(cget(thread2, 1)[::100])

            arghmm.delete_trans_emit_matrices(matrices)

        x = map(safelog, x)
        y = map(safelog, y)
        p = plot(dither(x, .1), dither(y, .1), xmin=5, ymin=5)
        p.plot([1, max(x)], [1, max(x)], style="lines")

        pause()



    def test_sample_local_trees(self):
        """
        Write local trees
        """

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
        util.tic("sample thread")
        path = arghmm.sample_posterior(model, length, verbose=True)
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

        k = 5
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

        data = zip(rx, ry)
        write_delim("tmp/recomb", data)
        
        pause()


    def test_sample_recomb2(self):
        """
        Test the sampling of thread and recombinations
        """

        k = 5
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000

        rx = []
        ry = []
        
        for i in range(20):
            arg = arglib.sample_arg(k, n, rho, start=0, end=length)
            arghmm.discretize_arg_recomb(arg)
            arg = arglib.smcify_arg(arg)
            arg.set_ancestral()
            muts = arglib.sample_arg_mutations(arg, mu)
            seqs = arglib.make_alignment(arg, muts)
            times = arghmm.get_time_points(ntimes=20)
            arghmm.discretize_arg(arg, times)
            
            # count initial recomb count
            nrecombs = ilen(arghmm.iter_visible_recombs(arg))
            
            # get new chrom thread
            new_name = "n%d" % (k-1)
            
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

            util.tic("sample recomb")
            for j in range(2):
                path = arghmm.sample_posterior(model, length, verbose=False)
                thread = list(arghmm.iter_thread_from_path(model, path))
                recombs = list(arghmm.sample_recombinations_thread(
                    model, thread))
                rx.append(new_recombs)
                ry.append(len(recombs))

                if ry[-1] - rx[-1] > 40:
                    print thread[0:length:length//20]
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
        length = 20000
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
        thread_clades = list(arghmm.iter_chrom_thread(
            arg, arg[new_name], by_block=True, use_clades=True))
        p = plot(cget(thread, 1), style="lines", ymin=10, ylog=10)

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

        # remake thread from clades
        thread = []
        for leaves, age, block in thread_clades:
            node = arghmm.arg_lca(arg, leaves, None, block[0]).name
            for i in range(block[0], block[1]):
                thread.append((node, age))

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
        arg = arghmm.add_arg_thread2(arg, new_name, thread, recombs, arg3=arg3)

        
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
        nrecombs1 = len(r)

        # sample a chrom thread
        util.tic("sample thread")
        path = arghmm.sample_posterior(model, length, verbose=False)
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
        nrecombs_new = len(recombs)


        #arg3 = arglib.read_arg("test/data/sample_recomb.arg")
        #arg = arghmm.add_arg_thread2(arg, new_name, thread2, recombs,
        #                             arg3=arg3)

        util.tic("add thread")
        arg = arghmm.add_arg_thread(arg, new_name, thread2, recombs)
        util.toc()

        nrecombs2 = ilen(arghmm.iter_visible_recombs(arg))
        print "recombs", nrecombs1, nrecombs_new, nrecombs1 + nrecombs_new, nrecombs2
        
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

        k = 16
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

            util.tic("sample thread")
            path = arghmm.sample_posterior(model, length)
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
            
            util.tic("sample thread")
            path = arghmm.sample_posterior(
                model, length, verbose=True)
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


    def test_sample_arg3(self):
        """
        Fully sample an ARG from stratch
        """

        k = 5
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 2000
        times = arghmm.get_time_points(ntimes=20)

        def sample_arg(seqs):
            def add_chrom(arg, new_name):
                util.tic("adding %s..." % new_name)

                model = arghmm.ArgHmm(arg, seqs, new_name=new_name,
                                      times=times, rho=rho, mu=mu)
                util.logger("states", len(model.states[0]))

                util.tic("sample thread")
                path = arghmm.sample_posterior(model, length)
                util.toc()

                util.tic("sample recombs")
                thread = list(arghmm.iter_thread_from_path(model, path))
                recombs = list(arghmm.sample_recombinations_thread(
                    model, thread))
                util.toc()

                util.tic("add thread")
                arg = arghmm.add_arg_thread(arg, new_name, thread, recombs)
                util.toc()

                util.toc()
                return arg, thread

            names = seqs.keys()
            arg = arghmm.make_trunk_arg(0, length, name=names[0])
            for j in xrange(1, len(names)):
                arg, thread2 = add_chrom(arg, names[j])

            return arg


        arg = arglib.sample_arg(k, n, rho, start=0, end=length)
        arghmm.discretize_arg_recomb(arg)
        arg = arglib.smcify_arg(arg)
        arg.set_ancestral()
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        seqs.names.sort()

        util.tic()
        arg2 = sample_arg(seqs)
        util.toc()
        


    def test_sample_arg_recomb(self):
        """
        Fully sample an ARG from stratch
        """

        k = 4
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20)

        def sample_arg(seqs):
            def add_chrom(arg, new_name):
                util.tic("adding %s..." % new_name)

                util.tic("create model")
                model = arghmm.ArgHmm(arg, seqs, new_name=new_name,
                                      times=times, rho=rho, mu=mu)
                util.logger("states", len(model.states[0]))
                util.toc()

                util.tic("sample thread")
                path = arghmm.sample_posterior(model, length)
                util.toc()

                util.tic("sample recombs")
                thread = list(arghmm.iter_thread_from_path(model, path))
                recombs = list(arghmm.sample_recombinations_thread(
                    model, thread))
                util.toc()

                util.tic("add thread")
                arg = arghmm.add_arg_thread(arg, new_name, thread, recombs)
                util.toc()

                util.toc()
                return arg, thread

            names = seqs.keys()
            arg = arghmm.make_trunk_arg(0, length, name=names[0])
            for j in xrange(1, len(names)):
                arg, thread2 = add_chrom(arg, names[j])

            return arg


        rx = []
        ry = []
        util.tic("plot")
        for i in range(20):

            arg = arglib.sample_arg(k, n, rho, start=0, end=length)
            arghmm.discretize_arg_recomb(arg)
            arg = arglib.smcify_arg(arg)
            arg.set_ancestral()
            muts = arglib.sample_arg_mutations(arg, mu)
            seqs = arglib.make_alignment(arg, muts)
            
            nrecombs = ilen(arghmm.iter_visible_recombs(arg))

            for j in range(3):
                util.tic("sample ARG %d, %d" % (i, j))
                arg2 = sample_arg(seqs)
                util.toc()
                
                nrecombs2 = ilen(arghmm.iter_visible_recombs(arg2))
                rx.append(nrecombs)
                ry.append(nrecombs2)
        util.toc()

        print rx, ry


        p = plot(rx, ry,
                 xlab="true # recombinations",
                 ylab="inferred # recombinations")
        p.plot([min(rx), max(rx)], [min(rx), max(rx)], style="lines")
        
        
        pause()




#=============================================================================
if __name__ == "__main__":

    test_main()


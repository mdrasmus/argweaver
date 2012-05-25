

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


    def test_sim_dsmc(self):
        """
        Simulate from DSMC
        """

        k = 5
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 1000
        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length)
        #arg = arglib.sample_arg_smc(k, 2*n, rho, start=0, end=length)
        makedirs("test/data")
        arglib.write_arg("test/data/sim_dsmc.arg", arg)


    def test_sim_dsmc_cmp_recomb(self):
        """
        Simulate from DSMC and compare to SMC
        """

        k = 10
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000

        r1 = []; r2 = []
        for i in range(1, 100):
            print i
            arg = arglib.sample_arg_smc(k, 2*n, i/100. * rho,
                                         start=0, end=length)
            arg2 = arghmm.sample_arg_dsmc(k, 2*n, i/100. * rho,
                                         start=0, end=length)
            r1.append(ilen(j for j in arg if j.event == "recomb"))
            r2.append(ilen(j for j in arg2 if j.event == "recomb"))

        p = plot(r1, r2, xlab="# recomb SMC", ylab="# recomb DSMC")
        p.plot([min(r1), max(r1)], [min(r1), max(r1)], style="lines")
        
        pause()
        


    def test_sim_dsmc_cmp_arglen(self):
        """
        Simulate from DSMC and compare to SMC
        """

        k = 5
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 20000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)

        l1 = []; l2 = []
        for i in range(1, 100):
            print i
            tic('smc')
            arg = arglib.sample_arg_smc(k, i/100. * 2*n, rho, 
                                        start=0, end=length)
            toc()
            tic('dsmc')
            arg2 = arghmm.sample_arg_dsmc(k, i/100. * 2*n, rho,
                                          start=0, end=length, times=times)
            toc()
            l1.append(arglib.arglen(arg))
            l2.append(arglib.arglen(arg2))
            
        p2 = plot(l1, l2, xlab="length SMC", ylab="length DSMC")
        p2.plot([min(l1), max(l1)], [min(l1), max(l1)], style="lines")
        
        pause()


    def test_sim_dsmc_seq(self):
        """
        Simulate from DSMC and compare to SMC
        """

        k = 5
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)

        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, times=times,
                                      start=0, end=length)
        arg.set_ancestral()

        lk = []
        for i in range(1, 100):
            print i
            muts = arghmm.sample_arg_mutations(arg, mu, times=times)
            seqs = arglib.make_alignment(arg, muts)
            lk.append(arghmm.calc_likelihood(arg, seqs, mu=mu, times=times))

        lk.sort()
        plot(lk)
        #p = plot(r1, r2, xlab="# recomb SMC", ylab="# recomb DSMC")
        #p.plot([min(r1), max(r1)], [min(r1), max(r1)], style="lines")
        write_list("test/data/lk.txt", lk)
        
        pause()
        
        

    def test_sample_thread(self):
        """
        Sample a thread from the posterior of the ArgHmm
        """

        k = 10
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000
        arg = arglib.sample_arg(k, 2*n, rho, start=0, end=length)
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

        fw = arghmm.forward_algorithm(model, length, verbose=True)
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
        """
        Sample many threads from the ArgHmm
        """

        k = 5
        n = 1e4
        rho = 1.5e-8 * 40
        mu = 2.5e-8 * 40
        length = 10000
        times = arghmm.get_time_points(ntimes=20)

        x = []
        y = []

        for i in range(20):
            arg = arglib.sample_arg(k, 2*n, rho, start=0, end=length)
            muts = arglib.sample_arg_mutations(arg, mu)
            seqs = arglib.make_alignment(arg, muts)
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

            for i in xrange(1):
                path = arghmm.sample_posterior(model, length, verbose=True)
                thread2 = list(arghmm.iter_thread_from_path(model, path))
                x.extend(cget(thread, 1)[::100])
                y.extend(cget(thread2, 1)[::100])

        x = map(safelog, x)
        y = map(safelog, y)
        p = plot(dither(x, .1), dither(y, .1), xmin=5, ymin=5)
        p.plot([1, max(x)], [1, max(x)], style="lines")

        pause()


    def test_sample_recomb(self):
        """
        Sample recombinations for a true thread
        """

        k = 8
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)


        rx = []
        ry = []

        for i in range(40):
            #arg = arglib.sample_arg(k, 2*n, rho, start=0, end=length)
            #arg = arglib.sample_arg_smc(k, 2*n, rho, start=0, end=length)
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            seqs = dict((l, "A"*length) for l in arg.leaf_names())
            #arghmm.discretize_arg(arg, times)

            #trees = list(arglib.iter_marginal_trees(arg))
            #for i in range(1, len(trees)):
            #    print trees[i].root.age == trees[i-1].root.age
            
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
                                  popsize=n, rho=rho, mu=mu)
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

        print "avg ratio:", mean([safediv(i, j, 0) for i, j in zip(ry, rx)])

        
        pause()


    def test_sample_recomb2(self):
        """
        Sample both a thread and its recombinations
        """

        k = 6
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)

        rx = []
        ry = []
        
        for i in range(20):
            util.tic("sim arg")
            #arg = arglib.sample_arg(k, 2*n, rho, start=0, end=length)
            #arg = arglib.sample_arg_smc(k, 2*n, rho, start=0, end=length)
            #arghmm.discretize_arg_recomb(arg)
            #arg = arglib.smcify_arg(arg)
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times)
            #arg.set_ancestral()
            #muts = arglib.sample_arg_mutations(arg, mu)
            seqs = arglib.make_alignment(arg, muts)
            arghmm.discretize_arg(arg, times)
            util.toc()
            
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
                #if ry[-1] - rx[-1] > 40:
                #    print thread[0:length:length//20]

                #arg2 = arghmm.sample_thread(model, length)
                #recombs = ilen(x for x in arg2 if x.event == "recomb")
                #rx.append(new_recombs)
                #ry.append(recombs - nrecombs2)

            util.toc()

        p = plot(dither(rx, .2), dither(ry, .2),
                 xlab="actual new recombs", ylab="sampled new recombs")
        p.plot([0, max(rx)], [0, max(rx)], style="lines")

        print "avg ratio:", mean([safediv(i, j, 0) for i, j in zip(ry, rx)])
        
        pause()


    def test_add_thread(self):
        """
        Add a thread to an ARG
        """

        k = 5
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000
        arg = arglib.sample_arg(k, 2*n, rho, start=0, end=length)
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

        
        arg = arghmm.add_arg_thread(arg, new_name, thread, recombs)
        
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

        k = 3
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000
        arg = arglib.sample_arg(k, 2*n, rho, start=0, end=length)
        arghmm.discretize_arg_recomb(arg)
        arg.set_ancestral()
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        times = arghmm.get_time_points(ntimes=20)
        arghmm.discretize_arg(arg, times)

        # save
        #arglib.write_arg("test/data/sample_recomb.arg", arg)
        #fasta.write_fasta("test/data/sample_recomb.fa", seqs)

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


    def test_add_thread_sample_c(self):
        """
        Test adding a sampled thread to an ARG using C code
        """

        k = 3
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000
        arg = arglib.sample_arg(k, 2*n, rho, start=0, end=length)
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
        
        # sample a chrom thread
        util.tic("sample thread")        
        arg = arghmm.sample_thread(model, length)
        util.toc()



    def test_sample_arg(self):
        """
        Fully sample an ARG from stratch
        """

        k = 5
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000
        arg = arglib.sample_arg(k, 2*n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        
        # save
        arglib.write_arg("test/data/sample_arg.arg", arg)
        fasta.write_fasta("test/data/sample_arg.fa", seqs)       

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
        
        arglib.write_arg("test/data/sample_arg_out.arg", arg)
        
        pause()


    def test_sample_arg2(self):
        """
        Fully sample an ARG from stratch using API
        """

        k = 100
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 50000
        times = arghmm.get_time_points(ntimes=20)
        refine = 0
        
        arg = arglib.sample_arg(k, 2*n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)
        arg.write("test/data/sample_arg2.arg")
        seqs.write("test/data/sample_arg2.fa")

        seqs.names.sort()

        util.tic("sample ARG")
        arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu, times=times,
                                 refine=refine)
        util.toc()

        print ilen(x for x in arg2 if x.event == "recomb")

        #arg2.write("test/data/sample_arg2_out.arg")
        


    def test_sample_arg_recomb(self):
        """
        Plot the recombinations from a fully sampled ARG
        """

        k = 3
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        refine = 0

        print "times", times

        rx = []
        ry = []
        util.tic("plot")
        for i in range(20):
            #arg = arglib.sample_arg(k, 2*n, rho, start=0, end=length)
            #arg = arglib.sample_arg_smc(k, 2*n, rho, start=0, end=length)
            #arg.set_ancestral()
            #muts = arglib.sample_arg_mutations(arg, mu)

            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times=times)
            seqs = arglib.make_alignment(arg, muts)
            
            nrecombs = ilen(arghmm.iter_visible_recombs(arg))

            for j in range(3):
                util.tic("sample ARG %d, %d" % (i, j))
                arg2 = arghmm.sample_arg(seqs, rho=rho2, mu=mu, times=times,
                                         refine=refine)
                util.toc()
                
                nrecombs2 = ilen(arghmm.iter_visible_recombs(arg2))
                rx.append(nrecombs)
                ry.append(nrecombs2)
        util.toc()

        print "avg ratio:", mean([safediv(i, j, 0) for i, j in zip(ry, rx)])

        p = plot(rx, ry,
                 xlab="true # recombinations",
                 ylab="inferred # recombinations")
        p.plot([min(rx), max(rx)], [min(rx), max(rx)], style="lines")
        
        pause()


    def test_sample_arg_recomb2(self):
        """
        Plot the recombinations from a fully sampled ARG over many Gibb iters
        """

        k = 4
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)

        #arg = arglib.sample_arg_smc(k, 2*n, rho, start=0, end=length)
        #arg.set_ancestral()
        #muts = arglib.sample_arg_mutations(arg, mu)

        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg, mu, times=times)
        seqs = arglib.make_alignment(arg, muts)
            
        nrecombs = ilen(arghmm.iter_visible_recombs(arg))
        print "real # recombs", nrecombs

        y = []
        
        util.tic("sample ARG")
        arg2 = arghmm.sample_arg(seqs, rho=rho2, mu=mu, times=times)
        util.toc()
        
        nrecombs2 = ilen(arghmm.iter_visible_recombs(arg2))
        y.append(nrecombs2)

        for i in range(50):
            util.tic("resample ARG %d" % i)
            arg2 = arghmm.resample_arg(arg2, seqs, rho=rho, mu=mu, times=times)
            util.toc()
            nrecombs2 = ilen(arghmm.iter_visible_recombs(arg2))
            y.append(nrecombs2)
            print nrecombs2

        
        p = plot(y)
        makedirs("data/sample_arg_recomb2/")
        write_list("data/sample_arg_recomb2/recombs.txt", [nrecombs] + y)
        p.plot([0, len(y)], [nrecombs, nrecombs], style="lines")
        
        pause()



    def test_sample_arg_recomb_core(self):
        """
        Plot the recombinations from a fully sampled ARG over many Gibb iters
        """

        k = 8
        corek = 4
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000
        refine = 10
        times = arghmm.get_time_points(ntimes=20, maxtime=160000)

        arg = arglib.sample_arg_smc(k, 2*n, rho, start=0, end=length)
        arg.set_ancestral()
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)
            
        nrecombs = ilen(arghmm.iter_visible_recombs(arg))
        print "real # recombs", nrecombs

        y = []
        
        util.tic("sample core ARG")
        core_seqs = seqs.get(seqs.keys()[:corek])
        arg2 = arghmm.sample_arg(core_seqs, rho=rho, mu=mu, times=times,
                                 refine=refine)
        util.toc()
        
        nrecombs2 = ilen(arghmm.iter_visible_recombs(arg2))

        for i in range(50):
            util.tic("resample ARG %d" % i)
            arg2 = arghmm.resample_arg(arg2, seqs, rho=rho, mu=mu, times=times,
                                       refine=10)
            util.toc()
            nrecombs2 = ilen(arghmm.iter_visible_recombs(arg2))
            y.append(nrecombs2)
            print nrecombs2

        
        p = plot(y)
        p.plot([0, len(y)], [nrecombs, nrecombs], style="lines")
        
        
        pause()


    def test_sample_arg_arglen(self):
        """
        Plot the ARG length from a fully sampled ARG
        """

        k = 2
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        refine = 0

        print "times", times

        rx = []
        ry = []
        util.tic("plot")
        for i in range(100):
            arg = arglib.sample_arg_smc(k, 2*n, rho, start=0, end=length)
            arg.set_ancestral()
            muts = arglib.sample_arg_mutations(arg, mu)
            seqs = arglib.make_alignment(arg, muts)
            
            arglen = arglib.arglen(arg)

            for j in range(3):
                util.tic("sample ARG %d, %d" % (i, j))
                arg2 = arghmm.sample_arg(seqs, rho=rho2, mu=mu, times=times,
                                         refine=refine)
                util.toc()
                
                rx.append(arglen)
                ry.append(arglib.arglen(arg2))
        util.toc()

        print "avg ratio:", mean([safediv(i, j, 0) for i, j in zip(ry, rx)])

        p = plot(rx, ry,
                 xlab="true ARG length",
                 ylab="inferred ARG length")
        p.plot([min(rx), max(rx)], [min(rx), max(rx)], style="lines")
        
        pause()


    def test_sample_arg_arglen2(self):
        """
        Plot the recombinations from a fully sampled ARG over many Gibb iters
        """

        k = 8
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=160000)

        arg = arglib.sample_arg_smc(k, 2*n, rho, start=0, end=length)
        arg.set_ancestral()
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)
        #arghmm.discretize_arg(arg, times=times)
            
        arglen = arglib.arglen(arg)
        print "real # arglen %e" % arglen

        y = []
        
        util.tic("sample ARG")
        arg2 = arghmm.sample_arg(seqs, rho=rho2, mu=mu, times=times)
        util.toc()
        
        arglen2 = arglib.arglen(arg2)
        y.append(arglen2)

        for i in range(50):
            util.tic("resample ARG %d" % i)
            arg2 = arghmm.resample_arg(arg2, seqs, rho=rho, mu=mu, times=times,
                                       refine=3)
            util.toc()
            arglen2 = arglib.arglen(arg2)
            y.append(arglen2)
            print "%e" % arglen2

        
        p = plot(y)
        makedirs("data/sample_arg_arglen2/")
        write_list("data/sample_arg_arglen2/arglen.txt", [arglen] + y)
        p.plot([0, len(y)], [arglen, arglen], style="lines")
        
        
        pause()



    def test_sample_arg_lk(self):
        """
        Plot the ARG likelihood from a fully sampled ARG
        """

        k = 2
        n = 1e4
        rho = 1.5e-8 * .5
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        refine = 0

        print "times", times

        rx = []
        ry = []
        util.tic("plot")
        for i in range(20):
            arg = arglib.sample_arg_smc(k, 2*n, rho, start=0, end=length)
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            #arg.set_ancestral()
            #muts = arglib.sample_arg_mutations(arg, mu)
            muts = arghmm.sample_arg_mutations(arg, mu, times=times)
            seqs = arglib.make_alignment(arg, muts)
            #arghmm.discretize_arg(arg, times=times)
            
            lk = arghmm.calc_likelihood(arg, seqs, mu=mu, times=times)

            for j in range(3):
                util.tic("sample ARG %d, %d" % (i, j))
                arg2 = arghmm.sample_arg(seqs, rho=rho2, mu=mu, times=times,
                                         refine=refine)
                util.toc()

                lk2 = arghmm.calc_likelihood(arg2, seqs, mu=mu, times=times)
                rx.append(lk)
                ry.append(lk2)
        util.toc()

        print "avg ratio:", mean([safediv(i, j, 0) for i, j in zip(ry, rx)])

        p = plot(rx, ry,
                 xlab="true ARG likelihood",
                 ylab="inferred ARG likelihood")
        p.plot([min(rx), max(rx)], [min(rx), max(rx)], style="lines")
        
        pause()



    def test_sample_arg_lk2(self):
        """
        Plot the recombinations from a fully sampled ARG over many Gibb iters
        """
        k = 6
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)

        #arg = arglib.sample_arg_smc(k, 2*n, rho, start=0, end=length)
        #arg.set_ancestral()
        #muts = arglib.sample_arg_mutations(arg, mu)
        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg, mu, times=times)
        seqs = arglib.make_alignment(arg, muts)
        #arghmm.discretize_arg(arg, times)
            
        lk = arghmm.calc_likelihood(arg, seqs, mu=mu, times=times)
        print "real # lk", lk

        y = []
        
        util.tic("sample ARG")
        arg2 = arghmm.sample_arg(seqs, rho=rho2, mu=mu, times=times)
        util.toc()
        
        lk2 = arghmm.calc_likelihood(arg2, seqs, mu=mu, times=times)
        y.append(lk2)


        for i in range(200):
            util.tic("resample ARG %d" % i)
            arg2 = arghmm.resample_arg(arg2, seqs, rho=rho, mu=mu, times=times,
                                       refine=1)
            util.toc()
            lk2 = arghmm.calc_likelihood(arg2, seqs, mu=mu, times=times)
            y.append(lk2)
            print lk2

        
        p = plot(y)
        makedirs("data/sample_arg_lk2/")
        write_list("data/sample_arg_lk2/lk.txt", [lk] + y)
        p.plot([0, len(y)], [lk, lk], style="lines")
        
        
        pause()


    def test_sample_arg_joint(self):
        """
        Plot the ARG joint prob from a fully sampled ARG
        """

        k = 6
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        refine = 0

        print "times", times

        rx = []
        ry = []
        util.tic("plot")
        for i in range(20):
            #arg = arglib.sample_arg_smc(k, 2*n, rho, start=0, end=length)
            #arg.set_ancestral()
            #muts = arglib.sample_arg_mutations(arg, mu)
            #arghmm.discretize_arg(arg, times=times)
            
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times=times)
            seqs = arglib.make_alignment(arg, muts)
            
            lk = arghmm.calc_joint_prob(arg, seqs, mu=mu, rho=rho, times=times)

            for j in range(3):
                util.tic("sample ARG %d, %d" % (i, j))
                arg2 = arghmm.sample_arg(seqs, rho=rho2, mu=mu, times=times,
                                         refine=refine)
                util.toc()

                lk2 = arghmm.calc_joint_prob(arg2, seqs, mu=mu, rho=rho,
                                             times=times)
                rx.append(lk)
                ry.append(lk2)
        util.toc()

        print "avg ratio:", mean([safediv(i, j, 0) for i, j in zip(ry, rx)])

        p = plot(rx, ry,
                 xlab="true ARG joint probability",
                 ylab="inferred ARG joint probability")
        p.plot([min(rx), max(rx)], [min(rx), max(rx)], style="lines")
        
        pause()



    def test_sample_arg_joint2(self):
        """
        Plot the recombinations from a fully sampled ARG over many Gibb iters
        """
        k = 12
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=160000)

        #arg = arglib.sample_arg_smc(k, 2*n, rho, start=0, end=length)
        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length)
        arg.set_ancestral()
        muts = arghmm.sample_arg_mutations(arg, mu, times=times)
        #muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)
        #arghmm.discretize_arg(arg, times)
        
        
        lk = arghmm.calc_joint_prob(arg, seqs, mu=mu, rho=rho, times=times)
        print "real joint", lk

        y = []
        
        util.tic("sample ARG")
        arg2 = arghmm.sample_arg(seqs, rho=rho2, mu=mu, times=times)
        util.toc()
        
        lk2 = arghmm.calc_joint_prob(arg2, seqs, mu=mu, rho=rho, times=times)
        y.append(lk2)

        for i in range(50):
            util.tic("resample ARG %d" % i)
            arg2 = arghmm.resample_arg(arg2, seqs, rho=rho, mu=mu, times=times,
                                       refine=1)
            util.toc()
            lk2 = arghmm.calc_joint_prob(arg2, seqs, mu=mu, rho=rho, times=times)
            y.append(lk2)
            print lk2

        
        p = plot(y)
        makedirs("data/sample_arg_joint2/")
        write_list("data/sample_arg_joint2/joint.txt", [lk] + y)
        p.plot([0, len(y)], [lk, lk], style="lines")
        
        
        pause()


    def test_treeset(self):
        """
        Test the treeset representation of an ARG
        """

        k = 5
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000
        arg = arglib.sample_arg(k, 2*n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)

        times = arghmm.get_time_points(ntimes=20)
        arghmm.discretize_arg(arg, times)

        # remove chrom
        new_name = "n%d" % (k-1)
        keep = ["n%d" % i for i in range(k-1)]
        arglib.subarg_by_leaf_names(arg, keep)
        arg = arglib.smcify_arg(arg)
        print list(x.pos for x in arg if x.event == "recomb")

        
        print arghmm.get_treeset(arg, times)



#=============================================================================
if __name__ == "__main__":

    test_main()


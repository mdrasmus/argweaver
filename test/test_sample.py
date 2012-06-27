

import unittest, random

import arghmm

from rasmus.common import *
from rasmus import stats, hmm
from rasmus.testing import *

from compbio import coal, arglib, fasta



def strip_tree(tree):
    arglib.remove_single_lineages(tree)
    return tree

def hash_tree(tree):
    return phylo.hash_tree(tree)

def get_tree(arg, pos):
    return strip_tree(arg.get_marginal_tree(pos)).get_tree()

def get_branch_correct(arg1, arg2, pos):
    tree1 = get_tree(arg1, pos)
    tree2 = get_tree(arg2, pos)
    return 1.0 - phylo.robinson_foulds_error(tree1, tree2)


#=============================================================================

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
        makedirs("test/data/sim_dsmc")
        arglib.write_arg("test/data/sim_dsmc/0.arg", arg)


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

        for i in range(10):
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times)
            seqs = arglib.make_alignment(arg, muts)

            new_name = "n%d" % (k-1)
            thread = list(arghmm.iter_chrom_thread(arg, arg[new_name],
                                                   by_block=False))

            # remove chrom
            arg = arghmm.remove_arg_thread(arg, new_name)

            for j in xrange(1):
                arg2 = arghmm.resample_arg(arg, seqs, popsize=n,
                                           times=times, rho=rho, mu=mu,
                                           refine=0)
                thread2 = list(arghmm.iter_chrom_thread(
                    arg2, arg2[new_name], by_block=False))
                x.extend(cget(thread, 1)[::100])
                y.extend(cget(thread2, 1)[::100])
                print len(x), len(y)

        x = map(safelog, x)
        y = map(safelog, y)
        p = plot(dither(x, .1), dither(y, .1), xmin=5, ymin=5)
        p.plot([1, max(x)], [1, max(x)], style="lines")

        pause()


    def test_max_thread(self):
        """
        Maximize a thread from the posterior of the ArgHmm
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

        new_name = "n%d" % (k-1)
        thread = list(arghmm.iter_chrom_thread(arg, arg[new_name],
                                               by_block=False))    
        p = plot(cget(thread, 1), style="lines", ymin=10,
                 ylog=10)

        # remove chrom
        arg = arghmm.remove_arg_thread(arg, new_name)

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


        arg2 = arghmm.max_thread(arg, seqs, rho=rho, mu=mu,
                                 popsize=n, times=times)
        thread2 = list(arghmm.iter_chrom_thread(arg2, arg2[new_name],
                                               by_block=False))
        p.plot(cget(thread2, 1), style="lines")
        #util.write_list("test/data/sample.thread", path)

        pause()


    def _test_sample_recomb(self):
        """
        Sample recombinations for a true thread
        """

        k = 2
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)


        rx = []
        ry = []

        for i in range(40):
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            seqs = dict((l, "A"*length) for l in arg.leaf_names())
            
            # count initial recomb count
            nrecombs = ilen(arghmm.iter_visible_recombs(arg))

            # get new chrom thread
            new_name = "n%d" % (k-1)
            thread_clades = list(arghmm.iter_chrom_thread(
                arg, arg[new_name], by_block=True, use_clades=True))

            # remove chrom
            util.tic("setup G_{n-1}")
            new_name = "n%d" % (k-1)
            arg = arghmm.remove_arg_thread(arg, new_name)
            nrecombs2 = ilen(arghmm.iter_visible_recombs(arg))
            new_recombs = nrecombs - nrecombs2
            util.toc()

            # setup model and convert thread
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
        print "avg ratio:", mean([safediv(i, j, 0) for i, j in zip(ry, rx)])

        
        pause()


    def test_sample_recomb2(self):
        """
        Sample both a thread and its recombinations
        """

        k = 5
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)

        rx = []
        ry = []
        
        for i in range(20):
            util.tic("sim arg")
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times)
            seqs = arglib.make_alignment(arg, muts)
            util.toc()
            
            # count initial recomb count
            nrecombs = ilen(arghmm.iter_visible_recombs(arg))
            
            # get new chrom thread
            new_name = "n%d" % (k-1)
            
            # remove chrom
            util.tic("setup G_{n-1}")
            new_name = "n%d" % (k-1)
            arg = arghmm.remove_arg_thread(arg, new_name)
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
        arg = arghmm.remove_arg_thread(arg, new_name)
        
        # setup model
        model = arghmm.ArgHmm(arg, seqs, new_name=new_name, times=times,
                              rho=rho, mu=mu)
        print "states", len(model.states[0])
        print "muts", len(muts)
        print "recomb", len(model.recomb_pos) - 2, model.recomb_pos[1:-1]
        
        # sample a chrom thread
        util.tic("sample thread")        
        arg = arghmm.sample_thread(arg, seqs, rho=rho, mu=mu,
                                   popsize=n, times=times)
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
        length = 500000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        refine = 0
        
        arg = arglib.sample_arg(k, 2*n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)
        #arg.write("test/data/sample_arg2.arg")
        #seqs.write("test/data/sample_arg2.fa")

        print times
        
        util.tic("sample ARG")
        arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu, times=times,
                                 refine=refine, verbose=True)
        util.toc()

        print ilen(x for x in arg2 if x.event == "recomb")

        #arg2.write("test/data/sample_arg2_out.arg")


    def test_sample_arg_region(self):
        """
        Resample a region within an ARG
        """

        k = 10
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20)
        refine = 0
        region = [1000, 1000+500]
        
        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg, mu, times=times)
        seqs = arglib.make_alignment(arg, muts)
        print ilen(x for x in arg if x.event == "recomb")
        
        util.tic("sample ARG")
        arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu, times=times,
                                 refine=refine, verbose=True)
        util.toc()
        print ilen(x for x in arg2 if x.event == "recomb")

        util.tic("sample ARG region")
        arg2 = arghmm.resample_arg_region(arg2, seqs, region[0], region[1],
                                          rho=rho, mu=mu, times=times,
                                          verbose=True)
        util.toc()
        print ilen(x for x in arg2 if x.event == "recomb")
        #arg2.write("test/data/sample_arg2_out.arg")


    def test_sample_arg_region2(self):
        """
        Resample tangled regions in an ARG
        """

        k = 12
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20)
        refine = 0

        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg, mu, times=times)
        seqs = arglib.make_alignment(arg, muts)
        print ilen(x for x in arg if x.event == "recomb")
        recomb_pos = list(x.pos for x in arg if x.event == "recomb")
        
        util.tic("sample ARG")
        arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu, times=times,
                                 refine=refine, verbose=True)
        util.toc()
        print ilen(x for x in arg2 if x.event == "recomb")
        recomb_pos2 = list(x.pos for x in arg2 if x.event == "recomb")
        recomb_pos3 = recomb_pos2

        for i in range(40):
            maxr = 0
            for i,j,a,b in iter_window_index(recomb_pos3, 500):
                r = j - i + 1
                if r > maxr:
                    maxr = r
                    region = [recomb_pos3[i]-10, recomb_pos3[j]+10]
            print i, region

            util.tic("sample ARG region")
            arg2 = arghmm.resample_arg_region(arg2, seqs, region[0], region[1],
                                              rho=rho, mu=mu, times=times,
                                              verbose=True)
            util.toc()
            print ilen(x for x in arg2 if x.event == "recomb")
            recomb_pos3 = list(x.pos for x in arg2 if x.event == "recomb")


        # plotting
        p = Gnuplot()
        p.enableOutput(False)
        p.plot(recomb_pos, [0] * len(recomb_pos))
        p.plot(recomb_pos2, [-10] * len(recomb_pos2))
        p.plot(recomb_pos3, [-20] * len(recomb_pos3))
        
        x = []; y = []
        for i,j,a,b in iter_window_index(recomb_pos, 1000):
            x.append((a+b)/2.0); y.append(j - i + 1)
        p.plot(x, y, style="lines")

        x = []; y = []
        for i,j,a,b in iter_window_index(recomb_pos2, 1000):
            x.append((a+b)/2.0); y.append(j - i + 1)
        p.plot(x, y, style="lines")

        x = []; y = []
        for i,j,a,b in iter_window_index(recomb_pos3, 1000):
            x.append((a+b)/2.0); y.append(j - i + 1)
        p.plot(x, y, style="lines")

        p.enableOutput(True)
        p.replot()
        

        pause()
        

    def test_sample_arg_recomb(self):
        """
        Plot the recombinations from a fully sampled ARG
        """

        k = 2
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 20000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        nremove = 1
        refine = 0

        print "times", times

        rx = []
        ry = []
        util.tic("plot")
        for i in range(40):
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times=times)
            seqs = arglib.make_alignment(arg, muts)
            
            nrecombs = ilen(arghmm.iter_visible_recombs(arg))

            for j in range(1):
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


    def test_sample_arg_recomb_core(self):
        """
        Plot the ARG joint prob from a fully sampled ARG
        """

        k = 12
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        refine = 0; nremove = 1; core=5

        names = []
        rx = []
        ry = []
        util.tic("plot")
        for i in range(20):
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times=times)
            seqs = arglib.make_alignment(arg, muts)
            r = ilen(arghmm.iter_visible_recombs(arg))

            for j in range(1):
                util.tic("sample ARG %d, %d" % (i, j))
                
                names = seqs.keys()
                random.shuffle(names)
                
                arg2 = arglib.subarg_by_leaf_names(arg, names[:core])
                arg2 = arglib.smcify_arg(arg2)
                arg2 = arghmm.resample_arg(arg2, seqs, rho=rho, mu=mu,
                                           times=times, refine=2)
                
                util.toc()

                r2 = ilen(arghmm.iter_visible_recombs(arg2))
                rx.append(r)
                ry.append(r2)
                names.append([i, j])
        util.toc()

        print "avg ratio:", mean([safediv(i, j, 0) for i, j in zip(ry, rx)])

        p = plot(rx, ry,
                 xlab="true ARG recombs",
                 ylab="inferred ARG recombs")
        p.plot([min(rx), max(rx)], [min(rx), max(rx)], style="lines")
        
        pause()


    def test_sample_arg_recomb_region(self):
        """
        Plot the ARG joint prob from a fully sampled ARG
        """

        k = 12
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        refine = 0; nremove = 1;

        names = []
        rx = []
        ry = []
        util.tic("plot")
        for i in range(20):
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times=times)
            seqs = arglib.make_alignment(arg, muts)
            r = ilen(arghmm.iter_visible_recombs(arg))

            for j in range(1):
                util.tic("sample ARG %d, %d" % (i, j))
                arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu,
                                         times=times, refine=refine)
                arg2 = arghmm.resample_arg_regions(
                    arg2, seqs, niters=20, width=1000,
                    rho=rho, mu=mu,popsize=n, times=times, verbose=True)
                util.toc()

                r2 = ilen(arghmm.iter_visible_recombs(arg2))
                rx.append(r)
                ry.append(r2)
                names.append([i, j])
        util.toc()

        print "avg ratio:", mean([safediv(i, j, 0) for i, j in zip(ry, rx)])

        p = plot(rx, ry,
                 xlab="true ARG recombs",
                 ylab="inferred ARG recombs")
        p.plot([min(rx), max(rx)], [min(rx), max(rx)], style="lines")
        
        pause()


    def test_sample_arg_recomb2(self):
        """
        Plot the recombinations from a fully sampled ARG over many Gibb iters
        """

        k = 12
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 20000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        write = False
        #nremove = 2; refine = 5
        nremove = 1; refine = 1

        makedirs("test/data/sample_arg_recomb2/")

        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg, mu, times=times)
        seqs = arglib.make_alignment(arg, muts)
        if write:
            arglib.write_arg("test/data/sample_arg_recomb2/arg.arg", arg)
            seqs.write("test/data/sample_arg_recomb2/seqs.fa")
            
        nrecombs = ilen(arghmm.iter_visible_recombs(arg))
        print "real # recombs", nrecombs

        y = []
        
        util.tic("sample ARG")
        arg2 = arghmm.sample_arg(seqs, rho=rho2, mu=mu, times=times)
        util.toc()
        
        nrecombs2 = ilen(arghmm.iter_visible_recombs(arg2))
        y.append(nrecombs2)
        print nrecombs2

        for i in range(50):
            util.tic("resample ARG %d" % i)
            arg2 = arghmm.resample_arg(arg2, seqs, rho=rho, mu=mu, times=times,
                                       refine=refine, nremove=nremove)
            util.toc()
            nrecombs2 = ilen(arghmm.iter_visible_recombs(arg2))
            y.append(nrecombs2)
            print nrecombs2

            if write:
                arglib.write_arg("test/data/sample_arg_recomb2/%d.arg" % i,
                                 arg2)

        
        p = plot(y)
        p.plot([0, len(y)], [nrecombs, nrecombs], style="lines")
        
        pause()


    def test_max_arg_recomb2(self):
        """
        Plot the recombinations from a fully sampled ARG over many Gibb iters
        """

        k = 5
        n = 1e4
        rho = 1.5e-8 * 40
        rho2 = rho
        mu = 2.5e-8 * 40
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        write = False
        #nremove = 2; refine = 5
        nremove=1; refine = 1

        makedirs("test/data/sample_arg_recomb2/")

        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg, mu, times=times)
        seqs = arglib.make_alignment(arg, muts)
        if write:
            arglib.write_arg("test/data/sample_arg_recomb2/arg.arg", arg)
            seqs.write("test/data/sample_arg_recomb2/seqs.fa")
            
        nrecombs = ilen(arghmm.iter_visible_recombs(arg))
        print "real # recombs", nrecombs

        y = []
        
        util.tic("sample ARG")
        arg2 = arghmm.sample_arg(seqs, rho=rho2, mu=mu, times=times)
        util.toc()
        
        nrecombs2 = ilen(arghmm.iter_visible_recombs(arg2))
        y.append(nrecombs2)
        print nrecombs2

        for i in range(40):
            util.tic("resample ARG %d" % i)
            arg2 = arghmm.remax_arg(arg2, seqs, rho=rho, mu=mu, times=times,
                                    refine=refine, nremove=nremove)
            util.toc()
            nrecombs2 = ilen(arghmm.iter_visible_recombs(arg2))
            y.append(nrecombs2)
            print nrecombs2

            if write:
                arglib.write_arg("test/data/sample_arg_recomb2/%d.arg" % i,
                                 arg2)

        
        p = plot(y)
        p.plot([0, len(y)], [nrecombs, nrecombs], style="lines")
        
        pause()


    def test_sample_arg_arglen(self):
        """
        Plot the ARG length from a fully sampled ARG
        """

        k = 12
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
        for i in range(50):
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times=times)
            seqs = arglib.make_alignment(arg, muts)
            
            arglen = arglib.arglen(arg)

            for j in range(1):
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



    def test_max_arg_arglen2(self):
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
        print "real arglen %e" % arglen

        y = []
        
        util.tic("sample ARG")
        arg2 = arghmm.sample_arg(seqs, rho=rho2, mu=mu, times=times)
        util.toc()
        
        arglen2 = arglib.arglen(arg2)
        y.append(arglen2)

        for i in range(20):
            util.tic("resample ARG %d" % i)
            arg2 = arghmm.remax_arg(arg2, seqs, rho=rho, mu=mu, times=times,
                                    refine=1)
            util.toc()
            arglen2 = arglib.arglen(arg2)
            y.append(arglen2)
            print "%e" % arglen2

        
        p = plot(y)
        makedirs("data/max_arg_arglen2/")
        write_list("data/max_arg_arglen2/arglen.txt", [arglen] + y)
        p.plot([0, len(y)], [arglen, arglen], style="lines")
        
        pause()

        

    def test_sample_arg_lk(self):
        """
        Plot the ARG likelihood from a fully sampled ARG
        """

        k = 2
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        refine = 0
        write = True
        if write:
            make_clean_dir("test/data/sample_arg_lk")

        print "times", times

        names = []
        rx = []
        ry = []
        util.tic("plot")
        for i in range(20):
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times=times)
            seqs = arglib.make_alignment(arg, muts)
            if write:
                arglib.write_arg("test/data/sample_arg_lk/%d.arg" % i, arg)
                seqs.write("test/data/sample_arg_lk/%d.fa" % i)
            
            lk = arghmm.calc_likelihood(arg, seqs, mu=mu, times=times)

            for j in range(4):
                util.tic("sample ARG %d, %d" % (i, j))
                arg2 = arghmm.sample_arg(seqs, rho=rho2, mu=mu, times=times,
                                         refine=refine)
                if write:
                    arglib.write_arg("test/data/sample_arg_lk/%d-%d.arg" %
                                     (i, j), arg2)
                util.toc()

                lk2 = arghmm.calc_likelihood(arg2, seqs, mu=mu, times=times)
                names.append([i, j])
                rx.append(lk)
                ry.append(lk2)
        util.toc()

        if write:
            data = [[i, j, a, b] for (i, j), a, b in zip(names, rx, ry)]
            write_delim("test/data/sample_arg_lk/data.txt", data)

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
        k = 8
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


        for i in range(50):
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

        k = 5
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        refine = 0; nremove = 1
        write = False
        if write:
            make_clean_dir("test/data/sample_arg_joint")

        names = []
        rx = []
        ry = []
        util.tic("plot")
        for i in range(20):
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times=times)
            seqs = arglib.make_alignment(arg, muts)
            if write:
                arglib.write_arg("test/data/sample_arg_joint/%d.arg" % i, arg)
                seqs.write("test/data/sample_arg_joint/%d.fa" % i)

            
            lk = arghmm.calc_joint_prob(arg, seqs, mu=mu, rho=rho, times=times)

            for j in range(1):
                util.tic("sample ARG %d, %d" % (i, j))
                arg2 = arghmm.sample_arg(seqs, rho=rho2, mu=mu, times=times,
                                         refine=refine, nremove=nremove)
                util.toc()

                lk2 = arghmm.calc_joint_prob(arg2, seqs, mu=mu, rho=rho,
                                             times=times)
                rx.append(lk)
                ry.append(lk2)
                names.append([i, j])
                if write:
                    arglib.write_arg("test/data/sample_arg_joint/%d-%d.arg" %
                                     (i, j), arg2)
        util.toc()

        if write:
            data = [[i, j, a, b] for (i, j), a, b in zip(names, rx, ry)]
            write_delim("test/data/sample_arg_joint/data.txt", data)
        

        print "avg ratio:", mean([safediv(i, j, 0) for i, j in zip(ry, rx)])

        p = plot(rx, ry,
                 xlab="true ARG joint probability",
                 ylab="inferred ARG joint probability")
        p.plot([min(rx), max(rx)], [min(rx), max(rx)], style="lines")
        
        pause()


    def test_sample_arg_joint_core(self):
        """
        Plot the ARG joint prob from a fully sampled ARG
        """

        k = 12
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        refine = 5; nremove = 1; core=5
        write = False
        if write:
            make_clean_dir("test/data/sample_arg_joint")


        names = []
        rx = []
        ry = []
        util.tic("plot")
        for i in range(50):
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times=times)
            seqs = arglib.make_alignment(arg, muts)
            if write:
                arglib.write_arg("test/data/sample_arg_joint/%d.arg" % i, arg)
                seqs.write("test/data/sample_arg_joint/%d.fa" % i)

            
            lk = arghmm.calc_joint_prob(arg, seqs, mu=mu, rho=rho, times=times)

            for j in range(1):
                util.tic("sample ARG %d, %d" % (i, j))
                
                names = seqs.keys()
                random.shuffle(names)
                
                #arg2 = arglib.subarg_by_leaf_names(arg, names[:core])
                #arg2 = arglib.smcify_arg(arg2)

                # infer core
                arg2 = arghmm.sample_arg(seqs.get(names[:core]),
                                         rho=rho, mu=mu, times=times,
                                         refine=refine)
                
                arg2 = arghmm.resample_arg(arg2, seqs, rho=rho, mu=mu,
                                           times=times,
                                           refine=0)
                
                util.toc()

                lk2 = arghmm.calc_joint_prob(arg2, seqs, mu=mu, rho=rho,
                                             times=times)
                rx.append(lk)
                ry.append(lk2)
                names.append([i, j])
                if write:
                    arglib.write_arg("test/data/sample_arg_joint/%d-%d.arg" %
                                     (i, j), arg2)
        util.toc()

        if write:
            data = [[i, j, a, b] for (i, j), a, b in zip(names, rx, ry)]
            write_delim("test/data/sample_arg_joint/data.txt", data)
        

        print "avg ratio:", mean([safediv(i, j, 0) for i, j in zip(ry, rx)])

        p = plot(rx, ry,
                 xlab="true ARG joint probability",
                 ylab="inferred ARG joint probability")
        p.plot([min(rx), max(rx)], [min(rx), max(rx)], style="lines")
        
        pause()


    def test_sample_arg_joint_region(self):
        """
        Plot the ARG joint prob from a fully sampled ARG
        """

        k = 12
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        refine = 0; nremove = 1;

        names = []
        rx = []
        ry = []
        util.tic("plot")
        for i in range(20):
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times=times)
            seqs = arglib.make_alignment(arg, muts)
            lk = arghmm.calc_joint_prob(arg, seqs, mu=mu, rho=rho, times=times)

            for j in range(1):
                util.tic("sample ARG %d, %d" % (i, j))
                arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu,
                                         times=times, refine=refine)
                arg2 = arghmm.resample_arg_regions(
                    arg2, seqs, niters=40, width=500,
                    rho=rho, mu=mu,popsize=n, times=times, verbose=True)
                util.toc()

                lk2 = arghmm.calc_joint_prob(arg2, seqs, mu=mu, rho=rho,
                                             times=times)
                
                rx.append(lk)
                ry.append(lk2)
                names.append([i, j])
        util.toc()

        print "avg ratio:", mean([safediv(i, j, 0) for i, j in zip(ry, rx)])

        p = plot(rx, ry,
                 xlab="true ARG joint probability",
                 ylab="inferred ARG joint probability")
        p.plot([min(rx), max(rx)], [min(rx), max(rx)], style="lines")
        p.save("test/data/sample_arg_joint_region.pdf")
        
        pause()



    def test_sample_arg_joint2(self):
        """
        Plot the recombinations from a fully sampled ARG over many Gibb iters
        """
        k = 7
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        nremove = 1; refine = 5
        write = False
        
        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg, mu, times=times)
        seqs = arglib.make_alignment(arg, muts)
        if write:
            arglib.write_arg("test/data/sample_arg_joint/%d.arg" % i, arg)
            seqs.write("test/data/sample_arg_joint/%d.fa" % i)
        
        lk = arghmm.calc_joint_prob(arg, seqs, mu=mu, rho=rho, times=times,
                                    popsize=n)
        print "real joint", lk

        y = []
        
        util.tic("sample ARG")
        arg2 = arghmm.sample_arg(seqs, rho=rho2, mu=mu, times=times,
                                 popsize=n)
        util.toc()
        
        lk2 = arghmm.calc_joint_prob(arg2, seqs, mu=mu, rho=rho, times=times,
                                     popsize=n)
        y.append(lk2)

        for i in range(50):
            util.tic("resample ARG %d" % i)
            arg2 = arghmm.resample_arg(arg2, seqs, rho=rho, mu=mu, times=times,
                                     popsize=n, refine=refine, nremove=nremove)
            util.toc()
            lk2 = arghmm.calc_joint_prob(arg2, seqs, mu=mu, rho=rho,
                                         times=times, popsize=n)
            y.append(lk2)
            print lk2

        
        p = plot(y)
        makedirs("data/sample_arg_joint2/")
        write_list("data/sample_arg_joint2/joint.txt", [lk] + y)
        p.plot([0, len(y)], [lk, lk], style="lines")
        
        
        pause()


    def test_sample_arg_joint3(self):
        """
        Plot the recombinations from a fully sampled ARG over many Gibb iters
        """
        k = 8
        n = 1e4
        rho = 1.5e-8 * 40
        rho2 = rho
        mu = 2.5e-8 * 40
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        nremove = 1; refine = 1
        #nremove = 1; refine = 1
        write = False
        
        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg, mu, times=times)
        seqs = arglib.make_alignment(arg, muts)
        if write:
            arglib.write_arg("test/data/sample_arg_joint3/%d.arg" % i, arg)
            seqs.write("test/data/sample_arg_joint3/%d.fa" % i)
        
        lk = arghmm.calc_joint_prob(arg, seqs, mu=mu, rho=rho, times=times,
                                    popsize=n)
        print "real joint", lk

        y = []
        
        util.tic("sample ARG")
        arg2 = arghmm.sample_arg(seqs, rho=rho2, mu=mu, times=times,
                                 popsize=n)
        lk2 = arghmm.calc_joint_prob(arg2, seqs, mu=mu, rho=rho, times=times,
                                     popsize=n)
        y.append(lk2)
        util.toc()

        for i in range(5):
            util.tic("remax ARG %d" % i)
            arg2 = arghmm.remax_arg(arg2, seqs, rho=rho, mu=mu, times=times,
                                    popsize=n, refine=1)
            lk2 = arghmm.calc_joint_prob(
                arg2, seqs, mu=mu, rho=rho, times=times, popsize=n)
            y.append(lk2)
            print lk2
            util.toc()

        for i in range(40):
            util.tic("resample ARG %d" % i)
            arg2 = arghmm.resample_arg(arg2, seqs, rho=rho, mu=mu, times=times,
                                     popsize=n, refine=refine, nremove=nremove)
            util.toc()
            lk2 = arghmm.calc_joint_prob(arg2, seqs, mu=mu, rho=rho,
                                         times=times, popsize=n)
            y.append(lk2)
            print lk2

        
        p = plot(y)
        makedirs("data/sample_arg_joint3/")
        write_list("data/sample_arg_joint3/joint.txt", [lk] + y)
        p.plot([0, len(y)], [lk, lk], style="lines")
        
        
        pause()



    def test_max_arg_joint(self):
        """
        Plot the ARG joint prob from a fully sampled ARG
        """

        k = 8
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        refine = 4; nremove = 1
        write = False
        if write:
            make_clean_dir("test/data/sample_arg_joint")

        print "times", times

        names = []
        rx = []
        ry = []
        util.tic("plot")
        for i in range(40):
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times=times)
            seqs = arglib.make_alignment(arg, muts)
            if write:
                arglib.write_arg("test/data/sample_arg_joint/%d.arg" % i, arg)
                seqs.write("test/data/sample_arg_joint/%d.fa" % i)

            
            lk = arghmm.calc_joint_prob(arg, seqs, mu=mu, rho=rho, times=times)

            for j in range(1):
                util.tic("sample ARG %d, %d" % (i, j))
                arg2 = arghmm.sample_arg(seqs, rho=rho2, mu=mu, times=times,
                                         refine=0)
                arg2 = arghmm.remax_arg(arg2, seqs, rho=rho, mu=mu,
                                        times=times, refine=5)
                arg2 = arghmm.resample_arg(arg2, seqs, rho=rho2, mu=mu,
                                           times=times, refine=refine)

                util.toc()

                lk2 = arghmm.calc_joint_prob(arg2, seqs, mu=mu, rho=rho,
                                             times=times)
                rx.append(lk)
                ry.append(lk2)
                names.append([i, j])
                if write:
                    arglib.write_arg("test/data/sample_arg_joint/%d-%d.arg" %
                                     (i, j), arg2)
        util.toc()

        if write:
            data = [[i, j, a, b] for (i, j), a, b in zip(names, rx, ry)]
            write_delim("test/data/sample_arg_joint/data.txt", data)
        

        print "avg ratio:", mean([safediv(i, j, 0) for i, j in zip(ry, rx)])

        p = plot(rx, ry,
                 xlab="true ARG joint probability",
                 ylab="inferred ARG joint probability")
        p.plot([min(rx), max(rx)], [min(rx), max(rx)], style="lines")
        
        pause()

        
    def test_max_arg_joint2(self):
        """
        Plot the recombinations from a fully sampled ARG over many Gibb iters
        """
        k = 8
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        nremove = 1; refine = 1
        #nremove = 1; refine = 1
        write = True
        if write:
            make_clean_dir("test/data/max_arg_joint2")
        
        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg, mu, times=times)
        seqs = arglib.make_alignment(arg, muts)
        if write:
            arglib.write_arg("test/data/max_arg_joint2/arg.arg", arg)
            seqs.write("test/data/max_arg_joint2/seqs.fa")
        
        lk = arghmm.calc_joint_prob(arg, seqs, mu=mu, rho=rho, times=times,
                                    popsize=n)
        print "real joint", lk

        y = []
        
        util.tic("sample ARG")
        arg2 = arghmm.sample_arg(seqs, rho=rho2, mu=mu, times=times,
                                 popsize=n)
        util.toc()
        if write:
            arglib.write_arg("test/data/max_arg_joint2/sample.arg", arg2)

        
        lk2 = arghmm.calc_joint_prob(arg2, seqs, mu=mu, rho=rho, times=times,
                                     popsize=n)
        y.append(lk2)
        print lk2

        for i in range(20):
            util.tic("remax ARG %d" % i)
            arg2 = arghmm.remax_arg(arg2, seqs, rho=rho, mu=mu, times=times,
                                    popsize=n, refine=refine, nremove=nremove)
            util.toc()
            lk2 = arghmm.calc_joint_prob(arg2, seqs, mu=mu, rho=rho,
                                         times=times, popsize=n)
            y.append(lk2)
            print lk2
            if write:
                arglib.write_arg("test/data/max_arg_joint2/%d.arg" % i, arg2)

        
        p = plot(y)
        p.plot([0, len(y)], [lk, lk], style="lines")
        
        
        pause()


        
    def test_sample_arg_joint_seq(self):
        """
        Plot the recombinations from a fully sampled ARG over sequential iters
        """
        k = 5
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 50000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)

        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg, mu, times=times)
        seqs = arglib.make_alignment(arg, muts)
                
        lk = arghmm.calc_joint_prob(arg, seqs, mu=mu, rho=rho,
                                    popsize=n, times=times)
        print "real joint", lk

        y = []
        for i in range(40):
            keys = seqs.keys()
            random.shuffle(keys)
            seqs = seqs.get(keys)
            
            util.tic("sample ARG %d" % i)
            arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu, popsize=n,
                                     times=times)
            util.toc()
            lk2 = arghmm.calc_joint_prob(arg2, seqs, mu=mu, rho=rho,
                                         popsize=n, times=times)
            y.append(lk2)
            print lk2

        
        p = plot(y)
        makedirs("data/sample_arg_joint_seq/")
        write_list("data/sample_arg_joint_seq/joint.txt", [lk] + y)
        p.plot([0, len(y)], [lk, lk], style="lines")
        
        
        pause()


    def test_branch_correct_plot(self):
        
        k = 12
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 30000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)

        arg1 = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg1, mu, times=times)
        seqs = arglib.make_alignment(arg1, muts)

        # make core
        util.tic("sample core ARG")
        names = seqs.keys()[:6]
        random.shuffle(names)
        core_seqs = seqs.get(names)
        core_arg = arghmm.sample_arg(core_seqs, rho=rho, mu=mu, times=times,
                                     refine=6)
        util.toc()
                
        util.tic("sample ARG")
        arg2 = arghmm.resample_arg(core_arg, seqs, rho=rho, mu=mu, times=times)
        util.toc()

        
        coords = range(0, arg1.end, 50)
        correct = []
        
        for pos in coords:
            correct.append(get_branch_correct(arg1, arg2, pos))
        print mean(correct)

        p = plot(coords, correct, style="lines")
        pause()


    def test_max_branch_correct_plot(self):
        
        k = 12
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 30000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)

        arg1 = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg1, mu, times=times)
        seqs = arglib.make_alignment(arg1, muts)

        # make core
        util.tic("sample core ARG")
        names = seqs.keys()#[:6]
        random.shuffle(names)
        core_seqs = seqs.get(names)
        core_arg = arghmm.sample_arg(core_seqs, rho=rho, mu=mu, times=times,
                                     refine=0)
        util.toc()
                
        util.tic("sample ARG")
        arg2 = arghmm.remax_arg(core_arg, seqs, rho=rho, mu=mu, times=times,
                                refine=4)

        arg2 = arghmm.resample_arg(arg2, core_seqs, rho=rho, mu=mu,
                                   times=times, refine=4)
        util.toc()

        
        coords = range(0, arg1.end, 50)
        correct = []
        
        for pos in coords:
            correct.append(get_branch_correct(arg1, arg2, pos))
        print mean(correct)

        p = plot(coords, correct, style="lines")
        pause()



    def test_branch_correct(self):
        
        k = 20
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)

        coords = range(0, length, 50)
        x = []
        for i in range(20):
            arg1 = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                          times=times)
            muts = arghmm.sample_arg_mutations(arg1, mu, times=times)
            seqs = arglib.make_alignment(arg1, muts)
            
            # make core
            util.tic("sample core ARG %d" % i)
            core_seqs = seqs.get(seqs.keys()[:5])
            core_arg = arghmm.sample_arg(core_seqs, rho=rho, mu=mu, times=times,
                                         refine=0)
            util.toc()

            arg2 = core_arg
            util.tic("sample ARG %d" % i)
            arg2 = arghmm.resample_arg(core_arg, seqs,
                                       rho=rho, mu=mu, times=times)
            util.toc()
            
            x.append(mean([get_branch_correct(arg1, arg2, pos)
                           for pos in coords]))
            print x[-1]

        print "avg", mean(x)
        p = plot(sorted(x))
        pause()


    def test_sample_arg_seq(self):
        """
        Plot the recombinations from a fully sampled ARG over sequential iters
        """
        k = 12
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        write = True

        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg, mu, times=times)
        seqs = arglib.make_alignment(arg, muts)
        if write:
            makedirs("test/data/sample_arg_seq/")
            arglib.write_arg("test/data/sample_arg_seq/arg.arg", arg)
            seqs.write("test/data/sample_arg_seq/seqs.fa")

                
        lk = arghmm.calc_joint_prob(arg, seqs, mu=mu, rho=rho, times=times)
        print "real joint", lk

        y = []
        for i in range(200):
            keys = seqs.keys()
            random.shuffle(keys)
            seqs = seqs.get(keys)
            
            util.tic("sample ARG %d" % i)
            arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu, times=times)
            util.toc()
            lk2 = arghmm.calc_joint_prob(arg2, seqs, mu=mu, rho=rho, times=times)
            y.append(lk2)
            print lk2
            if write:
                arglib.write_arg("test/data/sample_arg_seq/%d.arg" % i, arg2)
        
        p = plot(y)
        p.plot([0, len(y)], [lk, lk], style="lines")
        
        
        pause()


    def test_sample_gibbs_ages(self):
        """plot node ages across gibbs iterations"""

        k = 12
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = 10000
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        print times

        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg, mu, times=times)
        seqs = arglib.make_alignment(arg, muts)

        arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu, times=times)
        args = [arg2]
        
        for i in range(200):
            util.tic("sample ARG %d" % i)
            arg2 = arghmm.resample_arg(arg2, seqs, rho=rho, mu=mu, times=times)
            args.append(arg2)
            util.toc()

        trees = [a.get_marginal_tree(length/2) for a in args]
        x = []
        y = []
        for i, tree in enumerate(trees):
            for node in tree:
                if len(node.children) != 1:
                    x.append(i)
                    y.append(node.age)
        p = plot(x, y, ylog=10, ymin=10)

        # plot true node ages
        tree = arg.get_marginal_tree(length/2)
        y2 = [node.age for node in tree
              if len(node.children) != 1]
        p.plot([0] * len(y2), y2)

        pause()


    def test_sample_arg_popsizes(self):
        """
        Fully sample an ARG from stratch using API
        """

        k = 2
        rho = 1.5e-8
        mu = 2.5e-8
        length = 10000000
        times = arghmm.get_time_points(ntimes=50, maxtime=160000)
        times2 = arghmm.get_time_points(ntimes=50, maxtime=160000)
        popsizes = [1e4 * (61.-i)/60. for i in range(len(times))]
        refine = 0

        util.tic("sim ARG")
        #arg = arglib.sample_arg_smc(k, 2 * popsizes[0],
        #                            rho, start=0, end=length)        
        arg = arghmm.sample_arg_dsmc(k, [2*p for p in popsizes],
                                     rho, start=0, end=length, times=times)
        util.toc()
        popsizes2 = arghmm.est_arg_popsizes(arg, times=times2)
        
        print popsizes2
        p = plot(times, popsizes, xlog=10)
        p.plot(times2[1:], popsizes2)

        pause()


    def test_sample_arg_popsizes2(self):
        """
        Fully sample an ARG from stratch using API
        """

        k = 12
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 100000
        times = arghmm.get_time_points(ntimes=50, maxtime=160000)
        popsizes = [1e4 * (61.-i)/60. for i in range(len(times))]
        refine = 0

        util.tic("sim ARG")
        #arg = arglib.sample_arg_smc(k, 2 * popsizes[0],
        #                            rho, start=0, end=length)        
        arg = arghmm.sample_arg_dsmc(k, [2*p for p in popsizes],
                                     rho, start=0, end=length, times=times)
        util.toc()

        muts = arghmm.sample_arg_mutations(arg, mu, times=times)
        seqs = arglib.make_alignment(arg, muts)

        #arg2 = arghmm.sample_arg(seqs.get(seqs.keys()[:6]),
        #                         rho=rho, mu=mu, times=times, refine=5)
        #arg2 = arghmm.resample_arg(arg2, seqs, rho=rho, mu=mu, times=times)

        popsizes2 = [0] * (len(times) - 1)
        nsamples = 1
        for i in range(nsamples):
            arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu, times=times,
                                     popsizes=popsizes,
                                     refine=0, verbose=True)
            popsizes3 = arghmm.est_arg_popsizes(arg2, times=times)
            popsizes2 = vadd(popsizes2, popsizes3)
        popsizes2 = vdivs(popsizes2, float(nsamples))

        #popsizes2 = arghmm.est_arg_popsizes(arg, times=times)
        
        print popsizes2
        p = plot(popsizes)
        p.plot(popsizes2)
        #p.plot(popsizes3)        

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


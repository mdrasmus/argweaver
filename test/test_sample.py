

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

def get_arg_tmrca(arg, steps=100):
    x = []
    for pos in range(arg.start, arg.end, (arg.end - arg.start) / steps):
        tree = arg.get_marginal_tree(pos)
        #x.append(tree.root.age)
        x.append(max(get_tree_coals(tree)))
    return x

def get_args_tmrca(args, steps=100):
    xs = [get_arg_tmrca(arg, steps=steps) for arg in args]
    return [mean(row) for row in transpose(xs)]

def get_tree_coals(tree):
    leaves = list(tree.leaf_names())
    leaves.sort()
    pairs = [(leaves[i], leaves[j]) for i in range(len(leaves))
             for j in range(i+1, len(leaves))]
    tree2 = tree.get_tree()

    ages = []
    for pair in pairs:
        l = treelib.lca(mget(tree2, pair))
        ages.append(tree[l.name].age)
    
    return ages

def get_arg_coals(arg, steps=100):
    x = []
    for pos in range(arg.start, arg.end, (arg.end - arg.start) / steps):
        tree = arg.get_marginal_tree(pos)
        x.extend(get_tree_coals(tree))
    return x

def get_args_coals(args, steps=100):
    xs = [get_arg_coals(arg, steps=steps) for arg in args]
    return [mean(row) for row in transpose(xs)]




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

        k = 30
        n = 1e4
        rho = 1.5e-8
        mu = 2.5e-8
        length = int(200e3)
        #times = arghmm.get_time_points(ntimes=80, maxtime=200000)
        times = arghmm.get_time_points(ntimes=20, maxtime=180000)
        #times = [i/30. * 200000 for i in range(31)]

        r1 = []; r2 = []
        for i in range(1, 100):
            print i
            arg = arglib.sample_arg_smc(k, 2*n, i/100. * rho,
                                         start=0, end=length)
            arg2 = arghmm.sample_arg_dsmc(k, 2*n, i/100. * rho, times=times,
                                          start=0, end=length)
            r1.append(ilen(j for j in arg if j.event == "recomb"))
            r2.append(ilen(j for j in arg2 if j.event == "recomb"))

        p = plot(r1, r2, xlab="# recomb SMC", ylab="# recomb DSMC")
        p.plot([min(r1), max(r1)], [min(r1), max(r1)], style="lines")
        
        pause()

    def test_sim_smc_cmp_recomb(self):
        """
        Simulate from coal_recomb and compare to SMC
        """

        k = 20
        #n = 1e4
        n = 17241.0 / 2
        rho = 1.5e-8
        mu = 2.5e-8
        length = int(100e3)
        #times = arghmm.get_time_points(ntimes=80, maxtime=200000)
        times = arghmm.get_time_points(ntimes=20, maxtime=180000)
        #times = [i/30. * 200000 for i in range(31)]

        r1 = []; r2 = []
        nsamples = 50
        for i in range(1, nsamples):
            print i
            rho2 = i/float(nsamples) * rho
            arg = arglib.sample_arg(k, 2*n, rho2,
                                    start=0, end=length)
            #arg = arglib.smcify_arg(arg)
            arg2 = arghmm.sample_arg_dsmc(k, 2*n, rho2, times=times,
                                          start=0, end=length)
            #arg2 = arglib.sample_arg_smc(k, 2*n, rho2, 
            #                             start=0, end=length)
            r1.append(ilen(arghmm.iter_visible_recombs(arg)))
            r2.append(ilen(j for j in arg2 if j.event == "recomb"))

        p = plot(r1, r2, xlab="# recomb coal_recomb", ylab="# recomb DSMC")
        p.plot([min(r1), max(r1)], [min(r1), max(r1)], style="lines")

        print mean(r1), sdev(r1), mean(r2), sdev(r2)
        pause()


    def test_sim_recomb_coal(self):
        """
        Simulate from coal_recomb and compare to SMC
        """

        k = 20
        n = 1e4
        mu = 2.5e-8
        rho = mu / 2.0
        length = int(100e3)

        r = []; s = []
        nsamples = 40
        for i in range(1, nsamples):
            print i, 
            rho2 = rho
            arg = arglib.sample_arg(k, 2*n, rho2,
                                    start=0, end=length)
            muts = arglib.sample_arg_mutations(arg, mu)
            r.append(ilen(arghmm.iter_visible_recombs(arg)))            
            s.append(len(muts))
            print r[-1], s[-1]
        print mean(r), mean(s)
            



    def test_sim_dsmc_cmp_recomb3(self):
        """
        Simulate from DSMC and compare to SMC
        """

        k = 20
        #n = 1e4
        n = 17241.0 / 2
        rho = 1.5e-8
        mu = 2.5e-8
        length = int(500e3)
        #times = arghmm.get_time_points(ntimes=80, maxtime=200000)
        times = arghmm.get_time_points(ntimes=20, maxtime=180000)
        #times = [i/30. * 200000 for i in range(31)]

        r1 = []; r2 = []
        nsamples = 30
        for i in range(1, nsamples):
            print i
            rho2 = rho
            arg = arglib.sample_arg(k, 2*n, rho2,
                                    start=0, end=length)
            arg = arglib.smcify_arg(arg)
            arg2 = arghmm.sample_arg_dsmc(k, 2*n, rho2, times=times,
                                          start=0, end=length)
            r1.append(ilen(j for j in arg if j.event == "recomb"))
            r2.append(ilen(j for j in arg2 if j.event == "recomb"))

        p = plot(r1, r2, xlab="# recomb coal_recomb", ylab="# recomb DSMC")
        p.plot([min(r1), max(r1)], [min(r1), max(r1)], style="lines")

        print mean(r1), sdev(r1), mean(r2), sdev(r2)
        
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


    def test_sim_dsmc_cmp_muts(self):
        """
        Simulate from DSMC and compare to SMC
        """

        k = 20
        n = 1e4
        rho = 1.5e-8
        mu = 2.5e-8
        length = int(200e3)
        #times = arghmm.get_time_points(ntimes=80, maxtime=200000)
        times = arghmm.get_time_points(ntimes=20, maxtime=180000)
        #times = [i/30. * 200000 for i in range(31)]

        r1 = []; r2 = []
        nsamples = 20
        for i in range(1, nsamples+1):
            arg = arglib.sample_arg_smc(k, 2*n, rho,
                                         start=0, end=length)
            arg2 = arghmm.sample_arg_dsmc(k, 2*n, rho, times=times,
                                          start=0, end=length)
            mut = arglib.sample_arg_mutations(arg, i/float(nsamples) * mu)
            mut2 = arghmm.sample_arg_mutations(arg2, i/float(nsamples) * mu, times=times)
            r1.append(len(mut))
            r2.append(len(mut2))
            print i, r1[-1], r2[-1]

        p = plot(r1, r2, xlab="# mut SMC", ylab="# mut DSMC")
        p.plot([min(r1), max(r1)], [min(r1), max(r1)], style="lines")
        
        pause()

    #------------------------------------------


    def test_sample_thread(self):
        """
        Sample a thread from the posterior of the ArgHmm
        """

        k = 2
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = 1000
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
            path = arghmm.sample_posterior(model, length, verbose=True)
            thread2 = [times[model.states[pos][state][1]]
                       for pos, state in enumerate(path)]
            p.plot(thread2, style="lines")
            #util.write_list("test/data/sample.thread", path)

        pause()


    def test_sample_thread2(self):
        """
        Sample many threads from the ArgHmm
        """

        k = 30
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(1000e3) / 20
        times = arghmm.get_time_points(ntimes=20)

        x = []
        y = []

        for i in range(10):
            util.tic("sim")
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times)
            seqs = arglib.make_alignment(arg, muts)
            util.toc()

            new_name = "n%d" % (k-1)
            thread = list(arghmm.iter_chrom_thread(arg, arg[new_name],
                                                   by_block=False))

            # remove chrom
            util.tic("remove thread")
            arg = arghmm.remove_arg_thread(arg, new_name)
            util.toc()

            for j in xrange(1):
                arg2 = arghmm.resample_arg(arg, seqs, popsizes=n,
                                           times=times, rho=rho, mu=mu,
                                           refine=0, verbose=True)
                thread2 = list(arghmm.iter_chrom_thread(
                    arg2, arg2[new_name], by_block=False))
                x.extend(cget(thread, 1)[::100])
                y.extend(cget(thread2, 1)[::100])

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

        k = 6
        n = 1e4
        rho = 1.5e-8 * 20 / 3
        mu = 2.5e-8 * 20
        length = int(200e3) / 20
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        refine = 0
        
        arg = arglib.sample_arg(k, 2*n, rho, start=0, end=length)
        muts = arglib.sample_arg_mutations(arg, mu)
        seqs = arglib.make_alignment(arg, muts)
        arg.write("test/data/sample_arg2.arg")
        seqs.write("test/data/sample_arg2.fa")

        print len(muts)
        
        util.tic("sample ARG")
        arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu, times=times,
                                 refine=refine, verbose=True)
        util.toc()

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
                                 refine=refine, carg=True, verbose=True)
        util.toc()
        #print ilen(x for x in arg2 if x.event == "recomb")

        util.tic("sample ARG region")
        arg2 = arghmm.resample_arg_region(arg2, seqs, region[0], region[1],
                                          rho=rho, mu=mu, times=times,
                                          refine=10, verbose=True)
        util.toc()
        #print ilen(x for x in arg2 if x.event == "recomb")
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


    def test_sample_thread_recomb(self):
        """
        Plot the ARG joint prob from a fully sampled ARG
        """

        k = 6
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(1000e3) / 20
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

            nrecombs = ilen(arghmm.iter_visible_recombs(arg))

            
            util.tic("sample ARG %d" % i)
            arg = arghmm.remove_arg_thread(arg, "n%d" % (k-1))
            common = ilen(arghmm.iter_visible_recombs(arg))
            
            arg2 = arghmm.resample_arg(arg, seqs, rho=rho, mu=mu, times=times,
                                       refine=refine, nremove=nremove)
            util.toc()

            nrecombs2 = ilen(arghmm.iter_visible_recombs(arg2))
            rx.append(nrecombs - common)
            ry.append(nrecombs2 - common)
            if write:
                arglib.write_arg("test/data/sample_arg_joint/%d.arg" % i, arg2)
        util.toc()        

        print "avg ratio:", mean([safediv(i, j, 0) for i, j in zip(ry, rx)])

        p = plot(rx, ry,
                 xlab="true ARG # recombinations",
                 ylab="inferred ARG # recombinations")
        p.plot([min(rx), max(rx)], [min(rx), max(rx)], style="lines")
        
        pause()
        

    def test_sample_arg_recomb(self):
        """
        Plot the recombinations from a fully sampled ARG
        """

        k = 12
        n = 1e4
        #rho = 1.5e-8 * 20        
        mu = 2.5e-8 * 20
        rho = mu / 20.
        length = int(400e3 / 20)
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        nremove = 0
        refine = 0

        print "times", times

        rx = []
        ry = []
        util.tic("plot")
        for i in range(20):
            #arg = arglib.sample_arg_smc(k, 2*n, rho, start=0, end=length)
            #arg = arglib.smcify_arg(arg)
            
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times=times)
            seqs = arglib.make_alignment(arg, muts)
            
            nrecombs = ilen(arghmm.iter_visible_recombs(arg))

            for j in range(1):
                util.tic("sample ARG %d, %d" % (i, j))
                arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu, times=times,
                                         carg=True)
                arg2 = arghmm.resample_climb_arg(arg2, seqs,
                                                 rho=rho, mu=mu, times=times,
                                                 refine=100, carg=True)
                #arg2 = arghmm.resample_all_arg(arg, seqs, rho=rho, mu=mu,
                #                               times=times, refine=100)
                arg2 = arghmm.resample_mcmc_arg(arg2, seqs, rho=rho, mu=mu,
                    times=times, refine=refine, carg=True)
                util.toc()

                nrecombs2 = arghmm.get_local_trees_ntrees(arg2[0]) - 1
                rx.append(nrecombs)
                ry.append(nrecombs2)
                print i, rx[-1], ry[-1]
        util.toc()

        print "avg ratio:", mean([safediv(i, j, 0) for i, j in zip(ry, rx)])

        p = plot(rx, ry,
                 xlab="true # recombinations",
                 ylab="inferred # recombinations")
        p.plot([min(rx), max(rx)], [min(rx), max(rx)], style="lines")

        
        
        pause()



    def test_sample_arg_recomb_seq_gibbs(self):
        """
        Plot the recombinations from a fully sampled ARG
        """

        k = 8
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(200e3) / 20
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        seqiters=1; gibbsiters=12

        print "times", times

        rx = []
        ry = []
        util.tic("plot")
        for i in range(20):
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times=times)
            seqs = arglib.make_alignment(arg, muts)
            
            nrecombs = ilen(arghmm.iter_visible_recombs(arg))

            for j in range(1):
                util.tic("sample ARG %d, %d" % (i, j))
                arg2 = arghmm.sample_arg_seq_gibbs(
                    seqs, rho=rho, mu=mu, times=times,
                    seqiters=seqiters, gibbsiters=gibbsiters)
                util.toc()
                
                nrecombs2 = ilen(arghmm.iter_visible_recombs(arg2))
                rx.append(nrecombs)
                ry.append(nrecombs2)
                print nrecombs, nrecombs2, nrecombs2/float(nrecombs)
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

        k = 20
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho 
        mu = 2.5e-8 * 20
        length = int(200e3) / 20
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        write = False
        #nremove = 1; refine = 1
        nremove = 1; refine = 20

        makedirs("test/data/sample_arg_recomb2/")

        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        #arg = arglib.sample_arg(k, 2*n, rho, start=0, end=length)
        muts = arghmm.sample_arg_mutations(arg, mu, times=times)
        seqs = arghmm.make_alignment(arg, muts)
        if write:
            arglib.write_arg("test/data/sample_arg_recomb2/arg.arg", arg)
            seqs.write("test/data/sample_arg_recomb2/seqs.fa")
            
        nrecombs = ilen(arghmm.iter_visible_recombs(arg))
        print "real # recombs", nrecombs

        y = []
        
        util.tic("sample ARG")
        arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu, times=times,
                                 refine=0, carg=True)
        #arg2 = arghmm.sample_arg(seqs.get(seqs.keys()[:6]),
        #                         rho=rho, mu=mu, times=times,
        #                         refine=10, carg=True)
        #arg2 = arghmm.resample_arg(arg2, seqs, rho=rho, mu=mu,
        #                           times=times, refine=10, carg=True)
        util.toc()
        #arg2 = arg

        nrecombs2 = arghmm.get_local_trees_ntrees(arg2[0]) - 1
        #nrecombs2 = ilen(arghmm.iter_visible_recombs(arg2))
        y.append(nrecombs2)
        print nrecombs2

        for i in range(1000):
            util.tic("resample ARG %d" % i)
            arg2 = arghmm.resample_mcmc_arg(arg2, seqs, rho=rho, mu=mu,
                    times=times, refine=1, carg=True)
            #arg2 = arghmm.resample_arg(arg2, seqs, rho=rho, mu=mu, times=times,
            #                           refine=refine, nremove=nremove)
            #arg2 = arghmm.resample_all_arg(arg2, seqs, rho=rho, mu=mu,
            #                               times=times, refine=1, carg=True)
            #arg2 = arghmm.resample_climb_arg(arg2, seqs, rho=rho, mu=mu,
            #                                 times=times, refine=1,
            #                                 carg=True)
            util.toc()
            #nrecombs2 = ilen(arghmm.iter_visible_recombs(arg2))
            nrecombs2 = arghmm.get_local_trees_ntrees(arg2[0]) - 1
            y.append(nrecombs2)
            print nrecombs2, nrecombs

            if write:
                arglib.write_arg("test/data/sample_arg_recomb2/%d.arg" % i,
                                 arg2)

        
        p = plot(y, style='lines')
        p.plot([0, len(y)], [nrecombs, nrecombs], style="lines")
        
        pause()


    def test_max_arg_recomb2(self):
        """
        Plot the recombinations from a fully sampled ARG over many Gibb iters
        """

        k = 8
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = int(20e3) / 20
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

        k = 6
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(1000e3) / 20
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        refine = 0

        print "times", times

        rx = []
        ry = []
        util.tic("plot")
        for i in range(20):
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times=times)
            seqs = arghmm.make_alignment(arg, muts)
            
            arglen = arglib.arglen(arg)

            for j in range(1):
                util.tic("sample ARG %d, %d" % (i, j))
                arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu, times=times,
                                         carg=True)
                arg2 = arghmm.resample_climb_arg(arg2, seqs,
                                                 rho=rho, mu=mu, times=times,
                                                 refine=5, carg=True)
                arg2 = arghmm.resample_all_arg(arg2, seqs,
                                               rho=rho, mu=mu, times=times,
                                               refine=30)
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

        k = 12
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(200e3) / 20
        times = arghmm.get_time_points(ntimes=20, maxtime=160000)

        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        arg.set_ancestral()
        muts = arghmm.sample_arg_mutations(arg, mu, times)
        seqs = arghmm.make_alignment(arg, muts)
            
        arglen = arglib.arglen(arg)
        print "real # arglen %e" % arglen

        y = []
        
        util.tic("sample ARG")
        arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu, times=times)
        util.toc()
        
        arglen2 = arglib.arglen(arg2)
        y.append(arglen2)

        for i in range(50):
            util.tic("climb ARG %d" % i)
            arg2 = arghmm.resample_climb_arg(arg2, seqs, rho=rho, mu=mu,
                                             times=times)
            util.toc()
            arglen2 = arglib.arglen(arg2)
            y.append(arglen2)
            print "%e %e" % (arglen2, arglen)


        for i in range(300):
            util.tic("resample ARG %d" % i)
            arg2 = arghmm.resample_all_arg(arg2, seqs, rho=rho, mu=mu,
                                           times=times)
            util.toc()
            arglen2 = arglib.arglen(arg2)
            y.append(arglen2)
            print "%e %e" % (arglen2, arglen)

        
        p = plot(y)
        makedirs("test/data/sample_arg_arglen2/")
        write_list("test/data/sample_arg_arglen2/arglen.txt", [arglen] + y)
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
        makedirs("test/data/max_arg_arglen2/")
        write_list("test/data/max_arg_arglen2/arglen.txt", [arglen] + y)
        p.plot([0, len(y)], [arglen, arglen], style="lines")
        
        pause()

        

    def test_sample_arg_lk(self):
        """
        Plot the ARG likelihood from a fully sampled ARG
        """

        k = 12
        n = 1e4
        rho = 1.5e-8 * 20
        rho2 = rho
        mu = 2.5e-8 * 20
        length = int(200e3) / 20
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        refine = 20
        write = False
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

            for j in range(1):
                util.tic("sample ARG %d, %d" % (i, j))
                arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu, times=times,
                                         refine=1, carg=True)

                #arg2 = arghmm.sample_all_arg(seqs, rho=rho2, mu=mu, times=times,
                #                             refine=refine)

                
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
        length = int(400e3) / 20
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
        #arg2 = arg
        util.toc()
        
        lk2 = arghmm.calc_likelihood(arg2, seqs, mu=mu, times=times)
        y.append(lk2)


        for i in range(50):
            util.tic("resample ARG %d" % i)
            arg2 = arghmm.resample_arg(arg2, seqs, rho=rho, mu=mu, times=times,
                                       refine=1, carg=True)
            util.toc()
            lk2 = arghmm.calc_likelihood(arg2, seqs, mu=mu, times=times,
                                         delete_arg=False)
            y.append(lk2)
            print lk2

        
        p = plot(y)
        makedirs("test/data/sample_arg_lk2/")
        write_list("test/data/sample_arg_lk2/lk.txt", [lk] + y)
        p.plot([0, len(y)], [lk, lk], style="lines")
        
        
        pause()


    
    def test_sample_arg_prior(self):
        """
        Plot the ARG joint prob from a fully sampled ARG
        """

        k = 12
        n = 1e4
        rho = 1.5e-8 * 20 / 10
        mu = 2.5e-8 * 20
        length = int(100e3) / 20
        times = arghmm.get_time_points(ntimes=20, maxtime=180000)
        refine = 6*10; nremove = 1
        #refine = 0;
        write = False
        if write:
            make_clean_dir("test/data/sample_arg_joint")

        names = []
        rx = []
        ry = []
        util.tic("plot")
        for i in range(15):
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times=times)
            seqs = arglib.make_alignment(arg, muts)
            if write:
                arglib.write_arg("test/data/sample_arg_joint/%d.arg" % i, arg)
                seqs.write("test/data/sample_arg_joint/%d.fa" % i)

            
            lk = arghmm.calc_prior_prob(arg, seqs, rho=rho, times=times)

            for j in range(1):
                util.tic("sample ARG %d, %d" % (i, j))
                arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu, times=times,
                                         refine=10, carg=True)
                util.toc()

                lk2 = arghmm.calc_prior_prob(arg2, seqs, rho=rho, times=times)
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
                 xlab="true ARG prior probability",
                 ylab="inferred ARG prior probability")
        p.plot([min(rx), max(rx)], [min(rx), max(rx)], style="lines")
        
        pause()



    def test_sample_thread_joint(self):
        """
        Plot the ARG joint prob from a fully sampled ARG
        """

        k = 6
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(100e3) / 20
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

            
            util.tic("sample ARG %d" % i)
            arg = arghmm.remove_arg_thread(arg, "n%d" % (k-1))
            arg2 = arghmm.resample_arg(arg, seqs, rho=rho, mu=mu, times=times,
                                       refine=refine, nremove=nremove,
                                       carg=True)
            util.toc()

            lk2 = arghmm.calc_joint_prob(arg2, seqs, mu=mu, rho=rho,
                                         times=times)
            rx.append(lk)
            ry.append(lk2)
            if write:
                arglib.write_arg("test/data/sample_arg_joint/%d.arg" % i, arg2)
        util.toc()        

        print "avg ratio:", mean([safediv(i, j, 0) for i, j in zip(ry, rx)])

        p = plot(rx, ry,
                 xlab="true ARG joint probability",
                 ylab="inferred ARG joint probability")
        p.plot([min(rx), max(rx)], [min(rx), max(rx)], style="lines")
        
        pause()



    def test_sample_arg_joint(self):
        """
        Plot the ARG joint prob from a fully sampled ARG
        """

        k = 12
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        #mu = 1.5e-8 * 20
        length = int(200e3) / 20
        times = arghmm.get_time_points(ntimes=20, maxtime=180000)
        climb = 50; refine = 0
        write = False
        if write:
            make_clean_dir("test/data/sample_arg_joint")

        arghmm.setLogLevel(1)

        names = []
        rx = []
        ry = []
        util.tic("plot")
        for i in range(20):
            #arg = arglib.sample_arg(k, 2*n, rho, start=0, end=length)
            #arg = arglib.smcify_arg(arg)
            arg = arglib.sample_arg_smc(k, 2*n, rho, start=0, end=length)
            arghmm.discretize_arg(arg, times, round_age="closer")
            #arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
            #                             times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times=times)
            seqs = arghmm.make_alignment(arg, muts)
            if write:
                arglib.write_arg("test/data/sample_arg_joint/%d.arg" % i, arg)
                seqs.write("test/data/sample_arg_joint/%d.fa" % i)

            
            lk = arghmm.calc_joint_prob(arg, seqs, mu=mu, rho=rho, times=times)

            for j in range(1):
                util.tic("sample ARG %d, %d" % (i, j))
                
                arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu, times=times,
                                         refine=0, carg=True)
                arg2 = arghmm.resample_climb_arg(arg2,
                    seqs, rho=rho, mu=mu, times=times, popsizes=n,
                    refine=climb, recomb_pref=.9, carg=True)
                arg2 = arghmm.resample_mcmc_arg(arg2,
                    seqs, rho=rho, mu=mu, times=times, popsizes=n,
                    refine=refine, carg=True)
                #arg2 = arghmm.resample_all_arg(arg2,
                #    seqs, rho=rho, mu=mu, times=times, popsizes=n,
                #    refine=refine, carg=True)
                util.toc()

                lk2 = arghmm.calc_joint_prob(arg2, seqs, mu=mu, rho=rho,
                                             times=times)
                print lk, lk2
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
        length = int(400e3) / 20
        times = arghmm.get_time_points(ntimes=20, maxtime=180000)
        refine = 5; nremove = 1; core=6
        write = False
        if write:
            make_clean_dir("test/data/sample_arg_joint")


        names = []
        rx = []
        ry = []
        util.tic("plot")
        for i in range(30):
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times=times)
            seqs = arghmm.make_alignment(arg, muts)
            if write:
                arglib.write_arg("test/data/sample_arg_joint/%d.arg" % i, arg)
                seqs.write("test/data/sample_arg_joint/%d.fa" % i)

            
            lk = arghmm.calc_joint_prob(arg, seqs, mu=mu, rho=rho, times=times)

            for j in range(1):
                util.tic("sample ARG %d, %d" % (i, j))

                # infer core
                names = seqs.keys()
                random.shuffle(names)
                arg2 = arghmm.sample_arg(seqs.get(names[:core]),
                                         rho=rho, mu=mu, times=times,
                                         refine=20, carg=True)
                arg2 = arghmm.resample_arg(arg2, seqs, rho=rho, mu=mu,
                                           times=times, refine=3, carg=True)
                arg2 = arghmm.resample_climb_arg(
                    arg2, seqs, rho=rho, mu=mu, times=times,
                    popsizes=n, refine=50, nclimb=10, carg=True)
                arg2 = arghmm.resample_all_arg(arg2, seqs, rho=rho, mu=mu,
                                               times=times, refine=100,
                                               carg=True)
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
        length = int(200e3) / 20
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        refine = 10; nremove = 1;

        names = []
        rx = []
        ry = []
        util.tic("plot")
        for i in range(15):
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
                    arg2, seqs, niters=0, width=500,
                    rho=rho, mu=mu,popsize=n, times=times, carg=True,
                    verbose=True)
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
        k = 12
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(200e3) / 20
        times = arghmm.get_time_points(ntimes=20, maxtime=200e3)
        climb = 200; refine = 1000
        niters2 = 5
        window = int(200e3) / 20
        nremove = 1; 
        write = False

        
        
        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg, mu, times=times)
        seqs = arghmm.make_alignment(arg, muts)
        if write:
            arglib.write_arg("test/data/sample_arg_joint/%d.arg" % i, arg)
            seqs.write("test/data/sample_arg_joint/%d.fa" % i)
        
        lk = arghmm.calc_joint_prob(arg, seqs, mu=mu, rho=rho, times=times,
                                    popsizes=n)
        nrecombs = ilen(x for x in arg if x.event == "recomb")
        print "real joint", lk, nrecombs

        y = []
        
        util.tic("sample ARG")
        arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu, times=times,
                                 popsizes=n, refine=0, carg=True)
        util.toc()
        
        lk2 = arghmm.calc_joint_prob(arg2, seqs, mu=mu, rho=rho, times=times,
                                     popsizes=n, delete_arg=False)
        y.append(lk2)

        
        for i in range(climb):
            util.tic("climb ARG %d" % i)
            arg2 = arghmm.resample_climb_arg(
                arg2, seqs, rho=rho, mu=mu, times=times,
                popsizes=n, recomb_pref=.9, carg=True)
            util.toc()
            nrecombs2 = arghmm.get_local_trees_ntrees(arg2[0]) - 1
            lk2 = arghmm.calc_joint_prob(arg2, seqs, mu=mu, rho=rho,
                                         times=times, popsizes=n,
                                         delete_arg=False)
            y.append(lk2)
            print lk2, lk
            print nrecombs2, nrecombs

        arghmm.setLogLevel(0)
        
        for i in range(refine):
            util.tic("resample ARG %d" % i)
            #arg2 = arghmm.resample_arg(
            #    arg2, seqs, rho=rho, mu=mu, times=times,
            #    popsizes=n, refine=1, carg=True)

            #arg2 = arghmm.resample_all_arg(
            #    arg2, seqs, rho=rho, mu=mu, times=times,
            #    popsizes=n, refine=1, carg=True)
            
            arg2 = arghmm.resample_mcmc_arg(
                arg2, seqs, rho=rho, mu=mu, times=times,
                popsizes=n, refine=1, carg=True,
                window=window, niters2=niters2)
            util.toc()
            nrecombs2 = arghmm.get_local_trees_ntrees(arg2[0]) - 1
            lk2 = arghmm.calc_joint_prob(arg2, seqs, mu=mu, rho=rho,
                                         times=times, popsizes=n,
                                         delete_arg=False)
            y.append(lk2)
            print lk2, lk
            print nrecombs2, nrecombs
        
        p = plot(y, style='lines')
        makedirs("test/data/sample_arg_joint2/")
        write_list("test/data/sample_arg_joint2/joint.txt", [lk] + y)
        p.plot([0, len(y)], [lk, lk], style="lines")


        #arg2 = arghmm.resample_all_arg(
        #        arg2, seqs, rho=rho, mu=mu, times=times,
        #        popsizes=n)

        # plot block starts
        #x = [start for (start, end), tree in arglib.iter_tree_tracks(arg)]
        #x2 = [start for (start, end), tree in arglib.iter_tree_tracks(arg2)]
        #q = plot(x)
        #q.plot(x2)

        pause()



    def test_sample_arg_joint_copy(self):
        """
        Plot the recombinations from a fully sampled ARG over many Gibb iters
        """
        k = 30
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(200e3) / 20
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        nremove = 1; refine = 12 * 20
        write = False
        
        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg, mu, times=times)
        seqs = arghmm.make_alignment(arg, muts)
        if write:
            arglib.write_arg("test/data/sample_arg_joint/%d.arg" % i, arg)
            seqs.write("test/data/sample_arg_joint/%d.fa" % i)
        
        lk = arghmm.calc_joint_prob(arg, seqs, mu=mu, rho=rho, times=times,
                                    popsizes=n)
        nrecombs = ilen(x for x in arg if x.event == "recomb")
        print "real joint", lk, nrecombs

        y = []
        
        util.tic("sample ARG")
        core_seqs = seqs.get(seqs.keys()[:20])
        
        arg2 = arghmm.sample_arg(core_seqs, rho=rho, mu=mu, times=times,
                                 popsizes=n, refine=0, carg=True)
        util.toc()
        
        for i in range(50):
            util.tic("resample ARG %d" % i)
            arg2 = arghmm.resample_climb_arg(
                arg2, core_seqs, rho=rho, mu=mu, times=times,
                popsizes=n, refine=1, nclimb=10, carg=True)
            util.toc()
            lk2 = arghmm.calc_joint_prob(arg2, seqs, mu=mu, rho=rho,
                                         times=times, popsizes=n,
                                         delete_arg=False)
            y.append(lk2)
            print lk2, lk

        for i in range(500):
            util.tic("resample ARG %d" % i)
            
            arg2 = arghmm.resample_all_arg(
                arg2, core_seqs, rho=rho, mu=mu, times=times,
                popsizes=n, refine=1, carg=True)
            
            trees = arghmm.arghmm_copy_trees(arg2[0])
            arg3 = (trees, arg2[1])
            arg3 = arghmm.resample_arg(arg3, seqs, rho=rho, mu=mu, times=times,
                                       popsizes=n, refine=0, carg=True)
            
            
            util.toc()
            nrecombs2 = arghmm.get_local_trees_ntrees(arg3[0]) - 1

            lk2 = arghmm.calc_joint_prob(arg3, seqs, mu=mu, rho=rho,
                                         times=times, popsizes=n,
                                         delete_arg=True)
            y.append(lk2)
            print lk2, lk
            print nrecombs2, nrecombs

            


        
        p = plot(y)
        makedirs("test/data/sample_arg_joint2/")
        write_list("test/data/sample_arg_joint2/joint.txt", [lk] + y)
        p.plot([0, len(y)], [lk, lk], style="lines")
        
        
        pause()



    def test_sample_arg_joint3(self):
        """
        Plot the recombinations from a fully sampled ARG over many Gibb iters
        """
        k = 6
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(100e3) / 20
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        nremove = 1; refine = 12 * 20
        write = False
        
        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg, mu, times=times)
        seqs = arghmm.make_alignment(arg, muts)
        if write:
            arglib.write_arg("test/data/sample_arg_joint/%d.arg" % i, arg)
            seqs.write("test/data/sample_arg_joint/%d.fa" % i)
        
        lk = arghmm.calc_joint_prob(arg, seqs, mu=mu, rho=rho, times=times,
                                    popsizes=n)
        print "real joint", lk


        p = plot([0, 200], [lk, lk], style="lines")

        for j in range(10):
            util.tic("sample ARG")
            y = []

            arg2 = arghmm.sample_arg(seqs.get(seqs.keys()[:6]),
                                     rho=rho, mu=mu, times=times,
                                     popsizes=n, refine=10, carg=True)
            arg2 = arghmm.sample_arg_seq_gibbs(
                        seqs, rho=rho, mu=mu, times=times,
                        seqiters=2, gibbsiters=1, carg=True)
            util.toc()

            lk2 = arghmm.calc_joint_prob(arg2, seqs, mu=mu, rho=rho,
                                         times=times,
                                         popsizes=n, delete_arg=False)
            y.append(lk2)

            for i in range(200):
                util.tic("resample ARG %d" % i)
                arg2 = arghmm.resample_all_arg(
                    arg2, seqs, rho=rho, mu=mu, times=times,
                    popsizes=n, refine=1, carg=True)
                util.toc()
                lk2 = arghmm.calc_joint_prob(arg2, seqs, mu=mu, rho=rho,
                                             times=times, popsizes=n,
                                             delete_arg=False)
                y.append(lk2)
                print lk2, lk

            p .plot(y)
        
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
        for i in range(10):
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
        makedirs("test/data/sample_arg_joint_seq/")
        write_list("test/data/sample_arg_joint_seq/joint.txt", [lk] + y)
        p.plot([0, len(y)], [lk, lk], style="lines")
        
        
        pause()


    def test_sample_arg_tmrca(self):
        """
        Plot the ARG joint prob from a fully sampled ARG
        """

        k = 6
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(200e3) / 20
        times = arghmm.get_time_points(ntimes=30, maxtime=180000)
        refine = 6*10; nremove = 1
        steps = 200

        names = []
        rx = []
        ry = []
        util.tic("plot")
        for i in range(1):
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times=times)
            seqs = arghmm.make_alignment(arg, muts)

            rx.extend(get_arg_tmrca(arg, steps))

            for j in range(1):
                util.tic("sample ARG %d, %d" % (i, j))
                args = []
                for l in range(30):
                    arg2 = arghmm.sample_arg(
                        seqs, rho=rho, mu=mu, times=times, refine=10)
                    arg2 = arghmm.resample_all_arg(arg2,
                        seqs, rho=rho, mu=mu, times=times, refine=20)
                    args.append(arg2)
                util.toc()

                ry.extend(get_args_tmrca(args, steps))
        util.toc()

        print "avg ratio:", mean([safediv(i, j, 0) for i, j in zip(ry, rx)])
        rx = map(log10, rx)
        ry = map(log10, ry)
        p = plot(dither(rx, .03), dither(ry, .03),
                 xlab="true ARG TMRCA",
                 ylab="inferred ARG TMRCA")
        p.plot([min(rx), max(rx)], [min(rx), max(rx)], style="lines")
        
        pause()


    def test_sample_arg_coal(self):
        """
        Plot the ARG joint prob from a fully sampled ARG
        """

        k = 12
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(200e3) / 20
        times = arghmm.get_time_points(ntimes=30, maxtime=180000)
        refine = 6*10; nremove = 1
        steps = 40

        names = []
        rx = []
        ry = []
        util.tic("plot")
        for i in range(1):
            arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                         times=times)
            muts = arghmm.sample_arg_mutations(arg, mu, times=times)
            seqs = arghmm.make_alignment(arg, muts)

            rx.extend(get_arg_coals(arg, steps))

            for j in range(1):
                util.tic("sample ARG %d, %d" % (i, j))
                args = []
                for l in range(30):
                    arg2 = arghmm.sample_arg(
                        seqs, rho=rho, mu=mu, times=times, carg=True)
                    arg2 = arghmm.resample_climb_arg(arg2, 
                        seqs, rho=rho, mu=mu, times=times, refine=200)
                    args.append(arg2)
                util.toc()

                ry.extend(get_args_coals(args, steps))
        util.toc()

        print "avg ratio:", mean([safediv(i, j, 0) for i, j in zip(ry, rx)])
        rx = map(log10, rx)
        ry = map(log10, ry)
        p = plot(dither(rx, .03), dither(ry, .03),
                 xlab="true ARG coal ages",
                 ylab="inferred ARG coal ages")
        p.plot([min(rx), max(rx)], [min(rx), max(rx)], style="lines")
        
        pause()


    def test_branch_correct_plot(self):
        
        k = 12
        n = 1e4
        rho = 1.5e-8 * 20  / 4
        mu = 2.5e-8 * 20
        length = int(200e3) / 20
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)

        arg1 = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                      times=times)
        muts = arghmm.sample_arg_mutations(arg1, mu, times=times)
        seqs = arghmm.make_alignment(arg1, muts)

        # make core
        util.tic("sample ARG")
        arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu, times=times)
        arg2 = arghmm.resample_climb_arg(arg2, seqs, rho=rho, mu=mu,
                                         times=times, refine=200)
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


    def test_branch_correct2(self):
        
        k = 12
        n = 1e4
        rho = 1.5e-8 * 20 / 4
        mu = 2.5e-8 * 20
        length = int(200e3) / 20
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)

        coords = range(0, length, 50)
        x = []

        arg1 = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                      times=times)
        muts = arghmm.sample_arg_mutations(arg1, mu, times=times)
        seqs = arghmm.make_alignment(arg1, muts)

        arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu, times=times)

        for i in range(40):
            util.tic("sample ARG %d" % i)
            arg2 = arghmm.resample_climb_arg(arg2, seqs,
                                             rho=rho, mu=mu, times=times,
                                             refine=5)
            
            x.append(mean([get_branch_correct(arg1, arg2, pos)
                           for pos in coords]))
            print x[-1]
            util.toc()
            
        p = plot(x)
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

        k = 20
        rho = 1.5e-8
        mu = 2.5e-8
        length = int(100e3)
        times = arghmm.get_time_points(ntimes=20, maxtime=120000)
        a = 30.
        b = 15
        #popsizes = [1e4 * (a - b + abs(i-b))/a for i in range(len(times))]
        popsizes = [1e4 * (a - i)/a for i in range(len(times))]
        #popsizes = [1e4 for i in range(len(times))]
        refine = 0

        util.tic("sim ARG")
        #arg = arglib.sample_arg_smc(k, 2 * popsizes[0],
        #                            rho, start=0, end=length)
        arg = arghmm.sample_arg_dsmc(k, [2*p for p in popsizes],
                                     rho, start=0, end=length, times=times)
        util.toc()

        print "recombs", ilen(x for x in arg if x.event == "recomb")

        util.tic("estimate popsizes")
        popsizes2 = arghmm.est_arg_popsizes(arg, times=times,
                                            popsize_mu=1e4, popsize_sigma=1e4)
        util.toc()
        
        print popsizes2
        p = plot(times, popsizes, xlog=10, xmin=10, ymin=0, ymax=20000)
        #p = plot(times, popsizes, xmin=10, ymin=0, ymax=20000)
        p.plot(times[1:], popsizes2)
        
        pause()


    def test_sample_arg_popsizes_window(self):
        """
        Estimate population sizes per time and window
        """

        k = 50
        rho = 1.5e-8
        mu = 2.5e-8
        length = int(500e3)
        window = int(300e3)
        times = arghmm.get_time_points(ntimes=20, maxtime=120000)
        a = 30.
        b = 15
        #popsizes = [1e4 * (a - b + abs(i-b))/a for i in range(len(times))]
        popsizes = [1e4 * (a - i)/a for i in range(len(times))]
        #popsizes = [1e4 for i in range(len(times))]
        refine = 0

        util.tic("sim ARG")
        arg = arghmm.sample_arg_dsmc(k, [2*p for p in popsizes],
                                     rho, start=0, end=length, times=times)
        util.toc()

        print "recombs", ilen(x for x in arg if x.event == "recomb")

        util.tic("estimate popsizes")
        step = 50e3
        starts = range(arg.start, arg.end - window, step)
        popsizes2 = []
        for start in starts:
            print start
            arg2 = arglib.smcify_arg(arg, start, start+window)
            popsizes2.append(arghmm.est_arg_popsizes(
                arg2, times=times,
                popsize_mu=1e4, popsize_sigma=1e4))
        util.toc()

        popsizes2 = [popsizes[:len(popsizes2[0])]] + popsizes2
        #popsizes2 = transpose(popsizes2)
        
        pc(popsizes2)
        #heatmap(popsizes2, showVals=True)
        rp.assign("X", flatten(popsizes2))
        rp.assign("l", len(popsizes2[0]))
        rp("heatmap(matrix(X, l), Rowv=NA, Colv=NA, scale='none', col=rainbow(40, start=0, end=.3))")
        #rp.heatmap(rp.matrix(popsizes2, len(popsizes2)), reorderfunc=None)
        #rp.heatmap(rp.matrix(popsizes2, len(popsizes2)), reorderfunc=None)
        
        pause()



    def test_sample_arg_popsizes_trees(self):
        """
        Fully sample an ARG from stratch using API
        """

        k = 20
        rho = 1.5e-8
        mu = 2.5e-8
        length = int(100e3)
        times = arghmm.get_time_points(ntimes=20, maxtime=160000)
        popsizes = [1e4 * (61.-i)/60. for i in range(len(times))]
        #popsizes = [1e4 for i in range(len(times))]
        refine = 0

        util.tic("sim ARG")
        #arg = arglib.sample_arg_smc(k, 2 * popsizes[0],
        #                            rho, start=0, end=length)
        arg = arghmm.sample_arg_dsmc(k, [2*p for p in popsizes],
                                     rho, start=0, end=length, times=times)
        util.toc()

        util.tic("estimate popsizes")
        popsizes2 = arghmm.est_popsizes_trees(arg, times=times,
                                              step=length/1000, verbose=True)
        util.toc()
        
        print popsizes2
        p = plot(times, popsizes, xlog=10, xmin=10, ymin=0, ymax=20000)
        p.plot(times[1:], popsizes2)

        pause()


    def test_sample_arg_popsizes_infer(self):
        """
        Fully sample an ARG from stratch using API
        """

        k = 2
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(100e6) / 20
        times = arghmm.get_time_points(ntimes=30, maxtime=160000)
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
        
        popsizes2 = [0] * (len(times) - 1)
        nsamples = 1
        for i in range(nsamples):
            arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu, times=times,
                                     popsizes=popsizes,
                                     refine=0, verbose=True)
            popsizes3 = arghmm.est_arg_popsizes(arg2, times=times)
            #popsizes3 = arghmm.est_popsizes_trees(arg2, times, length/200)
            print popsizes3
            popsizes2 = vadd(popsizes2, popsizes3)
        popsizes2 = vdivs(popsizes2, float(nsamples))

        print popsizes2
        p = plot(times, popsizes, xlog=10, xmin=10)
        p.plot(times[1:], popsizes2)

        pause()


    def test_sample_arg_popsizes_trees_infer(self):
        """
        Fully sample an ARG from stratch using API
        """

        k = 6
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(10e3) / 20
        times = arghmm.get_time_points(ntimes=20, maxtime=160000)
        popsizes = [1e4 * (61.-i)/60. for i in range(len(times))]
        refine = 5

        util.tic("sim ARG")
        arg = arghmm.sample_arg_dsmc(k, [2*p for p in popsizes],
                                     rho, start=0, end=length, times=times)
        muts = arghmm.sample_arg_mutations(arg, mu, times=times)
        seqs = arglib.make_alignment(arg, muts)
        util.toc()        
        
        popsizes2 = [0] * (len(times) - 1)
        nsamples = 1
        for i in range(nsamples):
            arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu, times=times,
                                     popsizes=popsizes,
                                     refine=refine, verbose=True, carg=True)
            popsizes3 = arghmm.est_popsizes_trees(arg2, times, length/1000,
                                                  verbose=True)
            print popsizes3
            popsizes2 = vadd(popsizes2, popsizes3)
        popsizes2 = vdivs(popsizes2, float(nsamples))

        print popsizes2
        p = plot(times, popsizes, xlog=10, xmin=10)
        p.plot(times[1:], popsizes2)

        pause()


    #------------------------------------------------------------

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


    def test_sample_removal_path(self):
        """
        Test the sampling of a branch removal path
        """

        k = 12
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(500e3) / 20
        times = arghmm.get_time_points(ntimes=20)

        util.tic("sim")
        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        trees, names = arghmm.arg2ctrees(arg,times)
        nnodes = 2*k - 1
        ntrees = arghmm.get_local_trees_ntrees(trees)
        print ntrees, nnodes
        util.toc()

        x = []
        y = []
        for i in range(20):
            path = [0] * ntrees
            node = random.randint(0, nnodes - 1)
            arghmm.arghmm_sample_arg_removal_path(trees, node, path)

            for j in xrange(ntrees):
                x.append(j)
                y.append(path[j])
        
        arghmm.delete_local_trees(trees)

        p = plot(dither(x, .4), dither(y, .4))

        pause()


    def test_sample_removal_path2(self):
        """
        Test the sampling of a branch removal path
        """

        k = 12
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(500e3) / 20
        times = arghmm.get_time_points(ntimes=20)

        util.tic("sim")
        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        trees, names = arghmm.arg2ctrees(arg,times)
        nnodes = 2*k - 1
        ntrees = arghmm.get_local_trees_ntrees(trees)
        print ntrees, nnodes
        util.toc()

        x = []
        y = []
        for i in range(3):
            path = [0] * ntrees
            node = random.randint(0, nnodes - 1)
            #pos = random.randint(0, length)
            pos = length - 1
            arghmm.arghmm_sample_arg_removal_path2(trees, node, pos, path)

            for j in xrange(ntrees):
                x.append(j)
                y.append(path[j])
        
        arghmm.delete_local_trees(trees)

        p = plot(dither(x, .4), dither(y, .4))

        pause()


    def test_sample_removal_path_recomb(self):
        """
        Test the sampling of a branch removal path
        """

        k = 12
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(500e3) / 20
        times = arghmm.get_time_points(ntimes=20)

        util.tic("sim")
        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        trees, names = arghmm.arg2ctrees(arg,times)
        nnodes = 2*k - 1
        ntrees = arghmm.get_local_trees_ntrees(trees)
        print ntrees, nnodes
        util.toc()

        x = []
        y = []
        for i in range(1):
            path = [0] * ntrees
            arghmm.arghmm_sample_arg_removal_path_recomb(trees, .5, path)

            for j in xrange(ntrees):
                x.append(j)
                y.append(path[j])
        
        arghmm.delete_local_trees(trees)

        p = plot(dither(x, .4), dither(y, .4))

        pause()




    def test_remove_thread_path(self):
        """
        Test the sampling of a branch removal path
        """

        k = 8
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(200e3) / 20
        times = arghmm.get_time_points(ntimes=20)

        util.tic("sim")
        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        trees, names = arghmm.arg2ctrees(arg,times)
        nnodes = 2*k - 1
        ntrees = arghmm.get_local_trees_ntrees(trees)
        print ntrees, nnodes
        util.toc()


        path = [0] * ntrees
        node = random.randint(0, nnodes - 1)
        arghmm.arghmm_sample_arg_removal_path(trees, node, path)
        arghmm.arghmm_remove_arg_thread_path(trees, path, len(times)+1)


        arghmm.delete_local_trees(trees)


    def test_remove_thread_path2(self):
        """
        Test the sampling of a branch removal path
        """

        k = 12
        n = 1e4
        rho = 1.5e-8 * 20 / 5
        mu = 2.5e-8 * 20 * 5
        length = int(500e3) / 20
        times = arghmm.get_time_points(ntimes=80)

        util.tic("sim")
        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg, mu, times=times)
        seqs = arglib.make_alignment(arg, muts)

        trees, names = carg = arghmm.arg2ctrees(arg, times)
        nnodes = 2*k - 1
        ntrees = arghmm.get_local_trees_ntrees(trees)
        print ntrees, nnodes
        util.toc()


        path = [0] * ntrees
        original_path = [0] * length
        node = random.randint(0, nnodes - 1)
        pos = random.randint(0, length)
        arghmm.arghmm_sample_arg_removal_path2(trees, node, pos, path)
        arghmm.arghmm_remove_arg_thread_path2(trees, path, len(times)+1,
                                              original_path)


        ntimes = len(times)
        popsizes = [n] * ntimes
        seqs2 = arghmm.seqs2cseqs(seqs, names)
        sample_path = [0] * length
        arghmm.arghmm_sample_arg_thread_internal(
            trees, times, ntimes, popsizes, rho, mu,
            seqs2, len(seqs), seqs.alignlen(), sample_path)


        x = [0] * length
        arghmm.arghmm_get_thread_times(trees, len(times), original_path, x)
        y = [0] * length
        arghmm.arghmm_get_thread_times(trees, len(times), sample_path, y)

        #p = plot(x)
        #p.plot(y)
        x = x[::length/1000]
        y = y[::length/1000]
        p = plot(dither(x, .4), dither(y, .4))
        p.plot([0, len(times)], [0, len(times)], style="lines")

        pause()

        #times2 = arghmm.get_time_points(ntimes=22)
        #arg2 = arghmm.ctrees2arg(trees, names, times2)

        '''
        for (start, end), tree in arglib.iter_tree_tracks(arg):
            arglib.remove_single_lineages(tree)
            draw_tree_names(tree.get_tree(), maxlen=5, minlen=5)
            #tree.get_tree().write()

            last_state = None
            for i in range(start, end):
                if original_path[i] != last_state:
                    last_state = original_path[i]
                    print last_state
        '''



    def test_forward_internal_thread(self):
        """
        Test the sampling of a branch removal path
        """

        k = 12
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(1000e3) / 20
        times = arghmm.get_time_points(ntimes=20)

        util.tic("sim")
        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg, mu, times=times)
        seqs = arglib.make_alignment(arg, muts)
        
        trees, names = carg = arghmm.arg2ctrees(arg, times)
        nnodes = 2*k - 1
        ntrees = arghmm.get_local_trees_ntrees(trees)
        print ntrees, nnodes
        util.toc()


        path = [0] * ntrees
        node = random.randint(0, nnodes - 1)
        arghmm.arghmm_sample_arg_removal_path(trees, node, path)
        arghmm.arghmm_remove_arg_thread_path(trees, path, len(times)+1)


        # run forward algorithm
        arghmm.arghmm_forward_algorithm(carg, seqs, rho=rho,
                                        mu=mu, popsizes=n, times=times,
                                        verbose=True, internal=True)
        
        arghmm.delete_local_trees(trees)


    def test_resample_all_arg(self):
        """
        Test the sampling of a branch removal path
        """

        k = 8
        n = 1e4
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(200e3) / 20
        times = arghmm.get_time_points(ntimes=20)
        refine = 20

        util.tic("sim")
        arg = arghmm.sample_arg_dsmc(k, 2*n, rho, start=0, end=length,
                                     times=times)
        muts = arghmm.sample_arg_mutations(arg, mu, times=times)
        seqs = arglib.make_alignment(arg, muts)

        arg = arghmm.sample_arg(seqs, rho=rho, mu=mu, popsizes=n,
                                times=times, verbose=True, carg=True)
        arg = arghmm.resample_all_arg(arg, seqs, rho=rho,
                                      mu=mu, popsizes=n,
                                      refine=refine, times=times, verbose=True,
                                      carg=False)



#=============================================================================
if __name__ == "__main__":

    test_main()




import unittest, random

import arghmm

from rasmus.common import *
from rasmus import stats, hmm
from rasmus.testing import *

from compbio import coal, arglib, fasta

#=============================================================================

def mle_popsize_coal_times(k, times, mintime):
    s = 0
    i = k
    last = 0
    for t in times:
        s += i*(i-1) * max(t - last, mintime)
        i -= 1
        last = t
    return s / float(4 * k - 4)


def mle_popsize_tree(tree, mintime):
    times = sorted([x.age for x in tree if not x.is_leaf()])
    #print ilen(tree.leaves()), times
    return mle_popsize_coal_times(ilen(tree.leaves()), times, mintime)


#=============================================================================

class Popsize (unittest.TestCase):

    def test_est_popsize(self):
        """
        Fully sample an ARG from stratch using API
        """

        k = 50
        rho = 1.5e-8
        mu = 2.5e-8
        length = int(1e6)
        times = arghmm.get_time_points(ntimes=30, maxtime=200000)
        popsize = 1e4
        refine = 0

        util.tic("sim ARG")
        arg = arghmm.sample_arg_dsmc(k, 2 * popsize,
                                     rho, start=0, end=length, times=times)
        #arg = arglib.sample_arg_smc(k, 2 * popsize,
        #                            rho, start=0, end=length)
        #arg = arglib.sample_arg(k, 2 * popsize, rho, start=0, end=length)
        util.toc()

        x = []
        for tree in arglib.iter_marginal_trees(arg):
            arglib.remove_single_lineages(tree)
            x.append(mle_popsize_tree(tree, mintime=0))
        
        p = plot(x, ymin=0)
        p.plot([0, len(x)], [popsize, popsize], style='lines')
        
        pause()


    def test_est_popsize2(self):
        """
        Fully sample an ARG from stratch using API
        """

        k = 20
        rho = 1.5e-8
        mu = 2.5e-8
        length = int(4e6)
        popsize = 1e4
        popsize2 = 1e4 * .5
        a = int(.3 * length)
        b = int(.7 * length)
        refine = 0

        util.tic("sim ARG")
        arg = arglib.sample_arg_smc(k, 2 * popsize,
                                    rho, start=0, end=a)
        arg = arglib.sample_arg_smc(k, 2 * popsize2,
                                    rho, start=a, end=b,
                                    init_tree=arg)
        arg = arglib.sample_arg_smc(k, 2 * popsize,
                                    rho, start=b, end=length,
                                    init_tree=arg)

        util.toc()

        x = []; y = []
        for (start, end), tree in arglib.iter_tree_tracks(arg):
            arglib.remove_single_lineages(tree)
            x.append(start)
            y.append(mle_popsize_tree(tree, mintime=0))

        x2, y2 = stats.smooth2(x, y, 100e3)
        p = plot(x, y, ymin=0)
        p.plot(x2, y2, style='lines')
        p.plot([0, a, a, b, b, length],
               [popsize, popsize, popsize2, popsize2, popsize, popsize],
               style='lines')
        
        pause()


    def test_est_arg_popsize(self):
        """
        Fully sample an ARG from stratch using API
        """

        k = 20
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(2e6) / 20
        times = arghmm.get_time_points(ntimes=20, maxtime=200000)
        popsize = 1e4
        popsize2 = 1e4 * .5
        a = int(.3 * length)
        b = int(.7 * length)
        refine = 0

        util.tic("sim ARG")
        arg = arglib.sample_arg_smc(k, 2 * popsize,
                                    rho, start=0, end=a)
        arg = arglib.sample_arg_smc(k, 2 * popsize2,
                                    rho, start=a, end=b,
                                    init_tree=arg)
        arg = arglib.sample_arg_smc(k, 2 * popsize,
                                    rho, start=b, end=length,
                                    init_tree=arg)

        # sim seq
        mut = arghmm.sample_arg_mutations(arg, mu, times)
        seqs = arghmm.make_alignment(arg, mut)
        util.toc()

        # sample arg
        util.tic("sample arg")
        arg2 = arghmm.sample_arg(seqs, rho=rho, mu=mu, times=times,
                                 popsizes=1e4, carg=True)
        arg2 = arghmm.resample_climb_arg(arg2, seqs, popsizes=1e4, 
                                         rho=rho, mu=mu, times=times,
                                         refine=200)
        arg2 = arghmm.resample_all_arg(arg2, seqs, popsizes=1e4, 
                                       rho=rho, mu=mu, times=times,
                                       refine=200)
        util.toc()

        x = []; y = []
        for (start, end), tree in arglib.iter_tree_tracks(arg2):
            arglib.remove_single_lineages(tree)
            x.append(start)
            y.append(mle_popsize_tree(tree, mintime=0))

        # thin popsizes
        x2 = range(0, length, length//5000); y2 = []
        j = 0
        for i in range(len(x2)):
            while j < len(x) and x[j] < x2[i]:
                j += 1
            y2.append(y[min(j, len(y)-1)])

        x3, y3 = stats.smooth2(x2, y2, 100e3)
        p = plot(x, y, ymin=0)
        p.plot(x3, y3, style='lines')
        p.plot([0, a, a, b, b, length],
               [popsize, popsize, popsize2, popsize2, popsize, popsize],
               style='lines')
        
        pause()



    def test_popsizes_over_time(self):
        """
        Fully sample an ARG from stratch using API
        """

        k = 20
        rho = 1.5e-8 * 20
        mu = 2.5e-8 * 20
        length = int(1e6) / 20
        times = arghmm.get_time_points(ntimes=30, maxtime=160000)
        a = 60.
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

        util.tic("estimate popsizes")
        popsizes2 = arghmm.est_arg_popsizes(arg, times=times)
        util.toc()
        
        print popsizes2
        p = plot(times, popsizes, xlog=10, xmin=10, ymin=0, ymax=20000)
        p.plot(times[1:], popsizes2)
        
        pause()


    def test_sample_arg_popsizes_trees(self):
        """
        Fully sample an ARG from stratch using API
        """

        k = 2
        rho = 1.5e-8
        mu = 2.5e-8
        length = int(20e6)
        times = arghmm.get_time_points(ntimes=30, maxtime=160000)
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
        length = int(10e6) / 20
        times = arghmm.get_time_points(ntimes=20, maxtime=160000)
        popsizes = [1e4 * (61.-i)/60. for i in range(len(times))]
        refine = 5

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

#=============================================================================
if __name__ == "__main__":

    test_main()


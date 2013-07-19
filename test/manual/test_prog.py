
import unittest

import arghmm

from compbio import arglib
from compbio import coal
from compbio import fasta
from rasmus.common import *
from rasmus.testing import *
rplot_set_viewer("open")


def sites_split(names, col):
    part1 = []
    part2 = []

    c = col[0]
    for i in range(len(col)):
        if col[i] == c:
            part1.append(names[i])
        else:
            part2.append(names[i])
    return min([part1, part2], key=len)


def show_plots(arg_file, sites_file, stats_file, output_prefix,
               rho, mu, popsize):
    """
    Show plots of convergence.
    """

    # read true arg and seqs
    times = arghmm.get_time_points(ntimes=20, maxtime=200000)
    arg = arglib.read_arg(arg_file)
    arghmm.discretize_arg(arg, times, ignore_top=False, round_age="closer")
    arg = arglib.smcify_arg(arg)
    seqs = arghmm.sites2seqs(arghmm.read_sites(sites_file))

    # compute true stats
    arglen = arglib.arglen(arg)
    arg = arghmm.arg2ctrees(arg, times)
    nrecombs = arghmm.get_local_trees_ntrees(arg[0]) - 1
    lk = arghmm.calc_likelihood(
        arg, seqs, mu=mu, times=times,
        delete_arg=False)
    prior = arghmm.calc_prior_prob(
        arg, rho=rho, times=times, popsizes=popsize,
                        delete_arg=False)
    joint = lk + prior

    data = read_table(stats_file)

    # joint
    y2 = joint
    y = data.cget("joint")
    rplot_start(output_prefix + ".trace.joint.pdf", width=8, height=5)
    rp.plot(y, t="l", ylim=[min(min(y), y2), max(max(y), y2)],
            main="joint probability",
            xlab="iterations",
            ylab="joint probability")
    rp.lines([0, len(y)], [y2, y2], col="gray")
    rplot_end(True)

    # lk
    y2 = lk
    y = data.cget("likelihood")
    rplot_start(output_prefix + ".trace.lk.pdf", width=8, height=5)
    rp.plot(y, t="l", ylim=[min(min(y), y2), max(max(y), y2)],
            main="likelihood",
            xlab="iterations",
            ylab="likelihood")
    rp.lines([0, len(y)], [y2, y2], col="gray")
    rplot_end(True)

    # prior
    y2 = prior
    y = data.cget("prior")
    rplot_start(output_prefix + ".trace.prior.pdf", width=8, height=5)
    rp.plot(y, t="l", ylim=[min(min(y), y2), max(max(y), y2)],
            main="prior probability",
            xlab="iterations",
            ylab="prior probability")
    rp.lines([0, len(y)], [y2, y2], col="gray")
    rplot_end(True)

    # nrecombs
    y2 = nrecombs
    y = data.cget("recombs")
    rplot_start(output_prefix + ".trace.nrecombs.pdf",
                width=8, height=5)
    rp.plot(y, t="l", ylim=[min(min(y), y2), max(max(y), y2)],
            main="number of recombinations",
            xlab="iterations",
            ylab="number of recombinations")
    rp.lines([0, len(y)], [y2, y2], col="gray")
    rplot_end(True)

    # arglen
    y2 = arglen
    y = data.cget("arglen")
    rplot_start(output_prefix + ".trace.arglen.pdf",
                width=8, height=5)
    rp.plot(y, t="l", ylim=[min(min(y), y2), max(max(y), y2)],
            main="ARG branch length",
            xlab="iterations",
            ylab="ARG branch length")
    rp.lines([0, len(y)], [y2, y2], col="gray")
    rplot_end(True)


#=============================================================================

class Prog (unittest.TestCase):

    def test_prog(self):

        popsize = 1e4
        mu = 2.20e-8
        rho = 1.16e-8

        #if not os.path.exists("test/data/test_prog/0.sites"):
        if 1:
            make_clean_dir("test/data/test_prog")
            os.system("""bin/arg-sim \
                -k 8 -L 200000 --model dsmc \
                -N 1e4 -r 0.5e-8 -m 2.20e-8 \
                --ntimes 10 --maxtime 200e3 \
                -o test/data/test_prog/0""")

        if 1:
            make_clean_dir("test/data/test_prog/0.sample")
            os.system("""bin/arg-sample \
                -s test/data/test_prog/0.sites \
                -N 1e4 -r 0.5e-8 -m 2.20e-8 \
                --ntimes 20 --maxtime 200e3 -c 10 \
                -n 500 \
                -o test/data/test_prog/0.sample/out""")


        popsize = 1e4
        mu = 2.20e-8
        rho = 0.5e-8

        show_plots(arg_file="test/data/test_prog/0.arg",
                   sites_file="test/data/test_prog/0.sites",
                   stats_file="test/data/test_prog/0.sample/out.stats",
                   output_prefix="test/data/test_prog/0",
                   rho=rho, mu=mu, popsize=popsize)


    def test_prog_gibbs(self):

        #if not os.path.exists("test/data/test_prog/0.sites"):
        if 1:
            makedirs("test/data/test_prog")

            os.system("""bin/arg-sim \
            -k 6 -L 100000 --model dsmc \
            -N 1e4 -r 0.5e-8 -m 2.20e-8 --infsites \
            --ntimes 20 --maxtime 200e3 \
            -o test/data/test_prog/0""")

        if 1:
            make_clean_dir("test/data/test_prog/0.sample")
            os.system("""bin/arg-sample -q \
                -s test/data/test_prog/0.sites \
                -N 1e4 -r 0.5e-8 -m 2.20e-8 --infsites \
                --ntimes 20 --maxtime 200e3 -c 20 \
                -n 5000 --gibbs \
                -o test/data/test_prog/0.sample/out""")

        popsize = 1e4
        mu = 2.20e-8
        rho = 0.5e-8

        show_plots(arg_file="test/data/test_prog/0.arg",
                   sites_file="test/data/test_prog/0.sites",
                   stats_file="test/data/test_prog/0.sample/out.stats",
                   output_prefix="test/data/test_prog/0",
                   rho=rho, mu=mu, popsize=popsize)


    def test_lineages(self):

        popsize = 1e4
        mu = 2.5e-8
        rho = 1.5e-8

        if not os.path.exists("test/data/test_lineages/0.sites"):
            make_clean_dir("test/data/test_lineages")
            os.system("""arg-sim \
            -k 50 -L 1000 \
            -N 1e4 -r 1.5e-50 -m 2.5e-6 \
            --ntimes 20 --maxtime 400e3 \
            -o test/data/test_lineages/0""")
            os.system("""arg2smc --ntimes 50 --maxtime 400e3 \
            test/data/test_lineages/0.arg \
            test/data/test_lineages/0.smc""")

        #-a test/data/test_linegaes/0.smc

        make_clean_dir("test/data/test_lineages/0.sample")
        '''
        os.system("""subsites -n6 -s test/data/test_lineages/0.sites \
        > test/data/test_lineages/0.core.sites; \
        arg-sample \
    -s test/data/test_lineages/0.core.sites \
    -x 1 -N 1e4 -r 1.5e-50 -m 2.5e-6 \
    --ntimes 20 --maxtime 400e3 -c 20 \
    --climb 0 -n 100 \
    -o test/data/test_lineages/0.core/out""")
    '''

        #    -a test/data/test_lineages/0.smc \
        os.system("""arg-sample \
    -s test/data/test_lineages/0.sites \
    -x 1 -N 1e4 -r 1.5e-50 -m 2.5e-6 \
    --ntimes 20 --maxtime 400e3 -c 1 \
    --climb 0 -n 500 \
    -o test/data/test_lineages/0.sample/out""")

        data = read_delim(os.popen("""arg-diff-lineages \
        --ntimes 20 --maxtime 400e3 \
        test/data/test_lineages/0.arg \
        test/data/test_lineages/0.sample/out.%d.smc.gz"""), parse=True)
        write_delim("test/data/test_lineages/diff.txt", data)

        x, y = transpose(data)
        p = plot(x, y, style="lines")
        pause()



    def test_prog_infsites(self):

        popsize = 1e4
        mu = 2.5e-8
        rho = 1.5e-8

        if 1:
            make_clean_dir("test/data/test_prog_infsites")

            os.system("""arg-sim \
            -k 40 -L 200000 \
            -N 1e4 -r 1.5e-8 -m 2.5e-8 --infsites \
            --ntimes 20 --maxtime 400e3 \
            -o test/data/test_prog_infsites/0""")

            make_clean_dir("test/data/test_prog_infsites/0.sample")
            os.system("""arg-sample \
    -s test/data/test_prog_infsites/0.sites \
    -N 1e4 -r 1.5e-8 -m 2.5e-8 \
    --ntimes 5 --maxtime 100e3 -c 1 \
    --climb 0 -n 20 --infsites \
    -x 1 \
    -o test/data/test_prog_infsites/0.sample/out""")


        arg = arghmm.read_arg("test/data/test_prog_infsites/0.sample/out.0.smc.gz")
        sites = arghmm.read_sites("test/data/test_prog_infsites/0.sites")
        print "names", sites.names
        print

        noncompats = []
        for block, tree in arglib.iter_local_trees(arg):
            tree = tree.get_tree()
            treelib.remove_single_children(tree)
            phylo.hash_order_tree(tree)
            for pos, col in sites.iter_region(block[0]+1, block[1]+1):
                assert block[0]+1 <= pos <= block[1]+1, (block, pos)
                split = sites_split(sites.names, col)
                node = arglib.split_to_tree_branch(tree, split)
                if node is None:
                    noncompats.append(pos)
                    print "noncompat", block, pos, col
                    print phylo.hash_tree(tree)
                    print tree.leaf_names()
                    print "".join(col[sites.names.index(name)]
                                  for name in tree.leaf_names())
                    print split
                    print
        print "num noncompats", len(noncompats)
        #print histtab(noncompats)[:10]


    def test_prog_many(self):

        popsize = 1e4
        mu = 2.5e-8
        rho = 1.5e-8

        if not os.path.exists("test/data/test_prog_many/0.sites"):
            make_clean_dir("test/data/test_prog_many")
            os.system("""arg-sim \
            -k 200 -L 100000 \
            -N 1e4 -r 1.5e-8 -m 2.5e-8 \
            --ntimes 20 --maxtime 400e3 \
            -o test/data/test_prog_many/0""")

        make_clean_dir("test/data/test_prog_many/0.sample")
        os.system("""arg-sample \
    -s test/data/test_prog_many/0.sites \
    -N 1e4 -r 1.5e-8 -m 2.5e-8 \
    --ntimes 20 --maxtime 400e3 -c 20 \
    --climb 0 -n 100 \
    -x 1 \
    -o test/data/test_prog_many/0.sample/out""")



    def test_prog_resume(self):

        os.system("""arg-sample \
    -s test/data/test_prog_small/0.sites \
    -N 1e4 -r 1.5e-8 -m 2.5e-8 \
    --ntimes 20 --maxtime 400e3 -c 20 \
    --climb 10 -n 40 \
    -x 1 --resume \
    -o test/data/test_prog_small/0.sample/out""")


    def test_prog_mask(self):

        popsize = 1e4
        mu = 2.20e-8
        rho = 1.16e-8

        if not os.path.exists("test/data/test_prog_mask/0.sites"):
            makedirs("test/data/test_prog_mask")

            os.system("""arg-sim \
            -k 12 -L 10000 --model dsmc \
            -N 1e4 -r 1.16e-8 -m 2.20e-8 \
            --ntimes 20 --maxtime 200e3 \
            -o test/data/test_prog_mask/0""")

            mask = [["chr", 1000, 2000],
                    ["chr", 3000, 4000]]
            write_delim("test/data/test_prog_mask/mask.bed", mask)

        make_clean_dir("test/data/test_prog_mask/0.sample")
        os.system("""arg-sample \
    -s test/data/test_prog_mask/0.sites \
    -N 1e4 -r 1.16e-8 -m 2.20e-8 \
    --ntimes 20 --maxtime 200e3 -c 20 \
    -n 10 \
    --maskmap test/data/test_prog_mask/mask.bed \
    -o test/data/test_prog_mask/0.sample/out""")




    def test_prog_map(self):

        popsize = 1e4
        mu = 2.5e-8
        rho = 1.5e-8

        if not os.path.exists("test/data/test_prog_map/0.sites"):
            makedirs("test/data/test_prog_map")

            mutmap = [["chr", 0, 20000, mu],
                      ["chr", 20000, 50000, mu*.8],
                      ["chr", 50000, 120000, mu*.5],
                      ["chr", 120000, 200000, mu*.6]]
            write_delim("test/data/test_prog_map/mut.map.txt", mutmap)
            rmap = [["chr", -1000, 30000, rho],
                    ["chr", 30000, 60000, rho*.05],
                    ["chr", 60000, 100000, rho*.05]]
            write_delim("test/data/test_prog_map/recomb.map.txt", rmap)

            os.system("""arg-sim \
            -k 12 -L 100000 \
            -N 1e4 -r 1.5e-8 -m 2.5e-8 \
            --ntimes 20 --maxtime 400e3 \
            -o test/data/test_prog_map/0""")

        make_clean_dir("test/data/test_prog_map/0.sample")
        os.system("""arg-sample \
    -s test/data/test_prog_map/0.sites \
    -N 1e4 -V 3 \
    --mutmap test/data/test_prog_map/mut.map.txt \
    --recombmap test/data/test_prog_map/recomb.map.txt \
    --ntimes 20 --maxtime 400e3 -c 20 \
    --climb 0 -n 11 \
    -o test/data/test_prog_map/0.sample/out""")



if __name__ == '__main__':
    test_main()

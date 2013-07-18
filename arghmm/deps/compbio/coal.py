"""

   Coalescent methods


A note about population size.  In this code all population sizes N or n are
uncorrected.  If you need to compute a coalescent for a diploid species
you must multiply N by 2 before passing it to any of these functions.

"""

#=============================================================================
# imports

from __future__ import division

# python imports
import itertools
from itertools import chain, izip
from math import *
import random
from collections import defaultdict

# rasmus imports
from rasmus import treelib, stats, util, linked_list

try:
    from rasmus.symbolic import *
except ImportError:
    # only experimental functions need symbolic
    pass

# compbio imports
from . import birthdeath

# import root finder
try:
    from scipy.optimize import brentq
except ImportError:
    def brentq(f, a, b, disp=False):
        return stats.bisect_root(f, a, b)




#=============================================================================
# single coalescent PDFs, CDFs, and sampling functions


def prob_coal(t, k, n):
    """
    Returns the probability density of observing the first coalesce of 'k'
    individuals in a population size of 'n' at generation 't'
    """

    # k choose 2
    k2 = k * (k-1) / 2
    k2n = k2 / n
    return k2n * exp(- k2n * t)


def sample_coal(k, n):
    """
    Returns a sample coalescent time for 'k' individuals in a population 'n'
    """

    # k choose 2
    k2 = k * (k-1) / 2
    k2n = k2 / n
    return random.expovariate(k2n)


def sample_coal_times(k, n):
    """
    Returns a sampling of (k-1) coalescences for 'k' lineages in a
    population of size 'n'.
    """
    times = [0]
    for j in xrange(k, 1, -1):
        times.append(times[-1] + sample_coal(j, n))
    return times[1:]


def prob_coal_counts(a, b, t, n):
    """
    The probabiluty of going from 'a' lineages to 'b' lineages in time 't'
    with population size 'n'
    """

    C = stats.prod((b+y)*(a-y)/(a+y) for y in xrange(b))
    s = exp(-b*(b-1)*t/2.0/n) * C
    for k in xrange(b+1, a+1):
        k1 = k - 1
        C = (b+k1)*(a-k1)/(a+k1)/(b-k) * C
        s += exp(-k*k1*t/2.0/n) * (2*k-1) / (k1+b) * C

    return s / stats.factorial(b)


def prob_coal_counts_slow(a, b, t, n):
    """
    The probability of going from 'a' lineages to 'b' lineages in time 't'
    with population size 'n'

    Implemented more directly, but slower.  Good for testing against.
    """

    s = 0.0
    for k in xrange(b, a+1):
        i = exp(-k*(k-1)*t/2.0/n) * \
            (2*k-1)*(-1)**(k-b) / stats.factorial(b) / \
            stats.factorial(k-b) / (k+b-1) * \
            stats.prod((b+y)*(a-y)/(a+y) for y in xrange(k))
        s += i
    return s


def prob_coal_cond_counts(x, a, b, t, n):
    """
    Returns the probability density of a coalescent happening at time 'x'
    between 'a' lineages conditioned on there being 'b' lineages at time
    't'.  The population size is 'n'.
    """

    lama = -a*(a-1)/2.0/n
    C = stats.prod((b+y)*(a-1-y)/(a-1+y) for y in xrange(b))
    s = exp(-b*(b-1)/2.0/n*(t-x) + lama*x) * C
    for k in xrange(b+1, a):
        k1 = k - 1
        lam = -k*k1/2.0/n
        C = (b+k1)*(a-1-k1)/(a-1+k1)/(b-k) * C
        s += exp(lam*t + (lama-lam)*x) * (2*k-1) / (k1+b) * C

    return s / stats.factorial(b) * (-lama) / prob_coal_counts(a, b, t, n)


def prob_coal_cond_counts_simple(x, a, b, t, n):
    """
    Returns the probability density of a coalescent happening at time 'x'
    between 'a' lineages conditioned on there being 'b' lineages at time
    't'.  The population size is 'n'.
    """

    return (prob_coal_counts(a-1, b, t-x, n) * prob_coal(x, a, n) /
            prob_coal_counts(a, b, t, n))


def cdf_coal_cond_counts(x, a, b, t, n):
    """
    Returns the probability a coalescent happening *before* time 'x'
    between 'a' lineages conditioned on there being 'b' lineages at time
    't'.  The population size is 'n'.
    """

    lama = -a*(a-1)/2.0/n
    C = stats.prod((b+y)*(a-1-y)/(a-1+y) for y in xrange(b))
    c = -b*(b-1)/2.0/n
    s = exp(c*t) * (exp((lama-c)*x)-1.0) / (lama-c) * C
    for k in xrange(b+1, a):
        k1 = k - 1
        lam = -k*k1/2.0/n
        C = (b+k1)*(a-1-k1)/(a-1+k1)/(b-k) * C
        s += exp(lam*t) * (exp((lama-lam)*x) - 1.0) / (lama - lam) \
             * (2*k-1) / (k1+b) * C

    return s / stats.factorial(b) * (-lama) / prob_coal_counts(a, b, t, n)


def sample_coal_cond_counts(a, b, t, n):
    """
    Samples the next coalescent between 'a' lineages in a population size of
    'n', conditioned on there being 'b' lineages at time 't'.
    """

    # this code solves this equation for t
    #   cdf(t) - p = 0
    # where p ~ U(0, 1)

    p = random.random()

    # compute constants
    lama = -a*(a-1)/2.0/n
    C0 = stats.prod((b+y)*(a-1-y)/(a-1+y) for y in xrange(b))
    c = -b*(b-1)/2.0/n
    d = 1.0/stats.factorial(b) * (-lama) / prob_coal_counts(a, b, t, n)


    # CDF(t) - p
    def f(x):
        if x <= 0:
            return x - p
        if x >= t:
            return 1.0 - p + (x - t)

        C = C0
        s = exp(c*t) * (exp((lama-c)*x)-1.0) / (lama-c) * C
        for k in xrange(b+1, a):
            k1 = k - 1
            lam = -k*k1/2.0/n
            C = (b+k1)*(a-1-k1)/(a-1+k1)/(b-k) * C
            s += exp(lam*t) * (exp((lama-lam)*x) - 1.0) / (lama - lam) \
                 * (2*k-1) / (k1+b) * C

        return s * d - p

    return brentq(f, 0.0, t, disp=False)



def prob_mrca(t, k, n):
    """
    Probability density function of the age 't' of the most recent
    common ancestor (MRCA) of 'k' lineages in a population size 'n'
    """

    s = 0.0
    for i in xrange(1, k):
        lam = (i+1) * i / 2.0 / n
        s += lam * exp(- lam * t) * mrca_const(i, 1, k-1)
    return s


def cdf_mrca(t, k, n):
    """
    Cumulative probability density of the age 't' of the most recent common
    ancestor (MRCA) of 'k' lineages in a population size 'n'
    """

    if k == 1:
        return 1.0

    s = 0.0
    for i in xrange(1, k+1):
        lam = i * (i-1) / (2.0 * n)
        p = 1.0
        for y in xrange(1, i):
            p *= (y-k) / (k+y)
        s += exp(-lam * t) * (2*i - 1) * p
    return s


def mrca_const(i, a, b):
    """A constant used in calculating MRCA"""

    # i+1 choose 2
    y = (i+1) * i / 2.0
    prod = 1.0

    for j in xrange(a, b+1):
        if j == i:
            continue
        # j+1 choose 2
        x = (j+1) * j / 2.0
        prod *= x / (x - y)
    return prod



def prob_bounded_coal(t, k, n, T):
    """
    Probability density function of seeing a coalescence at 't' from
    'k' lineages in a population of size 'n' with bounding time 'T'
    """

    if t > T:
        return 0.0

    if k == 2:
        prob_coal(t, k, n)
    return prob_coal(t, k, n) * cdf_mrca(T-t, k-1, n) / \
           cdf_mrca(T, k, n)


def cdf_bounded_coal(t, k, n, T):
    """
    Cumalative density function of seeing a coalescence at 't' from
    'k' lineages in a population of size 'n' with bounding time 'T'
    """
    i = k - 1

    lam_i = (i+1)*i/2.0 / n
    C = [mrca_const(j, 1, i-1) for j in xrange(1, i)]
    A = lam_i / n / cdf_mrca(T, k, n)
    B = sum(C) / lam_i
    F = [C[j-1] * exp(-(j+1)*j/2.0/n * T) / ((j+1)*j/2.0/n - lam_i)
         for j in xrange(1, i)]

    return (lam_i / cdf_mrca(T, k, n) *
            (B * (1-exp(-lam_i * t))
             - sum(F[j-1] * (exp(((j+1)*j/2.0/n - lam_i)*t)-1)
                   for j in xrange(1, i))))


def sample_bounded_coal(k, n, T):
    """
    Sample a coalescent time 't' for 'k' lineages and population 'n'
    on the condition that the MRCA is before 'T'
    """

    # special case
    if k == 2:
        return sample_bounded_coal2(n, T)

    # this code solves this equation for t
    #   cdf(t) - p = 0
    # where p ~ U(0, 1)

    i = k - 1
    p = random.random()

    # compute constants
    lam_i = (i+1)*i/2.0 / n
    C = [mrca_const(j, 1, i-1) for j in xrange(1, i)]
    A = lam_i / cdf_mrca(T, k, n)
    B = sum(C) / lam_i
    F = [C[j-1] * exp(-(j+1)*j/2.0/n * T) / ((j+1)*j/2.0/n - lam_i)
         for j in xrange(1, i)]

    # CDF(t) - p
    def f(t):
        if t <= 0:
            return t - p
        if t >= T:
            return 1.0 - p + (t - T)

        return (A * (B * (1-exp(-lam_i * t))
             - sum(F[j-1] * (exp(((j+1)*j/2.0/n - lam_i)*t)-1)
                   for j in xrange(1, i)))) - p

    return brentq(f, 0.0, T, disp=False)


def sample_bounded_coal2(n, T):
    """
    Sample a coalescent time 't' for 'k=2' lineages and population 'n'
    on the condition that the MRCA is before 'T'
    """

    # sample from a truncated expontial distribution

    # k choose 2
    lam = 1 / n
    p = exp(-lam * T)
    return - log(random.uniform(p, 1.0)) / lam


def sample_bounded_coal_reject(k, n, T):
    """
    Sample a coalescent time 't' for 'k' lineages and population 'n'
    on the condition that the MRCA is before 'T'

    Uses rejection sampling.  It works but is very inefficient.
    """

    i = k - 1
    consts = [mrca_const(j, 1, i-1) for j in xrange(1, i)]
    x = sum(consts)

    while True:
        while True:
            t = sample_coal(k, n)
            if t < T:
                break

        if i == 1:
            return t

        y = sum(mrca_const(j, 1, i-1) * exp(-((j+1) * j / 2.0 / n) * (T - t))
                for j in xrange(1, i))

        r = 1 - y / x

        if random.random() < r:
            return t


def count_lineages_per_branch(tree, recon, stree):
    """
    Returns the count of gene lineages present at each node in the species
    tree 'tree' given a gene tree 'tree' and reconciliation 'recon'
    """

    # init lineage counts
    lineages = {}
    for snode in stree:
        lineages[snode] = [0, 0]

    for node in tree.postorder():
        snode = recon[node]
        if node.is_leaf():
            lineages[snode][0] += 1 # leaf lineage
        else:
            lineages[snode][1] -= 1 # coal

    for snode in stree.postorder():
        if not snode.is_leaf():
            lineages[snode][0] = sum(lineages[x][1] for x in snode.children)
        lineages[snode][1] += lineages[snode][0]

    return lineages


def get_topology_stats(tree, recon, stree):
    """
    The function computes terms necessary for many topology calculations
    """

    # How many gene nodes per species
    nodes_per_species = dict.fromkeys(stree, 0)

    # How many descendent nodes recon to the same species
    descend_nodes = {}

    # iterate through tree
    for node in tree.postorder():
        if len(node.children) > 1:
            nodes_per_species[recon[node]] += 1
            if not node.is_leaf():
                descend_nodes[node] = 1 + sum(descend_nodes.get(child, 0)
                                              for child in node.children
                                              if recon[child] == recon[node])

    return nodes_per_species, descend_nodes



def prob_multicoal_recon_topology(tree, recon, stree, n,
                                  lineages=None, top_stats=None):
    """
    Returns the log probability of a reconciled gene tree ('tree', 'recon')
    from the coalescent model given a species tree 'stree' and
    population sizes 'n'
    """

    popsizes = init_popsizes(stree, n)
    if lineages is None:
        lineages = count_lineages_per_branch(tree, recon, stree)
    if top_stats is None:
        top_stats = get_topology_stats(tree, recon, stree)

    # iterate through species tree branches
    lnp = 0.0 # log probability
    for snode in stree.postorder():
        if snode.parent:
            # non root branch
            a, b = lineages[snode]

            try:
                p = (util.safelog(prob_coal_counts(a, b, snode.dist,
                                                   popsizes[snode.name]))
                     + stats.logfactorial(top_stats[0].get(snode, 0))
                     - log(num_labeled_histories(a, b)))
            except:
                print (a, b, snode.dist, popsizes[snode.name],
                       prob_coal_counts(a, b, snode.dist,
                                        popsizes[snode.name]),
                       )

                raise

            #p = log(prob_coal_counts(a, b, snode.dist,
            #                            popsizes[snode.name]) *
            #           stats.factorial(top_stats[0].get(snode, 0))
            #           / num_labeled_histories(a, b))
            lnp += p
        else:
            a = lineages[snode][0]
            lnp += (stats.logfactorial(top_stats[0].get(snode, 0)) -
                    log(num_labeled_histories(a, 1)))

    for node, cnt in top_stats[1].iteritems():
        lnp -= log(cnt)

    return lnp



def cdf_mrca_bounded_multicoal(gene_counts, T, stree, n,
                               sroot=None, sleaves=None, stimes=None,
                               tree=None, recon=None):
    """
    What is the log probability that multispecies coalescent in species
    tree 'stree' with population sizes 'n' and extant gene counts 'gene_counts'
    will have a MRCA that occurs in branch 'sroot' before time 'T'.

    As a convenience, you can pass None for gene_counts and give a reconciled
    gene tree instead ('tree', 'recon').
    """

    # determine active part of species tree
    if sroot is None:
        sroot = stree.root
    if sleaves is None:
        sleaves = set(sroot.leaves())

    if len(sleaves) == 0:
        return 0.0

    # init gene counts
    if gene_counts is None:
        if tree is None:
            gene_counts = dict.fromkeys([x.name for x in sleaves], 1)
        else:
            gene_counts = dict.fromkeys([x.name for x in sleaves], 0)
            for leaf in tree.leaves():
                gene_counts[recon[leaf].name] += 1

    popsizes = init_popsizes(stree, n)

    # get time to MRCA above sroot
    if stimes is None:
        stimes = treelib.get_tree_timestamps(stree, sroot, sleaves)

    # use dynamic programming to calc prob of lineage counts
    prob_counts = calc_prob_counts_table(gene_counts, T, stree, popsizes,
                                         sroot, sleaves, stimes)
    return util.safelog(prob_counts[sroot][1][1])


def calc_prob_counts_table(gene_counts, T, stree, popsizes,
                            sroot, sleaves, stimes):

    # use dynamic programming to calc prob of lineage counts
    # format: prob_counts[node] = [a, b]
    prob_counts = {}
    def walk(node):
        if node in sleaves:
            # leaf case
            M = gene_counts[node.name]

            # populate starting lineage counts
            start = [0.0] * (M+1)
            start[M] = 1.0

        elif len(node.children) == 2:
            # internal node case with 2 children

            c1 = node.children[0]
            c2 = node.children[1]
            M1 = walk(c1)
            M2 = walk(c2)
            M = M1 + M2 # max lineage counts in this snode
            end1 = prob_counts[c1][1]
            end2 = prob_counts[c2][1]

            # populate starting lineage counts
            start = [0.0, 0.0]
            for k in xrange(2, M+1):
                start.append(sum(end1[i] * end2[k-i]
                                 for i in xrange(1, k)
                                 if i <= M1 and k-i <= M2))

        elif len(node.children) == 1:
            # single child case

            c1 = node.children[0]
            M1 = walk(c1)
            M = M1 # max lineage counts in this snode
            end1 = prob_counts[c1][1]

            # populate starting lineage counts with child's ending counts
            start = [0.0]
            for k in xrange(1, M+1):
                start.append(end1[k])

        else:
            # unhandled case
            raise Exception("not implemented")



        # populate ending lineage counts
        n = popsizes[node.name]
        ptime = stimes[node.parent] if node.parent else T
        if ptime is None:
            # unbounded end time, i.e. complete coalescence
            end = [0.0, 1.0] + [0.0] * (M-1)
        else:
            # fixed end time
            t = ptime - stimes[node]

            end = [0.0]
            for k in xrange(1, M+1):
                end.append(
                    sum(prob_coal_counts(i, k, t, n) * start[i]
                        for i in xrange(k, M+1)))

        prob_counts[node] = [start, end]

        assert abs(sum(start) - 1.0) < .001, (start, node.children)

        return M
    M = walk(sroot)

    return prob_counts


def prob_coal_bmc(t, u, utime, ucount, gene_counts, T, stree, n,
                  sroot=None, sleaves=None, stimes=None,
                  tree=None, recon=None):
    """
    The PDF of the waiting time 't' for the next coalescent event in species
    branch 'u' within a bounded multispecies coalescent (BMC) process.
    """

    # NOTE: not implemented efficiently

    if sroot is None:
        sroot = stree.root

    # find relevent leaves of stree (u should be treated as a leaf)
    if sleaves is None:
        sleaves = set()
        def walk(node):
            if node.is_leaf() or node == u:
                sleaves.add(node)
            else:
                for child in node.children:
                    walk(child)
        walk(sroot)


    # find timestamps of stree nodes
    if stimes is None:
        # modify timestamp of u to be that of the previous coal (utime)
        stimes = {u : utime}
        stimes = treelib.get_tree_timestamps(stree, sroot, sleaves, stimes)

    # init gene counts
    if gene_counts is None:
        if tree is None:
            gene_counts = dict.fromkeys([x.name for x in sleaves], 1)
        else:
            gene_counts = dict.fromkeys([x.name for x in sleaves], 0)
            for leaf in tree.leaves():
                gene_counts[recon[leaf].name] += 1

        # modify gene counts for species u
        gene_counts[u.name] = ucount

    popsizes = init_popsizes(stree, n)


    p = cdf_mrca_bounded_multicoal(gene_counts, T, stree, popsizes,
                                   sroot=sroot, sleaves=sleaves,
                                   stimes=stimes, tree=tree, recon=recon)

    gene_counts[u.name] = ucount - 1
    stimes[u] = utime + t

    p2 = cdf_mrca_bounded_multicoal(gene_counts, T, stree, popsizes,
                                    sroot=sroot, sleaves=sleaves,
                                    stimes=stimes, tree=tree, recon=recon)

    gene_counts[u.parent.name] = ucount
    stimes[u] = stimes[u.parent]

    p3 = cdf_mrca_bounded_multicoal(gene_counts, T, stree, popsizes,
                                    sroot=sroot, sleaves=sleaves,
                                    stimes=stimes, tree=tree, recon=recon)

    p4 = log(prob_coal(t, ucount, popsizes[u.name]))

    p5 = log(prob_coal_counts(ucount, ucount,
                              stimes[u.parent] - utime, popsizes[u.name]))

    return (p2 + p4) - stats.logsub(p, p3 + p5)


def prob_no_coal_bmc(u, utime, ucount, gene_counts, T, stree, n,
                     sroot=None, sleaves=None, stimes=None,
                     tree=None, recon=None):
    """
    Returns the log probability of no coalescent occurring in branch u
    of the species tree during a bounded multispecies coalescent (BMC).
    """

    if sroot is None:
        sroot = stree.root

    # find relevent leaves of stree (u should be treated as a leaf)
    if sleaves is None:
        sleaves = set()
        def walk(node):
            if node.is_leaf() or node == u:
                sleaves.add(node)
            else:
                for child in node.children:
                    walk(child)
        walk(sroot)


    # find timestamps of stree nodes
    if stimes is None:
        # modify timestamp of u to be that of the previous coal (utime)
        stimes = {u : utime}
        stimes = treelib.get_tree_timestamps(stree, sroot, sleaves, stimes)

    # init gene counts
    if gene_counts is None:
        if tree is None:
            gene_counts = dict.fromkeys([x.name for x in sleaves], 1)
        else:
            gene_counts = dict.fromkeys([x.name for x in sleaves], 0)
            for leaf in tree.leaves():
                gene_counts[recon[leaf].name] += 1

        # modify gene counts for species u
        gene_counts[u.name] = ucount

    popsizes = init_popsizes(stree, n)


    p = cdf_mrca_bounded_multicoal(gene_counts, T, stree, popsizes,
                                   sroot=sroot, sleaves=sleaves, stimes=stimes,
                                   tree=tree, recon=recon)

    gene_counts[u.parent.name] = ucount
    stimes[u] = stimes[u.parent]

    p2 = cdf_mrca_bounded_multicoal(gene_counts, T, stree, popsizes,
                                    sroot=sroot, sleaves=sleaves,
                                    stimes=stimes, tree=tree, recon=recon)

    p3 = log(prob_coal_counts(ucount, ucount,
                              stimes[u.parent] - utime, popsizes[u.name]))

    return p2 - p + p3



def num_labeled_histories(nleaves, nroots):
    n = 1.0
    for i in xrange(nroots + 1, nleaves + 1):
        n *= i * (i - 1) / 2.0
    return n

def log_num_labeled_histories(nleaves, nroots):
    n = 0.0
    for i in xrange(nroots + 1, nleaves + 1):
        n += log(i * (i - 1) / 2.0)
    return n


def prob_bounded_multicoal_recon_topology(tree, recon, stree, n, T,
                                          root=None, leaves=None,
                                          lineages=None, top_stats=None,
                                          stimes=None):
    """
    Returns the log probability of a reconciled gene tree ('tree', 'recon')
    from the coalescent model given a species tree 'stree' and
    population sizes 'n' and stopping time 'T'
    """

    # get input stats
    popsizes = init_popsizes(stree, n)
    if lineages is None:
        lineages = count_lineages_per_branch(tree, recon, stree)
    if top_stats is None:
        top_stats = get_topology_stats(tree, recon, stree)
    if stimes is None:
        stimes = treelib.get_tree_timestamps(stree)


    p = prob_multicoal_recon_topology(tree, recon, stree, popsizes,
                                      lineages=lineages, top_stats=top_stats)
    k_root = lineages[stree.root][0]
    T_root = T - stimes[stree.root]
    return log(cdf_mrca(T_root, k_root, popsizes[recon[tree.root].name])) + p \
           - cdf_mrca_bounded_multicoal(None, T, stree, popsizes,
                                      tree=tree, recon=recon, stimes=stimes)



#=============================================================================
# sampling coalescent trees
#
#  - normal kingman coalescent
#  - censored coalescent
#  - bounded coalescent (conditioned on completion before a fixed time)
#


def sample_coal_tree(k, n):
    """
    Returns a simulated coalescent tree for 'k' leaves from a population 'n'.
    """
    times = [0]
    for j in xrange(k, 1, -1):
        times.append(times[-1] + sample_coal(j, n))
    return make_tree_from_times(times)[0]


def sample_bounded_coal_tree(k, n, T, capped=False):
    """
    Returns a simulated coalescent tree for 'k' leaves from a populations 'n'
    with fixed maximum time 't'.  The simulation is conditioned on returning
    a tree that completely coaleces before time 'T'.

    capped -- if True an artificial root to the tree.  Used primarily by
              other methods.
    """
    times = [0]
    for j in xrange(k, 1, -1):
        times.append(times[-1] + sample_bounded_coal(j, n, T - times[-1]))
    return make_tree_from_times(times, t=T, capped=capped)[0]


def sample_bounded_coal_tree_reject(k, n, T, capped=False):
    """
    Returns a simulated coalescence tree for k leaves from a populations n
    with fixed maximum time t.  The simulation is conditioned on returning
    a tree that completely coaleces before time T.

    This works, but is very inefficient.  Use sample_coal_tree_bounded
    instead.
    """

    # sample times with rejection sampling
    while True:
        times = [0]
        for j in xrange(k, 1, -1):
            times.append(times[-1] + sample_coal(j, n))
        if times[-1] < t:
            break

    return make_tree_from_times(times, t=T, capped=capped)[0]


def sample_censored_coal_tree(k, n, t, capped=False):
    """
    Returns a simulated coalescence tree for 'k' leaves from a population size
    'n' with a fixed maximum time 't'.

    The return value is the tuple (tree, lineages) where lineages is a set
    of lineages that have not yet coalesced.

    capped -- if True, remaining lineages are added as children to a artificial
              tree root.
    """

    times = [0]
    for j in xrange(k, 1, -1):
        times.append(times[-1] + sample_coal(j, n))
        if times[-1] > t:
            times.pop()
            break

    return make_tree_from_times(times, k, t, capped=capped)


def sample_coal_cond_counts_tree(a, b, t, n, capped=False):
    """
    Returns a simulated coalescence tree for 'a' leaves from a population size
    'n', conditioned on their being 'b' lineages at time 't'.

    The return value is the tuple (tree, lineages) where lineages is a set
    of lineages that have not yet coalesced.

    capped -- if True, remaining lineages are added as children to a artificial
              tree root.
    """

    times = [0]
    for j in xrange(a, b, -1):
        times.append(times[-1] + sample_coal_cond_counts(j, b, t-times[-1], n))

    return make_tree_from_times(times, a, t, capped=capped)



def init_popsizes(stree, n):
    """
    Uses 'n' to initialize a population size dict for species tree 'stree'
    """

    if isinstance(n, (int, float)):
        return dict.fromkeys(stree.nodes.keys(), n)
    elif isinstance(n, dict):
        return n
    else:
        raise Exception("n must be a int or dict.")


# TODO: right now this assumes that there are at least 1 or more genes
# in each extant species

def sample_multicoal_tree(stree, n, leaf_counts=None,
                          namefunc=None, sroot=None, sleaves=None):
    """
    Returns a gene tree from a multi-species coalescence process

    stree       -- species tree
    n           -- population size (int or dict)
                   If n is a dict it must map from species name to
                   population size.
    leaf_counts -- dict of species names to a starting gene count.
                   Default is 1 gene per extant species.
    namefunc    -- a function that generates new gene names given a species
                   name.
    """

    if sleaves is None:
        sleaves = set(stree.leaves())
    if sroot is None:
        sroot = stree.root

    # initialize vector for how many genes per extant species
    if leaf_counts is None:
        leaf_counts = dict((l, 1) for l in stree.leaf_names())

    # initialize function for generating new gene names
    if namefunc is None:
        spcounts = dict((l, 1) for l in stree.leaf_names())
        def namefunc(sp):
            name = sp + "_" + str(spcounts[sp])
            spcounts[sp] += 1
            return name

    # initialize population sizes
    popsizes = init_popsizes(stree, n)

    # init gene counts
    counts = dict((n.name, 0) for n in stree)
    counts.update(leaf_counts)

    # init reconciliation
    recon = {}

    # subtrees
    subtrees = {}

    queue = MultiPushQueue(sleaves)

    # loop through species tree
    for snode in queue: #stree.postorder():
        # simulate population for one branch
        k = counts[snode.name]

        if snode != sroot:
            # non basal branch
            queue.push(snode.parent, len(snode.parent.children))
            subtree, lineages = sample_censored_coal_tree(
                k, popsizes[snode.name], snode.dist, capped=True)
            counts[snode.parent.name] += len(lineages)
        else:
            # basal branch
            subtree = sample_coal_tree(k, popsizes[snode.name])
            lineages = [subtree.root]
        subtrees[snode] = (subtree, lineages)

        for node in subtree:
            recon[node] = snode


    # stitch subtrees together
    tree = treelib.Tree()

    # add all nodes to total tree
    for subtree, lineages in subtrees.values():
        tree.merge_names(subtree)
        tree.remove(subtree.root)
        del recon[subtree.root]

    for snode in subtrees:
        if snode not in sleaves: # not snode.is_leaf():
            subtree, lineages = subtrees[snode]

            # get lineages from child subtrees
            lineages2 = chain(*[subtrees[child][1]
                                for child in snode.children])

            # ensure leaves are randomly attached
            leaves = subtree.leaves()
            random.shuffle(leaves)

            # stitch leaves of the subtree to children subtree lineages
            for leaf, lineage in izip(leaves, lineages2):
                tree.add_child(leaf, lineage)


    # set root
    tree.root = subtrees[sroot][0].root
    tree.add(tree.root)
    recon[tree.root] = sroot

    # name leaves
    for node in tree:
        if recon[node].is_leaf():
            tree.rename(node.name, namefunc(recon[node].name))

    #print "HERE"
    #treelib.draw_tree_names(tree, maxlen=8)
    #print "recon", [(x[0].name, x[1].name) for x in recon.items()]


    return tree, recon


def sample_bounded_multicoal_tree(stree, n, T, leaf_counts=None, namefunc=None,
                                  sroot=None, sleaves=None, stimes=None,
                                  gene_counts=None):
    """
    Returns a gene tree from a bounded multi-species coalescence process

    stree       -- species tree
    n           -- population size (int or dict)
                   If n is a dict it must map from species name to
                   population size.
    T           -- deadline for complete coalescence
    leaf_counts -- dict of species names to a starting gene count.
                   Default is 1 gene per extant species.
    namefunc    -- a function that generates new gene names given a species
                   name.
    sleaves     -- you can specify a subtree of the stree by giving the a
                   list 'sleaves' of leaf nodes of the stree
    sroot       -- you can specify a subtree of the stree by giving the
                   subtree root node 'sroot'
    """

    # initialize vector for how many genes per extant species
    if sleaves is None:
        sleaves = set(stree.leaves())
    if sroot is None:
        sroot = stree.root
    if leaf_counts is None:
        leaf_counts = dict((l.name, 1) for l in sleaves)

    # initialize function for generating new gene names
    if namefunc is None:
        spcounts = dict((l.name, 1) for l in sleaves)
        def namefunc(sp):
            name = sp + "_" + str(spcounts[sp])
            spcounts[sp] += 1
            return name

    # initialize population sizes
    popsizes = init_popsizes(stree, n)

    # init gene counts
    if gene_counts is None:
        gene_counts = dict.fromkeys([x.name for x in sleaves], 1)

    # init species tree timestamps
    if stimes is None:
        stimes = treelib.get_tree_timestamps(stree)


    # calc table
    prob_counts = calc_prob_counts_table(gene_counts, T, stree, popsizes,
                                         sroot, sleaves, stimes)

    # init lineage counts
    lineages = {sroot: [None, 1]}
    for node in sleaves:
        lineages[node] = [gene_counts[node.name], None]


    # sample lineage counts
    sample_lineage_counts(sroot, sleaves,
                          popsizes, stimes, T, lineages, prob_counts)

    # sample coal times
    tree, recon = coal_cond_lineage_counts(lineages, sroot, sleaves,
                                           popsizes, stimes, T, namefunc)
    return tree, recon


def sample_lineage_counts(node, leaves,
                          popsizes, stimes, T, lineages, prob_counts):
    """
    Sample lineage counts conditioned on counts at root and leaves
    of species tree
    """

    a, b = lineages[node]
    if node not in leaves:
        if len(node.children) == 2:
            # two child case

            c1 = node.children[0]
            c2 = node.children[1]
            probs1 = prob_counts[c1][1]
            probs2 = prob_counts[c2][1]

            if b is None:
                # special case where no ending count 'b' is conditioned
                k1 = stats.sample(probs1)
                k2 = stats.sample(probs2)
            else:
                # condition on ending count 'b'
                if node.parent:
                    t = stimes[node.parent] - stimes[node]
                else:
                    t = T - stimes[node]
                n = popsizes[node.name]

                reject = 0
                while True:
                    k1 = stats.sample(probs1)
                    k2 = stats.sample(probs2)
                    if random.random() < prob_coal_counts(k1 + k2, b, t, n):
                        # accept
                        break
                    reject += 1

            # set linages counts
            lineages[node][0] = k1 + k2
            if c1 not in lineages:
                lineages[c1] = [None, k1]
            else:
                lineages[c1][1] = k1
            if c2 not in lineages:
                lineages[c2] = [None, k2]
            else:
                lineages[c2][1] = k2

            # recurse
            sample_lineage_counts(c1, leaves,
                                  popsizes, stimes, T, lineages, prob_counts)
            sample_lineage_counts(c2, leaves,
                                  popsizes, stimes, T, lineages, prob_counts)

        elif len(node.children) == 1:
            # single child case

            c1 = node.children[0]
            probs1 = prob_counts[c1][1]

            if b is None:
                # special case where no ending count 'b' is conditioned
                k1 = stats.sample(probs1)
            else:
                # condition on ending count 'b'
                if node.parent:
                    t = stimes[node.parent] - stimes[node]
                else:
                    t = T - stimes[node]
                n = popsizes[node.name]

                reject = 0
                while True:
                    k1 = stats.sample(probs1)
                    if random.random() < prob_coal_counts(k1, b, t, n):
                        # accept
                        break
                    reject += 1

            # set linages counts
            lineages[node][0] = k1
            if c1 not in lineages:
                lineages[c1] = [None, k1]
            else:
                lineages[c1][1] = k1

            # recurse
            sample_lineage_counts(c1, leaves,
                                  popsizes, stimes, T, lineages, prob_counts)

        else:
            # unhandled case
            raise Excepiton("not implemented")


def coal_cond_lineage_counts(lineages, sroot, sleaves, popsizes, stimes, T,
                             namefunc):
    """Sample coalescent times conditioned on lineage counts"""

    # init reconciliation and subtree dicts
    recon = {}
    subtrees = {}
    caps = set()

    # sample coalescent times
    queue = MultiPushQueue(sleaves)

    # loop through species tree
    for snode in queue:
        # simulate population for one branch
        a, b = lineages[snode]

        if snode != sroot:
            t = stimes[snode.parent] - stimes[snode]
            queue.push(snode.parent, len(snode.parent.children))
        else:
            t = T - stimes[snode] if T is not None else None

        if t is None:
            subtree = sample_coal_tree(a, popsizes[snode.name])
            tops = [subtree.root]
        else:
            subtree, tops = sample_coal_cond_counts_tree(
                a, b, t, popsizes[snode.name], capped=True)

            caps.add(subtree.root)

        subtrees[snode] = (subtree, tops)
        for node in subtree:
            recon[node] = snode

    tree = join_subtrees(subtrees, recon, caps, sroot)

    # set name leaves
    for leaf in tree.leaves():
        tree.rename(leaf.name, namefunc(recon[leaf].name))

    return tree, recon


def join_subtrees(subtrees, recon, caps, sroot):
    """Join several subtrees together into one subtree"""

    # stitch subtrees together
    tree = treelib.Tree()

    # add all nodes to total tree
    for snode, (subtree, tops) in subtrees.iteritems():
        tree.merge_names(subtree)

    # remove cap nodes
    for node in caps:
        # remove cap node
        tree.remove(node)
        del recon[node]

    for snode in subtrees:
        subtree, tops = subtrees[snode]

        # get lineages from child subtrees
        lineages2 = list(chain(*[subtrees[child][1]
                                 for child in snode.children]))

        if len(lineages2) == 0:
            # noting to connect
            continue

        # ensure leaves are randomly attached
        leaves = subtree.leaves()
        random.shuffle(leaves)

        # stitch leaves of the subtree to children subtree lineages
        for leaf, lineage in izip(leaves, lineages2):
            tree.add_child(leaf, lineage)


    # set root
    tree.root = subtrees[sroot][0].root
    if tree.root in caps and len(tree.root.children) == 1:
        tree.root = tree.root.children[0]

    return tree




def sample_bounded_multicoal_tree_reject(stree, n, T, leaf_counts=None,
                                         namefunc=None, sleaves=None,
                                         sroot=None):
    """
    Returns a gene tree from a bounded multi-species coalescence process

    stree       -- species tree
    n           -- population size (int or dict)
                   If n is a dict it must map from species name to
                   population size.
    T           -- deadline for complete coalescence
    leaf_counts -- dict of species names to a starting gene count.
                   Default is 1 gene per extant species.
    namefunc    -- a function that generates new gene names given a species
                   name.
    sleaves     -- you can specify a subtree of the stree by giving the a
                   list 'sleaves' of leaf nodes of the stree
    sroot       -- you can specify a subtree of the stree by giving the
                   subtree root node 'sroot'
    """

    # initialize vector for how many genes per extant species
    if sleaves is None:
        sleaves = set(stree.leaves())
    if sroot is None:
        sroot = stree.root
    if leaf_counts is None:
        leaf_counts = dict((l.name, 1) for l in sleaves)

    # initialize function for generating new gene names
    if namefunc is None:
        spcounts = dict((l.name, 1) for l in sleaves)
        def namefunc(sp):
            name = sp + "_" + str(spcounts[sp])
            spcounts[sp] += 1
            return name

    # initialize population sizes
    popsizes = init_popsizes(stree, n)

    reject = 0
    while True:
        queue = MultiPushQueue(sleaves)

        # init gene counts
        counts = dict((n.name, 0) for n in stree)
        counts.update(leaf_counts)

        # init reconciliation
        recon = {}

        # subtrees
        subtrees = {}

        # loop through species tree
        for snode in queue:
            # simulate population for one branch
            k = counts[snode.name]

            if snode != sroot:
                # non basal branch
                subtree, lineages = sample_censored_coal_tree(
                    k, popsizes[snode.name], snode.dist, capped=True)
                queue.push(snode.parent, len(snode.parent.children))
            else:
                # basal branch
                subtree = sample_coal_tree(k, popsizes[snode.name])
                lineages = subtree.root
            subtrees[snode] = (subtree, lineages)
            if snode != sroot:
                counts[snode.parent.name] += len(lineages)
            for node in subtree:
                recon[node] = snode

        # stitch subtrees together
        tree = treelib.Tree()

        # add all nodes to total tree
        for subtree, lineages in subtrees.values():
            tree.merge_names(subtree)
            tree.remove(subtree.root)
            del recon[subtree.root]

        for snode in subtrees:
            if not snode.is_leaf():
                subtree, lineages = subtrees[snode]

                # get lineages from child subtrees
                lineages2 = chain(*[subtrees[child][1]
                                    for child in snode.children])

                # ensure leaves are randomly attached
                leaves = subtree.leaves()
                random.shuffle(leaves)

                # stitch leaves of the subtree to children subtree lineages
                for leaf, lineage in izip(leaves, lineages2):
                    tree.add_child(leaf, lineage)


        # set root
        tree.root = subtrees[sroot][0].root
        tree.add(tree.root)
        recon[tree.root] = sroot


        # reject tree if basal branch goes past deadline
        times = treelib.get_tree_timestamps(tree)
        if times[tree.root] < T:
            break
        else:
            reject += 1
            #print "reject", reject, times[tree.root], T


    # name leaves
    for leaf in tree.leaves():
        tree.rename(leaf.name, namefunc(recon[leaf].name))

    return tree, recon




def make_tree_from_times(times, k=None, t=None, leaves=None, capped=False):
    """
    Returns a Tree from a list of divergence times.

    The topology is choosen by randomly choosing pairs of leaves.
    """

    # initialize k
    if k is None:
        if leaves is not None:
            k = len(leaves)
        else:
            k = len(times)

    tree = treelib.Tree()

    # initialize k children
    if leaves is None:
        children = set(treelib.TreeNode(tree.new_name()) for i in xrange(k))
    else:
        children = set(treelib.TreeNode(name) for name in leaves)
    for child in children:
        tree.add(child)
        child.data["time"] = 0.0

    # perform random merges
    for i in xrange(1, len(times)):
        # make new parent and merge children
        parent = treelib.TreeNode(tree.new_name())
        parent.data["time"] = times[i]
        a, b = random.sample(children, 2)

        tree.add_child(parent, a)
        tree.add_child(parent, b)

        # adjust children set
        children.remove(a)
        children.remove(b)
        children.add(parent)


    # set branch lengths
    for node in tree:
        if not node.parent:
            if t is not None:
                node.dist = t - node.data["time"]
            else:
                node.dist = 0.0
        else:
            node.dist = node.parent.data["time"] - node.data["time"]

    # for convenience cap the tree for easy drawing/manipulation
    if capped:
        tree.make_root()
        for node in children:
            tree.add_child(tree.root, node)
    else:
        # set root
        if len(children) == 1:
            tree.root = list(children)[0]

    # return tree and remaining lineages
    return tree, children



#=============================================================================
# popsize inference

def mle_popsize_coal_times(k, times):
    s = 0
    i = k
    last = 0
    for t in times:
        s += i*(i-1) * (t - last)
        i -= 1
        last = t
    return s / float(2* k - 2)


def mle_popsize_many_coal_times(k, times):

    ntrees = len(times)
    s = 0

    for times2 in times:
        i = k
        last = 0
        for t in times2:
            s += i*(i-1) * (t - last)
            i -= 1
            last = t
    return s / float(2*ntrees*(k - 1))


def mle_popsize_tree(tree):
    timestamps = treelib.get_tree_timestamps(tree)
    times = sorted([timestamps[node] for node in tree.postorder()
                    if len(node.children) == 2])
    k = len(tree.leaves())
    return mle_popsize_coal_times(k, times)


#=============================================================================
# helper data structures

class MultiPushQueue (object):
    """
    A queue that requires multiple pushes before item is queued
    """

    def __init__(self, lst):
        self._lst = linked_list.LinkedList(lst)
        self._count = {}

    def __iter__(self):
        return self

    def push(self, item, needed):
        count = self._count.setdefault(item, 0)

        # must be queued 'needed' times
        if count + 1 == needed:
            self._lst.append(item)
        else:
            self._count[item] += 1

    def next(self):
        if len(self._lst) == 0:
            raise StopIteration
        else:
            return self._lst.pop_front()





#=============================================================================
# allele frequency

def sample_allele_freq(p, n):
    """
    Sample a new allele frequency using starting allele frequency p and
    population size n
    """

    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0

    if p < 0.05:
        return min(float(stats.poissonvariate(p*n))/n, n)
    if p > 0.95:
        return 1.0 - min(float(stats.poissonvariate((1-p)*n))/n, n)

    mu = p * n
    sigma = sqrt(n * p*(1 - p))
    p1 = random.normalvariate(mu, sigma) / n

    if p1 < 0:
        return 0.0
    if p1 > 1:
        return 1.0
    return p1


def freq_CDF(p, N, t, T, k=50):
    """
    Evaluates the CDF derived from Kimura.
    p is initial frequency of the allele in the population
    N is the population size
    t is time (units?)
    T is the upper limit of the CDF (int from 0 to T)
    k is approximation for the upper limit in the (supposed to be) infinite sum
    """
    return freq_CDF-_legs_ends(legendre(1.0-2*p), legendre(1.0-2*T),
      N, t, k=k)


def freq_CDF_legs_noends(leg_r,leg_T,N,t,k=50):
    """
    Evaluates the CDF derived from Kimura using two Legendre polynomials.
    This does not include the probabilities at 0 and 1 (partial CDF).
    leg_r is the legendre_lambda associated with r
    leg_T is the legendre_lambde associated with T (T', really)
    N is the population size
    t is the time elapsed
    k is the upper limit to approximate the infinite sum
    """
    s = 0.0
    expconst = float(t) / 4.0 / N
    for i in xrange(1,k+1):
        newterm = .5 * (leg_r(i-1) - leg_r(i+1))
        newterm *= exp(- i * (i+1) * expconst)
        newterm *= 1 - leg_T(i)
        s += newterm
    return s


def freq_CDF_legs_ends(leg_r,leg_T,N,t,k=50):
    """
    Evaluates the CDF derived from Kimura using two Legendre polynomials.
    This includes the probabilities at 0 and 1 (full CDF).
    leg_r is the legendre_lambda associated with r
    leg_T is the legendre_lambde associated with T (T', really)
    N is the population size
    t is the time elapsed
    k is the upper limit to approximate the infinite sum
    """
    s = prob_fix(1.0-leg_r(True),N,t) # leg_r(True) currently returns p, so this is probability of extinction
    expconst = float(t) / 4.0 / N
    for i in xrange(1,k+1):
        newterm = .5 * (leg_r(i-1) - leg_r(i+1))
        newterm *= exp(- i * (i+1) * expconst)
        newterm *= 1 - leg_T(i)
        s += newterm
    return s if leg_T(True) < 1.0 else s + prob_fix(leg_r(True),N,t) # add fixation probability if T==1


def freq_prob_range(p, N, t, T1, T2, k=50):
    leg_r = legendre(1.0-2*p)
    leg_T1 = legendre(1.0-2*T1)
    leg_T2 = legendre(1.0-2*T2)
    return (freq_CDF_legs_noends(leg_r, leg_T2, N, t, k=k) - \
      freq_CDF_legs_noends(leg_r, leg_T1, N, t, k=k))
    # uses noends because probabilities at 0 and 1 may be determined using other methods


def sample_freq_CDF(p, N, t):
    """
    Takes an allele frequency p, a population size N, and a time period t.
    Samples from the CDF derived from Kimura to get a new allele frequency.
    N.B.: The current version fails sometimes (on some N, t pairs), presumably
     due to errors in freq_CDF_leg.  These need to be fixed.
    """

    # special cases
    if p == 0.0:
        return 0.0
    elif p == 1.0:
        return 1.0
    elif t == 0.0:
        return p

    y = random.random()
    leg_r = legendre(1.0-2*p)
    extinction = prob_fix(1.0-p, N, t) # probability of allele extinction

    if y < extinction:
        return 0.0 # sample an extinction event
    elif y > 1.0 - prob_fix_leg(leg_r, N, t): #prob_fix(p, N, t):
        return 1.0 # sample a fixation event
    else:
        def f(T):
            return freq_CDF_legs_noends(leg_r, legendre(1.0-2*T), N, t) \
              - y + extinction  # trims extinction probability, assures brentq works

        try:
            return brentq(f, 0.0, 1.0, disp=False)
        except:
            print p, N, t
            raise



# new function for determining Legendre polynomial evaluations
def legendre(r):
    """
    Returns a lambda that calculates the Legendre polynomial based on a
     recursive formula (43) from
     http://mathworld.wolfram.com/LegendrePolynomial.html.
    As the value r is constant, results to calls for different n are cached,
     which reduces runtime for repeated calls.
    The old legendre_old(n,r) function below is intractible for n>~10.
    This function can run with n as high as one million in a fraction of a
     second (using isolated calls, so no caching to build higher values of n).
    """
    def cacheleg(i,d):
        if type(i) == bool:
            return (1.0-d[1])/2.0 if i else d[1] # utility function; may need to be removed
        assert (type(i) == int and i >= 0) # if i is not type bool
        m = d['max']
        if i <= m:
            return d[i]
        x = d[1]
        for n in xrange(m+1,i+1):
            d[n] = 1.0 * ( (2*n-1)*x*d[n-1] - (n-1)*d[n-2] ) / n
        d['max'] = i
        return d[i]
    d = {0:1.0, 1:r, 'max':1}
    assert -1.0 <= r and r <= 1.0 # ensure r in reasonable range
    return lambda n: cacheleg(n,d)


def gegenbauer(i, r):
    return ((i * (i+1)) / 2.0 * hypergeo(i+2, 1 - i, 2, (1 - r) / 2.0))

# this should be the fastest gegenbauer method now (21 July 2010)
def gegenbauer2(i, r):
    leg = legendre(r)
    return ((i * (i+1)) / float((2*i+1)*(1-r*r)) *
            (leg(i-1) - leg(i+1)))

def gegenbauer3(n, a, z):

    tot = 0
    for k in xrange(int(n/2)+1):
        tot += ((-1)**k * stats.gamma(n - k + a) / (
                stats.gamma(a) * stats.factorial(k) * stats.factorial(n - 2*k))
                * ((2*z) ** (n - 2*k)))
    return tot


# TODO: determine proper k and esp values
def prob_fix(p, n, t, k=50, esp=0.000001):
    """Probability of fixation"""
    r = 1 - 2*p
    leg = legendre(r)
    prob = p
    for i in xrange(1, k+1):
        term = (.5 * (-1)**i * (leg(i-1) - leg(i+1)) *
                 exp(-t * i * (i+1) / (4 * n)))
        if term != 0.0 and abs(term) < esp:
            return prob + term
        prob += term

    return prob


# added 02 August 2010
# saves information to leg_r
def prob_fix_leg(leg_r, n, t, k=50, esp=0.000001):
    """Probability of fixation"""
    leg = leg_r
    prob = leg(True) # gets p
    for i in xrange(1, k+1):
        term = (.5 * (-1)**i * (leg(i-1) - leg(i+1)) *
                 exp(-t * i * (i+1) / (4 * n)))
        if term != 0.0 and abs(term) < esp:
            return prob + term
        prob += term

    return prob


def hypergeo(a, b, c, z, k=100):
    """Hypergeometric function"""
    terms = [0.0]
    signs = [1.0]
    for i in xrange(1, k+1):
        term = float((i+a-1)*(i+b-1)*z)/(i+c-1)/i
        signs.append(util.sign(term) * signs[-1])
        if term == 0.0:
            break
        terms.append(log(abs(term)) + terms[i-1])
    return sum(s*exp(i) for s, i in zip(signs, terms))


def loghypergeo(a, b, c, z, k=100):
    """
    Hypergeometric function

    Performs computation in log-space
    """
    terms = [0.0]
    signs = [1.0]
    for i in xrange(1, k+1):
        term = float((i+a-1)*(i+b-1)*z)/(i+c-1)/i
        signs.append(util.sign(term) * signs[-1])
        if term == 0.0:
            break
        terms.append(log(abs(term)) + terms[i-1])

    sgn = 1
    tot = -util.INF

    for s, t in zip(signs, terms):
        sgn, tot = stats.logadd_sign(sgn, tot, s, t)
    return sgn, tot


def hypergeo_mult(i, z1, z2, k=100):

     h1 = hypergeo(1-i, i+2, 2, z1, k)
     h2 = hypergeo(1-i, i+2, 2, z2, k)
     return h1 * h2


def freq_pdf(x, p, n, t, k=8):

    if x > 0.5:
        return freq_pdf(1.0-x, 1.0-p, n, t, k)

    q = 1.0 - p
    prob = -util.INF
    sgn = 1
    t4n = t / (4*n)

    for i in xrange(1, k+1):
        #term = (p * q * i * (i+1) * (2*i+1) *
        #        hypergeo(1-i,i+2,2,p) * hypergeo(1-i,i+2,2,x) *
        #        exp(-t * i * (i+1) / (4*n)))

        lcoff = log(p * q * i * (i+1) * (2*i+1))
        s1, h1 = loghypergeo(1-i,i+2,2,p, i+2)
        s2, h2 = loghypergeo(1-i,i+2,2,x, i+2)
        sgn2 = s1 * s2
        term = (lcoff + h1 + h2 - (i * (i+1) * t4n))

        sgn, prob = stats.logadd_sign(sgn, prob, sgn2, term)

    return sgn * exp(prob)




#=============================================================================

if __name__ == "__main__":
    from rasmus.common import plotfunc

    if 0:
        for i in range(5):
            print "P_%d(x) = " % i, legendre_poly(i)
            print


    #========================
    # hypergeo speed

    a, b, c, z, k = 30, 20, 12, .3, 40

    util.tic("hypergeo_fast")
    for i in range(100):
        hypergeo_fast(a, b, c, z, k)
    util.toc()


    util.tic("hypergeo")
    for i in range(100):
        hypergeo(a, b, c, z, k)
    util.toc()

    util.tic("loghypergeo")
    for i in range(100):
        loghypergeo(a, b, c, z, k)
    util.toc()


    if 0:
        p0 = .5
        k=30

        p = plotfunc(lambda x: freq_pdf(x, p0, 1000, 100, k=k),
                     .01, .99, .01, style="lines")
        p.plotfunc(lambda x: freq_pdf(x, p0, 1000, 200, k=k),
                   .01, .99, .01, style="lines")
        p.plotfunc(lambda x: freq_pdf(x, p0, 1000, 500, k=k),
                   .01, .99, .01, style="lines")
        p.plotfunc(lambda x: freq_pdf(x, p0, 1000, 1000, k=k),
                   .01, .99, .01, style="lines")
        p.plotfunc(lambda x: freq_pdf(x, p0, 1000, 2000, k=k),
                   .01, .99, .01, style="lines")
        p.plotfunc(lambda x: freq_pdf(x, p0, 1000, 3000, k=k),
                   .01, .99, .01, style="lines")
        p.enableOutput(True)
        p.replot()

        #p.plotfunc(lambda x: normalPdf(x, (.5, .1135)),
        #           .01, .99, .01, style="lines")



    if 0:
        p0 = .1

        p = plotfunc(lambda x: freq_pdf(x, p0, 1000, 100, k=25),
                     .01, .99, .01, style="lines")
        p.plotfunc(lambda x: freq_pdf(x, p0, 1000, 200, k=25),
                   .01, .99, .01, style="lines")
        p.plotfunc(lambda x: freq_pdf(x, p0, 1000, 500, k=25),
                   .01, .99, .01, style="lines")
        p.plotfunc(lambda x: freq_pdf(x, p0, 1000, 1000, k=25),
                   .01, .99, .01, style="lines")
        p.plotfunc(lambda x: freq_pdf(x, p0, 1000, 2000, k=25),
                   .01, .99, .01, style="lines")
        p.plotfunc(lambda x: freq_pdf(x, p0, 1000, 3000, k=25),
                   .01, .99, .01, style="lines")
        p.enableOutput(True)
        p.replot()

        #p.plotfunc(lambda x: freq_pdf3(x, .5, 1000, 1000/10, k=40),
        #             .01, .99, .01, style="lines")


    if 0:
        p0 = .5
        k=30

        p = plotfunc(lambda x: freq_pdf(x, p0, 1000, 30, k=k),
                     .01, .99, .01, style="lines")
        p.enableOutput(True)
        p.replot()




#=============================================================================
# old versions


# Legendre polynomial
# this function should be depreciated
def legendre_poly(n):

    """ \frac{1}{2^n n!} d^n/dx^n [(x^2 - 1)^n] """

    return simplify(('mult', ('scalar', 1.0 / (2 ** n * stats.factorial(n))),
                    derivate(('power', ('add', ('power', ('var', 'x'),
                                                         ('scalar', 2)),
                                               ('scalar', -1)),
                                       ('scalar', n)),
                             'x', n)))


# this function should be depreciated
def legendre_old(n, r):
    l = simplify(assign_vars(legendre_poly(n), {'x': r}))
    assert l[0] == 'scalar'
    return l[1]


### TODO: distribution problems arise from probability masses at 0 and 1
def freq_CDF_leg_old(leg, N, t, T, k=50):
    """
    Evaluates the CDF derived from Kimura.
    N.B.: Appears to fail sometimes; this needs to be fixed
    leg is a Legendre (lambda) for evaluating the CDF
    N is the population size
    t is time (units?)
    T is the upper limit of the CDF (int from 0 to T)
    k is approximation for the upper limit in the (supposed to be) infinite sum
    """
    def innersum(i, T, j=0, s=0.0, c=1.0):
        if T == 0.0:
            return 1.0
        if j > i:
            return s
        newc = 1.0 if j == 0 else c * (-T) * (i+j) * (i-j+1) / j / j
        return innersum(i, T, j+1, s+newc, newc)
#    if p == 0.0: # none have the allele
#        return 1.0 # all weight is at 0, so CDF is equal to 1
#    if p == 1.0: # all have the allele
#        return 1.0 if T == 1.0 else 0.0
    s = 0.0
    for i in xrange(1,k+1):
        newterm = leg(i-1) - leg(i+1)
        newterm *= exp(- i * (i+1) / 4.0 * t / N)
        newterm *= .5 - .5 * innersum(i,T)
        s += newterm
    return s



def hypergeo_old(a, b, c, z, k=100):
    """Hypergeometric function"""
    terms = [1.0]
    for i in xrange(1, k+1):
        terms.append(float((i+a-1)*(i+b-1)*z)/(i+c-1)/i * terms[i-1])
    return sum(terms)


# this is depreciated; replaced by an equivalent (but faster) gegenbauer method
def gegenbauer2_old(i, r):
    return ((i * (i+1)) / float((2*i+1)*(1-r*r)) *
            (legendre_old(i-1, r) - legendre_old(i+1, r)))


def freq_pdf_old(x, p, n, t, k=8):

    if x > 0.5:
        return freq_pdf2(1.0-x, 1.0-p, n, t, k)

    q = 1.0 - p
    prob = -util.INF
    sgn = 1
    t4n = t / (4*n)

    for i in xrange(1, k+1):
        #term = (p * q * i * (i+1) * (2*i+1) *
        #        hypergeo(1-i,i+2,2,p) * hypergeo(1-i,i+2,2,x) *
        #        exp(-t * i * (i+1) / (4*n)))

        lcoff = log(p * q * i * (i+1) * (2*i+1))
        h1 = hypergeo(1-i,i+2,2,p, i+2)
        h2 = hypergeo(1-i,i+2,2,x, i+2)
        sgn2 = util.sign(h1) * util.sign(h2)

        if sgn2 != 0:
            term = (lcoff + log(abs(h1)) + log(abs(h2)) +
                    (- i * (i+1) * t4n))
            sgn, prob = stats.logadd_sign(sgn, prob, sgn2, term)

    return sgn * exp(prob)



def freq_pdf2(x, p, n, t, k=8):
    r = 1 - 2*p
    z = 1 - 2*x

    prob = 0.0
    for i in xrange(1, k+1):
        term = ((2*i + 1) * (i - r*r) / float(i * (i+1)) *
                gegenbauer(i, r) * gegenbauer(i, z) *
                exp(-t * i * (i+1) / (4*n)))
        print term
        prob += term

    return prob


def freq_pdf3(x, p, n, t, k=8):
    q = 1.0 - p
    prob = 0.0
    for i in xrange(1, k+1):
        term = (p * q * i * (i+1) * (2*i+1) *
                hypergeo(1-i,i+2,2,p,40) * hypergeo(1-i,i+2,2,x,40) *
                exp(-t * i * (i+1) / (4*n)))
        prob += term

    return prob


def freq_pdf4(x, p, n, t, k=8):
    q = 1.0 - p
    prob = 0.0
    for i in xrange(1, k+1):
        term = (p * q * i * (i+1) * (2*i+1) *
                hypergeo_mult(i, p, x, 100) *
                exp(-t * i * (i+1) / (4*n)))
        prob += term

    return prob


def cdf_mrca2(t, k, n):
    """
    Cumulative probability density of the age 't' of the most recent common
    ancestor (MRCA) of 'k' lineages in a population size 'n'
    """

    if k == 1:
        return 1.0

    s = 0.0
    for i in xrange(1, k):
        lam = (i+1) * i / 2.0 / n
        s += (1 - exp(- lam * t)) * mrca_const(i, 1, k-1)
    return s



def prob_multicoal_recon_topology_old(tree, recon, stree, n,
                                  root=None, leaves=None,
                                  lineages=None, top_stats=None):
    """
    Returns the log probability of a reconciled gene tree ('tree', 'recon')
    from the coalescent model given a species tree 'stree' and
    population sizes 'n'

    This definately has a bug, that the current code fixes.
    """

    popsizes = init_popsizes(stree, n)
    if lineages is None:
        lineages = count_lineages_per_branch(tree, recon, stree)
    if top_stats is None:
        top_stats = get_topology_stats(tree, recon, stree)

    # iterate through species tree branches
    lnp = 0.0 # log probability
    for snode in stree.postorder():
        if snode.parent:
            # non root branch
            a, b = lineages[snode]

            lnp += log(prob_coal_counts(a, b, snode.dist,
                                        popsizes[snode.name]))
            lnp -= log(num_labeled_histories(a, b))
        else:
            a = lineages[snode][0]
            lnp -= log(num_labeled_histories(a, 1))


    # correct for topologies H(T)
    # find connected subtrees that are in the same species branch
    subtrees = []
    subtree_root = {}
    for node in tree.preorder():
        if node.parent and recon[node] == recon[node.parent]:
            subtree_root[node] = subtree_root[node.parent]
        else:
            subtrees.append(node)
            subtree_root[node] = node

    # find leaves through recursion
    def walk(node, subtree, leaves):
        if node.is_leaf():
            leaves.append(node)
        elif (subtree_root[node.children[0]] != subtree and
              subtree_root[node.children[1]] != subtree):
            leaves.append(node)
        else:
            for child in node.children:
                walk(child, subtree, leaves)

    # apply correction for each subtree
    for subtree in subtrees:
        leaves = []
        for child in subtree.children:
            walk(subtree, subtree, leaves)
        if len(leaves) > 2:
            lnp += log(birthdeath.num_topology_histories(subtree, leaves))

    return lnp


def calc_prob_counts_table_old(gene_counts, T, stree, popsizes,
                               sroot, sleaves, stimes):

    root_time = T - stimes[sroot]

    # use dynamic programming to calc prob of lineage counts
    prob_counts = {}
    def walk(node):
        if node in sleaves:
            # leaf case
            M = gene_counts[node.name]
            prob_counts[node] = [0.0] * (M+1)
            prob_counts[node][M] = 1.0
            return M

        else:
            # internal node case
            assert len(node.children) == 2
            c1 = node.children[0]
            c2 = node.children[1]
            ptime = stimes[node]
            t1 = ptime - stimes[c1]   # c1.dist
            t2 = ptime - stimes[c2]   # c2.dist
            M1 = walk(c1)
            M2 = walk(c2)
            M = M1 + M2 # max lineage counts in this snode
            n1 = popsizes[c1.name]
            n2 = popsizes[c2.name]

            prob_counts[node] = [0, 0]
            for k in xrange(2, M+1):
                prob_counts[node].append(sum(
                    sum(prob_coal_counts(i, m, t1, n1) *
                        prob_counts[c1][i]
                        for i in xrange(m, M1+1)) *
                    sum(prob_coal_counts(i, k-m, t2, n2) *
                        prob_counts[c2][i]
                        for i in xrange(k-m, M2+1))
                    for m in xrange(1, k)))

            assert abs(sum(prob_counts[node]) - 1.0) < .001
            return M
    M = walk(sroot)

    return prob_counts



def count_lineages_per_branch_old(tree, recon, stree, rev_recon=None):
    """
    Returns the count of gene lineages present at each node in the species
    tree 'tree' given a gene tree 'tree' and reconciliation 'recon'
    """

    # init reverse reconciliation
    if rev_recon is None:
        rev_recon = get_rev_recon(tree, recon, stree)

    # init lineage counts
    lineages = {}
    for snode in stree:
        if snode.is_leaf():
            lineages[snode] = [len([x for x in rev_recon[snode]
                                   if x.is_leaf()]), 0]
        else:
            lineages[snode] = [0, 0]

    # iterate through species tree branches
    for snode in stree.postorder():
        if snode.parent:
            # non root branch
            a = lineages[snode][0]

            # subtract number of coals in branch
            b = a - len([x for x in rev_recon.get(snode, [])
                         if len(x.children) > 1])
            lineages[snode][1] = b
            lineages[snode.parent][0] += b
        else:
            lineages[snode][1] = 1

    lineages2 = count_lineages_per_branch(tree, recon, stree)
    assert lineages == lineages2

    return lineages


def get_rev_recon(tree, recon, stree):
    """
    Returns a reverse reconciliation

    A reverse reconciliation is a mapping from nodes in the species tree to
    lists of nodes in the gene tree.
    """
    rev_recon = {}
    nodes = set(tree.postorder())
    for node, snode in recon.iteritems():
        if node not in nodes:
            raise Exception("node '%s' not in tree" % node.name)
        rev_recon.setdefault(snode, []).append(node)
    return rev_recon


def get_topology_stats_old(tree, recon, stree, rev_recon=None):
    """
    The function computes terms necessary for many topology calculations
    """

    nodes_per_species = {} # How many gene nodes per species
    descend_nodes = {} # How many descendent nodes recon to the same species

    nodes_per_species = dict.fromkeys(stree, 0)

    # init reverse reconciliation
    if rev_recon is None:
        rev_recon = get_rev_recon(tree, recon, stree)

    # iterate through species tree
    for snode, nodes in rev_recon.iteritems():
        nodes_per_species[snode] = len([x for x in nodes
                                        if len(x.children) > 1])

    # iterate through tree
    for node in tree.postorder():
        if not node.is_leaf() and len(node.children) > 1:
            descend_nodes[node] = 1 + sum(descend_nodes.get(child, 0)
                                          for child in node.children
                                          if recon[child] == recon[node])

    return nodes_per_species, descend_nodes


# this is depreciated; replaced by prob_fix method using new legendre method
def prob_fix_old(p, n, t, k=8, esp=0.001):
    """Probability of fixation"""
    r = 1 - 2*p

    prob = p
    for i in xrange(1, k+1):
        term = (.5 * (-1)**i * (legendre_old(i-1, r) - legendre_old(i+1, r)) *
                 exp(-t * i * (i+1) / (4 * n)))
        if term != 0.0 and abs(term) < esp:
            return prob + term
        prob += term

    return prob

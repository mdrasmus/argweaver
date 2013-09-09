from collections import defaultdict
from itertools import izip
from math import exp
from math import log
from math import sqrt
import random

try:
    import scipy.optimize
except ImportError:
    pass

import argweaver

from compbio import arglib
from rasmus import stats
from rasmus import util


#=============================================================================
# numerically stable summations


def kahan_sum(vals):
    tot = 0.0
    c = 0.0
    for val in vals:
        y = val - c
        t = tot + y
        c = (t - tot) - y
        tot = t
    return tot


def safesum(x):

    n = len(x)

    while n > 2:
        x.sort(key=abs)
        y = []
        for i in range(0, n, 2):
            y.append(sum(x[i:i+2]))
        x = y
        n = len(x)

    return sum(x)


def safelogsum(x):

    n = len(x)

    while n > 1:
        x.sort(key=lambda z: z[1])
        y = []
        for i in range(0, n, 2):
            if i+1 < n:
                y.append(stats.logadd_sign(x[i][0], x[i][1],
                                           x[i+1][0], x[i+1][1]))
            else:
                y.append(x[i])
        x = y
        n = len(x)

    return x[0]


#=============================================================================
# coal counts


def prob_coal_counts_slow(a, b, t, n):
    """
    The probabiluty of going from 'a' lineages to 'b' lineages in time 't'
    with population size 'n'

    Implemented more directly, but slower.  Good for testing against.
    """

    s = 0.0
    for k in xrange(b, a+1):
        i = exp(-k*(k-1)*t/2.0/n) * \
            float(2*k-1)*(-1)**(k-b) / stats.factorial(b) / \
            stats.factorial(k-b) / (k+b-1) * \
            stats.prod((b+y)*(a-y)/float(a+y) for y in xrange(k))
        s += i
    return s


def sample_coal_count(a, t, n):
    """
    Sample the number lineages present starting after generations 't' starting
    with 'a' lineages and population size 'n'
    """

    t2 = 0.0
    b = a
    while b > 1:
        rate = b * (b - 1) / 2.0 / n
        t2 += random.expovariate(rate)
        if t2 < t:
            b -= 1
        else:
            break
    return b


def prob_coal_counts(a, b, t, n):
    """
    The probability of going from 'a' lineages to 'b' lineages in time 't'
    with population size 'n'
    """

    try:
        terms = []
        C = stats.prod((b+y)*(a-y)/float(a+y) for y in xrange(b)) \
            / float(stats.factorial(b))
        terms.append(exp(-b*(b-1)*t/2.0/n) * C)
        for k in xrange(b+1, a+1):
            k1 = k - 1
            C = (b+k1)*(a-k1)/float(a+k1)/float(b-k) * C
            terms.append(exp(-k*k1*t/2.0/n) * (2*k-1) / float(k1+b) * C)

        terms.sort(key=abs)
        return kahan_sum(terms)
    except:
        print a, b, t, n
        raise


def smooth_boundary(x, boundary1, boundary2):
    if x < boundary2:
        width = boundary2 - boundary1
        x = width * 1.0 / (boundary2 - x) + boundary1
    return x


def log_prob_coal_counts(a, b, t, n, minprob=1e-200,
                         boundary1=50, boundary2=100):

    # provide a smooth boundary for n
    n = smooth_boundary(n, boundary1, boundary2)

    #p = max(abs(prob_coal_counts(a, b, t, n)), minprob)
    #p = max(abs(prob_coal_counts_frac(a, b, t, n)), minprob)
    p = max(prob_coal_counts(a, b, t, n), minprob)
    return log(p)


def log_prob_many_coal_counts(As, Bs, t, n, Cs=None):
    if Cs is None:
        Cs = [1] * len(As)
    return sum(c * log_prob_coal_counts(a, b, t, n)
               for a, b, c in izip(As, Bs, Cs))


#=============================================================================
# maxmimum likelihood estimation


def mle_prob_coal_counts(a, b, t, n0):
    def f(x):
        if x[0] < 100:
            x[0] = 50 * 1.0 / (101 - n0) + 50
        return - log_prob_coal_counts(a, b, t, x[0])
    return scipy.optimize.fmin(f, n0, disp=False)[0]


def mle_prob_many_coal_counts(As, Bs, t, n0):
    """
    Find the maximum likelihood estimate of going from 'As' lineages
    to 'Bs' lineages over generations 't'.
    """

    # count unique (a,b) pairs
    counts = defaultdict(lambda: 0)
    for a, b in zip(As, Bs):
        counts[(a, b)] += 1
    As = []
    Bs = []
    Cs = []
    for (a, b), c in counts.items():
        As.append(a)
        Bs.append(b)
        Cs.append(c)

    def f(x):
        return - log_prob_many_coal_counts(As, Bs, t, x[0], Cs=Cs)
    return scipy.optimize.fmin(f, n0, disp=False)[0]


#=============================================================================
# trees


def count_tree_lineages(tree, times):
    """
    Counts the lineages present in a tree
    """
    ntimes = len(times)
    starts = []
    ends = []

    # get time steps
    midpoints = [0.0] + [sqrt((times[i+1]+1.0)*(times[i]+1.0))
                         for i in range(ntimes-1)]
    time_steps = [midpoints[i+1] - midpoints[i] for i in range(ntimes-1)]

    # count lineages
    nbranches = argweaver.get_nlineages(tree, times)
    nleaves = util.ilen(tree.leaves())

    for j in range(ntimes-1):
        starts.append(nbranches[j-1] if j > 0 else nleaves)
        ends.append(nbranches[j])

    return starts, ends, time_steps


def count_trees_lineages(trees, times):
    """
    Counts the lineages present in a set of trees
    """
    ntimes = len(times)
    starts = [[] for i in range(ntimes)]
    ends = [[] for i in range(ntimes)]

    for tree in trees:
        s, e, t = count_tree_lineages(tree, times)
        for j in range(ntimes-1):
            starts[j].append(s[j])
            ends[j].append(e[j])
        time_steps = t

    return starts, ends, time_steps


def est_popsize_trees(trees, times, n0=1e4):
    """
    Estimate population size from a set of independent trees
    """
    ntimes = len(times)
    starts, ends, time_steps = count_trees_lineages(trees, times)

    popsize = []
    for j in range(ntimes-1):
        popsize.append(mle_prob_many_coal_counts(
            starts[j], ends[j], time_steps[j], n0))

    return popsize


#=============================================================================
# more sophisticated popsize estimation


class PopsizeEstimator (object):

    def __init__(self, times):

        ntimes = len(times)

        self.times = times
        self.time_steps = [times[i] - times[i - 1] for i in range(1, ntimes)]
        self.ncoals = [0] * ntimes
        self.k_lineages = [0] * ntimes

        self.init_trees = []

        midpoints = [0.0] + [(times[i+1] + times[i]) / 2.0
                             for i in range(ntimes-1)]
        self.time_steps2 = [midpoints[i+1] - midpoints[i]
                            for i in range(ntimes-1)]

    def add_arg(self, arg):

        nleaves = len(list(arg.leaves()))
        times = self.times
        assert times
        eps = 1e-3

        def get_local_children(node, pos, local):
            return set(child for child in arg.get_local_children(node, pos)
                       if child in local)

        def get_parent(node, pos, local):
            parent = arg.get_local_parent(node, pos)
            while len(get_local_children(parent, pos, local)) == 1:
                parent = arg.get_local_parent(parent, pos)
            return parent

        # add initial tree
        tree = arg.get_marginal_tree(arg.start)
        starts, ends, time_steps = count_tree_lineages(tree, times)
        self.init_trees.append({"starts": starts,
                                "ends": ends,
                                "time_steps": time_steps})

        # loop through sprs
        for recomb_pos, (rnode, rtime), (cnode, ctime), local in \
                arglib.iter_arg_sprs(arg, use_local=True):
            i, _ = util.binsearch(times, ctime)
            self.ncoals[i] += 1

            recomb_node = arg[rnode]
            broken_node = get_parent(recomb_node, recomb_pos-eps, local)
            coals = [0.0] + [
                node.age for node in local
                if len(get_local_children(node, recomb_pos-eps, local)) == 2]

            coals.sort()
            nlineages = range(nleaves, 0, -1)
            assert len(nlineages) == len(coals)

            # subtract broken branch
            r = coals.index(recomb_node.age)
            r2 = coals.index(broken_node.age)
            for i in range(r, r2):
                nlineages[i] -= 1

            # get average number of branches in the time interval
            data = zip(coals, nlineages)
            for t in times[1:]:
                data.append((t, "time step"))
            data.sort()

            lineages_per_time = []
            counts = []
            last_lineages = 0
            last_time = 0.0
            for a, b in data:
                if b != "time step":
                    if a > last_time:
                        counts.append((last_lineages, a - last_time))
                    last_lineages = b
                else:
                    counts.append((last_lineages, a - last_time))
                    s = sum(u * v for u, v in counts)
                    total_time = sum(v for u, v in counts)
                    if s == 0.0:
                        lineages_per_time.append(last_lineages)
                    else:
                        lineages_per_time.append(s / total_time)
                    counts = []
                last_time = a

            assert len(lineages_per_time) == len(self.time_steps)

            r, _ = util.binsearch(times, rtime)
            c, _ = util.binsearch(times, ctime)
            for j in range(r, c):
                self.k_lineages[j] += lineages_per_time[j]

    def calc_prob(self, i, n):

        if n < 100:
            n = 50 * 1.0 / (101 - n) + 50
        p = (- self.time_steps2[i] * self.k_lineages[i] / (2.0*n)
             - self.ncoals[i] * log(n))

        p = 0.0
        for init_tree in self.init_trees:
            p += log_prob_coal_counts(init_tree["starts"][i],
                                      init_tree["ends"][i],
                                      self.time_steps2[i], 2*n)

        return p

    def mle_popsize(self, i, n0=1e4):
        def f(n):
            return - self.calc_prob(i, n)
        return scipy.optimize.fmin(f, n0, disp=False)[0]

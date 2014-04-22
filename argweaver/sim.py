#
# Simulating ARG and alignments from DSMC model
#


# python imports
from collections import defaultdict
from math import exp
import random

# rasmus, compbio imports
from compbio import arglib
from compbio import fasta
from rasmus import stats
from rasmus import util

# argweaver imports
import argweaver


def find_region(pos, track):
    """
    Returns the region in 'track' that containts position 'pos'
    """
    for i, region in enumerate(track):
        if pos >= region[1] and pos < region[2]:
            return i, region
    return len(track), None


def sample_tree(k, popsizes, times, start=0, end=1,
                names=None, make_names=True):
    """
    Samples a coalescent tree using multiple population sizes
    """
    ntimes = len(times)
    coal_times = []
    events = []

    timei = 0
    n = popsizes[timei]
    t = 0.0
    k2 = k
    while k2 > 1:
        coal_rate = (k2 * (k2-1) / 2) / float(n)
        t2 = random.expovariate(coal_rate)

        if timei < ntimes-2 and t + t2 > times[timei+1]:
            # advance to next time segment
            timei += 1
            t = times[timei]
            n = popsizes[timei]
            continue

        t += t2
        coal_times.append(t)
        events.append("coal")
        k2 -= 1

    arg = arglib.make_arg_from_times(k, coal_times, events,
                                     start=start, end=end,
                                     names=names, make_names=make_names)
    return arg


def sample_next_recomb(treelen, rho, pos=None, recombmap=None, minlen=1):
    """
    Sample when next recombination will occur along the genome
    """

    if recombmap:
        pos2 = 0
        i, region = find_region(pos, recombmap)
        if not region:
            # no more regions
            return sample_next_recomb(treelen, recombmap[-1][3])

        while True:
            rho = region[3]

            # sample next recomb
            while True:
                pos2 += random.expovariate(max(treelen * rho, rho))
                blocklen = pos2 - pos
                if pos2 < region[2] and blocklen > minlen:
                    return blocklen
                if blocklen >= region[2]:
                    break

            # recomb is not in this block, keep going
            pos2 = region[2]
            i += 1
            if i >= len(recombmap):
                # no more regions
                return sample_next_recomb(treelen, recombmap[-1][3])
            region = recombmap[i]

    else:
        blocklen = 0
        while blocklen < minlen:
            blocklen = random.expovariate(max(treelen * rho, rho))
        return blocklen


def get_coal_times(times):
    # get midpoints
    ntimes = len(times) - 1
    times2 = []
    for i in range(ntimes):
        times2.append(times[i])
        times2.append(((times[i+1]+1)*(times[i]+1))**.5)
    times2.append(times[ntimes])
    return times2


def get_coal_time_steps(times):
    # get midpoints
    ntimes = len(times) - 1
    times2 = []
    for i in range(ntimes):
        times2.append(times[i])
        times2.append(((times[i+1]+1)*(times[i]+1))**.5)
    times2.append(times[ntimes])

    coal_time_steps = []
    for i in range(0, len(times2), 2):
        coal_time_steps.append(times2[min(i+1, len(times2)-1)] -
                               times2[max(i-1, 0)])
    return coal_time_steps


def sample_dsmc_sprs(
        k, popsize, rho, recombmap=None, start=0.0, end=0.0, times=None,
        init_tree=None, names=None, make_names=True):
    """
    Sample ARG using Discrete Sequentially Markovian Coalescent (SMC)

    k          -- chromosomes
    popsize    -- effective population size (haploid)
    rho        -- recombination rate (recombinations / site / generation)
    recombmap  -- map for variable recombination rate
    start      -- staring chromosome coordinate
    end        -- ending chromsome coordinate
    t          -- initial time (default: 0)
    names      -- names to use for leaves (default: None)
    make_names -- make names using strings (default: True)
    """

    assert times is not None
    ntimes = len(times) - 1
    time_steps = [times[i] - times[i-1] for i in range(1, ntimes+1)]
    times2 = get_coal_times(times)

    if hasattr(popsize, "__len__"):
        popsizes = popsize
    else:
        popsizes = [popsize] * len(time_steps)

    # yield initial tree first
    if init_tree is None:
        init_tree = sample_tree(k, popsizes, times, start=start, end=end,
                                names=names, make_names=make_names)
        argweaver.discretize_arg(init_tree, times, ignore_top=True)
    yield init_tree

    # sample SPRs
    pos = start
    tree = init_tree.copy()
    while True:
        # sample next recomb point
        treelen = sum(x.get_dist() for x in tree)
        blocklen = int(sample_next_recomb(treelen, rho, pos=pos,
                                          recombmap=recombmap, minlen=1))
        pos += blocklen
        if pos >= end - 1:
            break

        root_age_index = times.index(tree.root.age)

        # choose time interval for recombination
        states = set(argweaver.iter_coal_states(tree, times))
        nbranches, nrecombs, ncoals = argweaver.get_nlineages_recomb_coal(
            tree, times)
        probs = [nbranches[i] * time_steps[i]
                 for i in range(root_age_index+1)]
        recomb_time_index = stats.sample(probs)
        recomb_time = times[recomb_time_index]

        # choose branch for recombination
        branches = [x for x in states if x[1] == recomb_time_index and
                    x[0] != tree.root.name]
        recomb_node = tree[random.sample(branches, 1)[0][0]]

        # choose coal time
        j = recomb_time_index
        last_kj = nbranches[max(j-1, 0)]
        while j < ntimes - 1:
            kj = nbranches[j]
            if ((recomb_node.name, j) in states and
                    recomb_node.parents[0].age > times[j]):
                kj -= 1
            assert kj > 0, (j, root_age_index, states)

            A = (times2[2*j+1] - times2[2*j]) * kj
            if j > recomb_time_index:
                A += (times2[2*j] - times2[2*j-1]) * last_kj
            coal_prob = 1.0 - exp(-A/float(popsizes[j]))
            if random.random() < coal_prob:
                break
            j += 1
            last_kj = kj
        coal_time_index = j
        coal_time = times[j]

        # choose coal node
        # since coal points collapse, exclude parent node, but allow sibling
        exclude = []

        def walk(node):
            exclude.append(node.name)
            if node.age == coal_time:
                for child in node.children:
                    walk(child)

        walk(recomb_node)
        exclude2 = (recomb_node.parents[0].name,
                    times.index(recomb_node.parents[0].age))
        branches = [x for x in states if x[1] == coal_time_index and
                    x[0] not in exclude and x != exclude2]
        coal_node = tree[random.sample(branches, 1)[0][0]]

        # yield SPR
        rleaves = list(tree.leaf_names(recomb_node))
        cleaves = list(tree.leaf_names(coal_node))

        yield pos, (rleaves, recomb_time), (cleaves, coal_time)

        # apply SPR to local tree
        broken = recomb_node.parents[0]
        recoal = tree.new_node(age=coal_time,
                               children=[recomb_node, coal_node])

        # add recoal node to tree
        recomb_node.parents[0] = recoal
        broken.children.remove(recomb_node)
        if coal_node.parents:
            recoal.parents.append(coal_node.parents[0])
            util.replace(coal_node.parents[0].children, coal_node, recoal)
            coal_node.parents[0] = recoal
        else:
            coal_node.parents.append(recoal)

        # remove broken node
        broken_child = broken.children[0]
        if broken.parents:
            broken_child.parents[0] = broken.parents[0]
            util.replace(broken.parents[0].children, broken, broken_child)
        else:
            broken_child.parents.remove(broken)

        del tree.nodes[broken.name]
        tree.set_root()


def sample_arg_dsmc(k, popsize, rho, recombmap=None,
                    start=0.0, end=0.0, times=None,
                    init_tree=None, names=None, make_names=True):
    """
    Returns an ARG sampled from the Discrete Sequentially Markov Coalescent

    k   -- chromosomes
    popsize -- effective population size
    rho -- recombination rate (recombinations / site / generation)
    recombmap -- map for variable recombination rate
    start -- staring chromosome coordinate
    end   -- ending chromsome coordinate

    names -- names to use for leaves (default: None)
    make_names -- make names using strings (default: True)
    """
    if times is None:
        maxtime = 160000
        delta = .01
        ntimes = 20
        times = argweaver.get_time_points(ntimes, maxtime, delta)

    it = sample_dsmc_sprs(
        k, popsize, rho, recombmap=recombmap,
        start=start, end=end, times=times,
        init_tree=init_tree, names=names, make_names=make_names)
    tree = it.next()
    arg = arglib.make_arg_from_sprs(tree, it)

    return arg


def sample_arg_mutations(arg, mu, times=None):
    """
    Simulate mutations on an ARG.

    Mutations are represented as (node, parent, site, time).

    arg -- ARG on which to simulate mutations
    mu -- mutation rate (mutations/site/gen)
    times -- optional list of discretized time points
    """
    mutations = []
    minlen = times[1] * .1 if times else 0.0

    for (start, end), tree in arglib.iter_local_trees(arg):
        arglib.remove_single_lineages(tree)
        for node in tree:
            if not node.parents:
                continue
            blen = max(node.get_dist(), minlen)
            rate = blen * mu
            i = start
            while i < end:
                i += random.expovariate(rate)
                if i < end:
                    t = random.uniform(node.age, node.age + blen)
                    mutations.append((node, node.parents[0], int(i), t))
    return mutations


def make_sites(arg, mutations, chrom="chr"):
    """
    Make Sites from ARG and sampled mutations.
    """
    leaves = list(arg.leaf_names())
    sites = argweaver.Sites(names=leaves, chrom=chrom,
                            region=[arg.start+1, arg.end])

    # sort mutations by position
    mutations.sort(key=lambda x: x[2])

    for mut_group in util.iter_groups(mutations, key=lambda x: int(x[2])):
        ancestral = "ACGT"[random.randint(0, 3)]
        pos = int(mut_group[0][2])

        # count mutations per branch
        mut_count = defaultdict(int)
        for mut in mut_group:
            mut_count[mut[0].name] += 1

        tree = arg.get_marginal_tree(pos-.5)
        bases = {tree.root.name: ancestral}

        for node in tree.preorder():
            if not node.parents:
                continue

            ancestral = bases[node.parents[0].name]
            if node.name in mut_count:
                c = mut_count[node.name]
                i = 0
                while True:
                    derived = ancestral
                    while derived == ancestral:
                        derived = "ACGT"[random.randint(0, 3)]
                    i += 1
                    if i == c:
                        break
                    ancestral = derived

                bases[node.name] = derived
            else:
                bases[node.name] = ancestral

        col = "".join(bases[l] for l in leaves)
        sites.append(pos+1, col)

    return sites


def make_alignment(arg, mutations, infsites=False):
    """
    Make FASTA alignment from ARG and sampled mutations.
    """
    aln = fasta.FastaDict()
    alnlen = int(arg.end - arg.start)
    leaves = list(arg.leaf_names())
    nleaves = len(leaves)

    # sort mutations by position
    mutations.sort(key=lambda x: x[2])

    # make align matrix
    mat = []

    muti = 0
    for i in xrange(alnlen):
        ancestral = "ACGT"[random.randint(0, 3)]

        if muti >= len(mutations) or i < int(mutations[muti][2]):
            # no mut
            mat.append(ancestral * nleaves)
        else:
            # mut
            mut_count = defaultdict(int)
            while muti < len(mutations) and i == int(mutations[muti][2]):
                mut_count[mutations[muti][0].name] += 1
                muti += 1

            # enforce infinite sites
            if infsites:
                mut_count = {random.sample(mut_count.items(), 1)[0][0]: 1}

            tree = arg.get_marginal_tree(i-.5)
            bases = {tree.root.name: ancestral}

            for node in tree.preorder():
                if not node.parents:
                    continue

                ancestral = bases[node.parents[0].name]
                if node.name in mut_count:
                    c = mut_count[node.name]
                    i = 0
                    while True:
                        derived = ancestral
                        while derived == ancestral:
                            derived = "ACGT"[random.randint(0, 3)]
                        i += 1
                        if i == c:
                            break
                        ancestral = derived

                    bases[node.name] = derived
                else:
                    bases[node.name] = ancestral

            mat.append("".join(bases[l] for l in leaves))

    # make fasta
    for i, leaf in enumerate(leaves):
        aln[leaf] = "".join(x[i] for x in mat)

    return aln

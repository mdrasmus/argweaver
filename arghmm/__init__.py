#
# Ancestral Recombination Graph Hidden Markov Model (ArgHmm)
#

# python libs
from math import exp, log
import random
from itertools import chain, izip
import heapq

# rasmus combio libs
from rasmus import hmm, util, stats, treelib
from rasmus.stats import logadd
from compbio import arglib, fasta, phylo

# arghmm libs
from arghmmc import *
from . import emit
from sample import *


#=============================================================================
# constants

PROGRAM_NAME = u"arghmm"
PROGRAM_VERSION_MAJOR = 0
PROGRAM_VERSION_MINOR = 1
PROGRAM_VERSION_RELEASE = 0
PROGRAM_VERSION = (PROGRAM_VERSION_MAJOR,
                   PROGRAM_VERSION_MINOR,
                   PROGRAM_VERSION_RELEASE)

if PROGRAM_VERSION_RELEASE != 0:
    PROGRAM_VERSION_TEXT = "%d.%d.%d" % (PROGRAM_VERSION_MAJOR,
                                         PROGRAM_VERSION_MINOR,
                                         PROGRAM_VERSION_RELEASE)
else:
    PROGRAM_VERSION_TEXT = "%d.%d" % (PROGRAM_VERSION_MAJOR,
                                      PROGRAM_VERSION_MINOR)


#=============================================================================
# discretization


def get_time_point(i, ntimes, maxtime, delta=10):
    """Returns a discretized time point"""
    return (exp(i/float(ntimes) * log(1 + delta * maxtime)) - 1) / delta


def get_time_points(ntimes=30, maxtime=80000, delta=.01):
    """Returns a list of discretized time points"""
    return [get_time_point(i, ntimes, maxtime, delta)
            for i in range(ntimes+1)]


def iter_coal_states(tree, times):
    """Iterates through the coalescent states of a local tree"""
    
    # NOTE: do not use top time
    ntimes = len(times) - 1
    seen = set()
    time_lookup = dict((t, i) for i, t in enumerate(times))
    
    for node in tree.preorder():
        if len(node.children) == 1:
            continue
        i = time_lookup[node.age]
        
        if node.parents:
            parent = node.parents[0]
            while parent and parent not in seen:
                parent = parent.parents[0]
            
            while i < ntimes and times[i] <= parent.age:
                yield (node.name, i)
                i += 1
        else:
            while i < ntimes:
                yield (node.name, i)
                i += 1

        seen.add(node)


def get_nlineages_recomb_coal(tree, times):
    """
    Count the number of lineages at each time point that can coal and recomb
    """

    # TODO: is nrecombs including basal point?  It shouldn't
    
    nbranches = [0 for i in times]
    nrecombs = [0 for i in times]
    ncoals = [0 for i in times]

    for name, timei in iter_coal_states(tree, times):
        node = tree[name]

        # find parent node
        if node.parents:
            parent = node.parents[0]
            while len(parent.children) == 1:
                parent = parent.parents[0]
        else:
            parent = None

        # count who passes through this time segment
        if not parent or times[timei] < parent.age:
            nbranches[timei] += 1

        # count as recomb and coal point
        nrecombs[timei] += 1
        ncoals[timei] += 1
    nbranches[-1] = 1
    
    return nbranches, nrecombs, ncoals


def discretize_arg(arg, times, ignore_top=True, round_age="down"):
    """
    Round node ages to the nearest time point

    If 'ignore_top' is True, then do not use last time point for rounding
    """

    if ignore_top:
        times = times[:-1]
    
    for node in arg:
        i, j = util.binsearch(times, node.age)
        if j is None: j = len(times) - 1
        if i is None: i = 0
        
        if round_age == "up":
            node.age = times[j]
        elif round_age == "down":
            node.age = times[i]
        elif round_age == "closer":
            if node.age - times[i] < times[j] - node.age:
                node.age = times[i]
            else:
                node.age = times[j]
        else:
            raise Exception("unknown round_age '%s'" % round_age)


    recombs = [node for node in arg if node.event == "recomb"]
    recombs.sort(key=lambda x: x.pos)

    last = 0
    for node in recombs:
        intpos = int(node.pos)
        if intpos > last:
            node.pos = intpos
        else:
            node.pos = last + 1
        last = node.pos

    # ensure no duplicate recombinations
    seen = set()
    for node in arg:
        if node.event == "recomb":
            assert node.pos not in seen, (node.pos, sorted(seen))
            seen.add(node.pos)
            

def discretize_arg_recomb(arg):
    """Round recomb node to the nearest integer"""
    
    recombs = [node for node in arg if node.event == "recomb"]
    recombs.sort(key=lambda x: x.pos)

    last = 0
    for node in recombs:
        intpos = int(node.pos)
        if intpos > last:
            node.pos = intpos
        else:
            node.pos = last + 1
        last = node.pos

    # ensure no duplicate recombinations
    seen = set()
    for node in arg:
        if node.event == "recomb":
            assert node.pos not in seen, (node.pos, sorted(seen))
            seen.add(node.pos)
            

def get_treelen(tree, times, use_basal=True):
    """Calculate tree length"""
    treelen = sum(x.get_dist() for x in tree)
    if use_basal:
        rooti = times.index(tree.root.age)
        root_time = times[rooti+1] - times[rooti]
        treelen += root_time
    return treelen


def get_treelen_branch(tree, times, node, time, use_basal=True):
    """Calculate tree length with an extra branch"""

    treelen = sum(x.get_dist() for x in tree)

    blen = time
    treelen2 = treelen + blen
    if node == tree.root.name:
        treelen2 += blen - tree.root.age
        rooti = times.index(time)
        root_time = times[rooti+1] - times[rooti]
    else:
        rooti = times.index(tree.root.age)
        root_time = times[rooti+1] - times[rooti]

    if use_basal:
        treelen2 += root_time
        
    return treelen2


def get_basal_length(tree, times, node=None, time=None):
    """
    Get basal branch length
    
    NOTE: 'node' can be None
    """
    
    if node == tree.root.name:
        rooti = times.index(time)
        root_time = times[rooti+1] - times[rooti]
    else:
        rooti = times.index(tree.root.age)
        root_time = times[rooti+1] - times[rooti]
    
    return root_time



def sample_dsmc_sprs(k, popsize, rho, start=0.0, end=0.0, times=None, 
                     init_tree=None,
                     names=None, make_names=True):
    """
    Sample ARG using Discrete Sequentially Markovian Coalescent (SMC)

    k   -- chromosomes
    popsize  -- effective population size (haploid)
    rho -- recombination rate (recombinations / site / generation)
    start -- staring chromosome coordinate
    end   -- ending chromsome coordinate
    t   -- initial time (default: 0)
    names -- names to use for leaves (default: None)
    make_names -- make names using strings (default: True)
    """

    assert times is not None
    ntimes = len(times) - 1
    time_steps = [times[i] -  times[i-1] for i in range(1, ntimes+1)]
    if hasattr(popsize, "__len__"):
        popsizes = popsize
    else:
        popsizes = [popsize] * len(time_steps)


    # yield initial tree first
    if init_tree is None:
        init_tree = arglib.sample_arg(k, popsizes[0],
                                      rho=0.0, start=start, end=end,
                                      names=names, make_names=make_names)
        discretize_arg(init_tree, times, ignore_top=True)
    yield init_tree


    # sample SPRs
    pos = start
    tree = init_tree.copy()
    while True:
        # sample next recomb point
        treelen = sum(x.get_dist() for x in tree)
        blocklen = 0
        while blocklen == 0:
            blocklen = int(random.expovariate(max(treelen * rho, rho)))
        pos += blocklen
        if pos > end:
            break

        root_age_index = times.index(tree.root.age)

        # choose time interval for recombination
        states = set(iter_coal_states(tree, times))
        nbranches, nrecombs, ncoals = get_nlineages_recomb_coal(tree, times)
        probs = [nbranches[k] * time_steps[k]
                 for k in range(root_age_index+1)]
        #assert sum(probs) == get_treelen(tree, times)
        recomb_time_index = stats.sample(probs)
        recomb_time = times[recomb_time_index]

        # choose branch for recombination
        branches = [x for x in states if x[1] == recomb_time_index and
                    x[0] != tree.root.name]
        recomb_node = tree[random.sample(branches, 1)[0][0]]

        #treelib.draw_tree_names(tree.get_tree(), minlen=5, maxlen=5)

        # choose coal time
        j = recomb_time_index
        while j < ntimes - 1:
            kj = nbranches[j]
            if ((recomb_node.name, j) in states and
                recomb_node.parents[0].age > times[j]):
                kj -= 1
            assert kj > 0, (j, root_age_index, states)
            coal_prob = 1.0 - exp(- time_steps[j] * kj / float(popsizes[j]))
            if random.random() < coal_prob:
                break
            j += 1
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

        #print
        #treelib.draw_tree(tree.get_tree(), maxlen=5)
        #print rleaves, recomb_time
        #print cleaves, coal_time
        
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

        #print root_age_index == times.index(tree.root.age)


def sample_arg_dsmc(k, popsize, rho, start=0.0, end=0.0, times=None, 
                    init_tree=None,
                    names=None, make_names=True):
    """
    Returns an ARG sampled from the Discrete Sequentially Markovian Coalescent (SMC)
    
    k   -- chromosomes
    popsize -- effective population size
    rho -- recombination rate (recombinations / site / generation)
    start -- staring chromosome coordinate
    end   -- ending chromsome coordinate
    
    names -- names to use for leaves (default: None)
    make_names -- make names using strings (default: True)
    """

    if times is None:
        maxtime = 160000
        delta = .01
        ntimes = 20
        times = get_time_points(ntimes, maxtime, delta)
    
    it = sample_dsmc_sprs(
        k, popsize, rho, start=start, end=end, times=times, 
        init_tree=init_tree, names=names, make_names=make_names)
    tree = it.next()
    arg = arglib.make_arg_from_sprs(tree, it)
    #discretize_arg(arg, times, ignore_top=True)
    
    return arg


def sample_arg_mutations(arg, mu, times):
    """
    mu -- mutation rate (mutations/site/gen)
    """

    mutations = []
    minlen = times[1]

    for (start, end), tree in arglib.iter_tree_tracks(arg):
        arglib.remove_single_lineages(tree)
        for node in tree:
            if not node.parents:
                continue
            blen = max(node.get_dist(), minlen)
            rate = blen * mu
            i = start
            while i < end:
                i += int(min(random.expovariate(rate), 2*end))
                if i < end:
                    t = random.uniform(node.age, node.age + blen)
                    mutations.append((node, node.parents[0], i, t))
    return mutations



#=============================================================================
# recombination


def find_tree_next_recomb(arg, pos, tree=False):
    """Returns the next recombination node in a local tree"""

    recomb = None
    nextpos = util.INF

    if tree:
        nodes = iter(arg)
    else:
        nodes = arg.postorder_marginal_tree(pos-.5)

    for node in nodes:
        if node.event == "recomb" and node.pos >= pos and node.pos < nextpos:
            recomb = node
            nextpos = node.pos

    return recomb


def iter_visible_recombs(arg, start=None, end=None):
    """Iterates through visible recombinations in an ARG"""
    
    pos = start if start is not None else 0
    while True:
        recomb = find_tree_next_recomb(arg, pos+1)
        if recomb:
            yield recomb
            pos = recomb.pos
        else:
            break

        

#=============================================================================
# chromosome threads


def iter_chrom_thread(arg, node, by_block=True, use_clades=False):

    start = 0
    recombs = chain((x.pos for x in iter_visible_recombs(arg)),
                     [arg.end-1])

    for recomb_pos in recombs:
        #print recomb_pos
        if start >= arg.end:
            continue
        tree = arg.get_marginal_tree(recomb_pos-.5)
        block = [start, recomb_pos+1]
        start = recomb_pos+1
        
        # find parent
        node2 = tree[node.name]
        last = node2
        parent = node2.parents[0]
        while len(parent.children) == 1:
            last = parent
            parent = parent.parents[0]

        # find sibling
        c = parent.children
        sib = c[1] if last == c[0] else c[0]
        while len(sib.children) == 1:
            sib = sib.children[0]

        if use_clades:
            branch = list(tree.leaf_names(sib))
        else:
            branch = sib.name

        if by_block:
            yield (branch, parent.age, block)
        else:
            for i in range(block[0], block[1]):
                yield (branch, parent.age)


def get_coal_point(arg, node, pos):

    tree = arg.get_marginal_tree(pos-.5)

    # find parent
    node2 = tree[node.name]
    last = node2
    parent = node2.parents[0]
    while len(parent.children) == 1:
        last = parent
        parent = parent.parents[0]

    # find sibling
    c = parent.children
    sib = c[1] if last == c[0] else c[0]
    while len(sib.children) == 1:
        sib = sib.children[0]

    return sib.name, parent.age



def iter_chrom_timeline(arg, node, by_block=True):

    for node, time, block in iter_chrom_thread(arg, node, by_block=True):
        if by_block:
            yield (block[0]+1, time)
            yield (block[1], time)
        else:
            for i in range(block[0]+1, block[1]+1):
                yield time
            

        
def iter_posterior_times(model, probs, perc=.5):

    times = model.times

    for pos, probcol in enumerate(probs):
        col = [0.0] * len(times)

        for j, p in enumerate(probcol):
            node, timei = model.states[pos][j]
            col[timei] += exp(p)

        tot = 0.0
        j = 0
        while j < len(times) and tot < perc:
            tot += col[j]
            j += 1
        yield times[j-1]


def iter_thread_from_path(model, path):
    times = model.times
    states = model.states

    for pos, state in enumerate(path):
        node, timei = states[pos][state]
        yield node, times[timei]


def get_clade_point(arg, node_name, time, pos):
    """Returns a point along a branch in the ARG in terms of a clade and time"""

    if node_name in arg:
        tree = arg.get_marginal_tree(pos - .5)
        if (time > tree.root.age or
            (time == tree.root.age and node_name not in tree)):
            return (list(tree.leaf_names()), time)
        return (list(tree.leaf_names(tree[node_name])), time)
    else:
        return ([node_name], time)



#=============================================================================
# ARG operations


def make_trunk_arg(start, end, name="ind1"):
    """
    Returns a trunk genealogy
    """
    
    arg = arglib.ARG(start=start, end=end)
    node = arg.new_node(name, event="gene", age=0)
    return arg


def remove_arg_thread(arg, *chroms):
    """
    Remove a thread(s) from an ARG
    """
    remove_chroms = set(chroms)
    keep = [x for x in arg.leaf_names() if x not in remove_chroms]
    arglib.subarg_by_leaf_names(arg, keep)
    return arglib.smcify_arg(arg)

    '''
    def prune(sprs):
        for recomb_pos, (rleaves, rtime), (cleaves, ctime) in sprs:
            rleaves = [x for x in rleaves if x not in remove_chroms]
            cleaves = [x for x in cleaves if x not in remove_chroms]
            spr = (recomb_pos, (rleaves, rtime), (cleaves, ctime))
            print spr
            yield spr
    
    arg2 = arg.get_marginal_tree(-.5)
    keep = [x for x in arg2.leaf_names() if x not in remove_chroms]
    arglib.subarg_by_leaf_names(arg2, keep)
    arg2.write()
    
    #arglib.remove_single_lineages(arg2)
    sprs = prune(arglib.iter_arg_sprs(arg, use_leaves=True))
    arglib.make_arg_from_sprs(arg2, sprs, ignore_self=True)
    '''

    return arg2


def add_arg_thread(arg, new_name, thread, recombs):
    """Add a thread to an ARG"""

    def is_local_coal(arg, node, pos, local):
        return (len(node.children) == 2 and
                node.children[0] in local and
                arg.get_local_parent(node.children[0], pos-.5) == node and
                node.children[1] in local and
                arg.get_local_parent(node.children[1], pos-.5) == node and
                node.children[0] != node.children[1])



    def walk_up(arg, leaves, time, pos, ignore=None):

        order = dict((node, i) for i, node in enumerate(
            arg.postorder_marginal_tree(pos-.5)))
        local = set(order.keys())
        if ignore is not None and ignore in arg:
            ptr = arg[ignore]
            if ptr in local:
                local.remove(ptr)
                ptr = arg.get_local_parent(ptr, pos-.5)
            else:
                ptr = None
            
            while ptr and ptr in local:
                if (len(ptr.children) == 2 and
                    ((ptr.children[0] in local and
                      arg.get_local_parent(ptr.children[0], pos-.5) == ptr) or
                     (ptr.children[1] in local and
                      arg.get_local_parent(ptr.children[1], pos-.5) == ptr))):
                    break
                local.remove(ptr)
                ptr = arg.get_local_parent(ptr, pos-.5)

        queue = [(order[arg[x]], arg[x]) for x in leaves]
        seen = set(x[1] for x in queue)
        heapq.heapify(queue)

        while len(queue) > 1:
            i, node = heapq.heappop(queue)
            parent = arg.get_local_parent(node, pos-.5)
            if parent and parent not in seen:
                seen.add(parent)
                heapq.heappush(queue, (order[parent], parent))
        node = queue[0][1]
        parent = arg.get_local_parent(node, pos-.5)

        
        while parent and parent.age <= time:
            if is_local_coal(arg, parent, pos, local):
                break
            node = parent
            parent = arg.get_local_parent(node, pos-.5)

        if parent:
            if parent.age < time:
                print leaves, parent.age, time, ignore
                tree = arg.get_marginal_tree(pos-.5).get_tree()
                tree.write()
                treelib.draw_tree_names(tree, maxlen=8, minlen=8)
                assert False

        return node


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


    arg_recomb = dict((x.pos, x) for x in iter_visible_recombs(arg))
    recomb_clades = [
        (pos-1, None) + get_clade_point(arg, rnode, rtime, pos-1)
        for pos, rnode, rtime in recombs] + [
        (node.pos, node.name) +
        get_clade_point(arg, node.name, node.age, node.pos)
        for node in iter_visible_recombs(arg)]
    recomb_clades.sort()

    # make initial tree
    arg2 = arg.get_marginal_tree(-1)
    arglib.remove_single_lineages(arg2)

    start = get_clade_point(arg, thread[0][0], thread[0][1], 0)
    node = walk_up(arg2, start[0], start[1], -1)
    node2 = add_node(arg2, node, start[1], -1, "coal")
    leaf = arg2.new_node(name=new_name, event="gene", age=0)
    leaf.parents.append(node2)
    node2.children.append(leaf)
    

    # add each recomb and re-coal
    for rpos, rname, rleaves, rtime in recomb_clades:
        if rpos in arg_recomb:
            # find re-coal for existing recomb

            if thread[rpos][1] != thread[rpos+1][1]:
                if rtime > min(thread[rpos][1], thread[rpos+1][1]):
                    print ">>", rtime, thread[rpos], thread[rpos+1]
                    treelib.draw_tree_names(
                        arg.get_marginal_tree(rpos-.5).get_tree(),
                        maxlen=8, minlen=8)
                    treelib.draw_tree_names(
                        arg.get_marginal_tree(rpos+.5).get_tree(),
                    maxlen=8, minlen=8)
                    assert False
            
            node = arg_recomb[rpos]
            local1 = set(arg.postorder_marginal_tree(rpos-.5))
            local2 = set(arg.postorder_marginal_tree(rpos+.5))
            last = node
            node = arg.get_local_parent(node, rpos+.5)
            while (not is_local_coal(arg, node, rpos+1, local2)):
                last = node
                node = arg.get_local_parent(node, rpos+.5)
            c = node.children
            child = c[0] if c[1] == last else c[1]
            recoal = node
            
            cleaves, ctime = get_clade_point(
                arg, child.name, node.age, rpos-.5)

            # get local tree T^{n-1}_i and add new branch
            tree = arg.get_marginal_tree(rpos+.5)
            arglib.remove_single_lineages(tree)            
            node_name, time = thread[rpos+1]
            node = tree[node_name]

            # add new branch
            node2 = add_node(tree, node, time, rpos+1, "coal")
            if not node2.parents:
                tree.root = node2
            leaf = tree.new_node(name=new_name, event="gene", age=0)
            leaf.parents.append(node2)
            node2.children.append(leaf)
            
            recomb = walk_up(tree, rleaves, rtime, rpos+1, new_name)

            if recomb == node2 and rtime == node2.age:
                # recomb and new coal-state are near each other
                # we must decide if recomb goes above or below coal-state

                # if this is a mediated SPR, then recomb goes below.
                # otherwise it goes above.

                # SPR is mediated if previous coal state is not recomb branch
                node_name, time = thread[rpos]
                if node2.children[0].name != node_name:
                    # this is a mediated coal
                    recomb = node2.children[0]
            
            coal = recomb.parents[0]
            c = coal.children
            child = c[0] if c[1] == recomb else c[1]

            # get coal point in T^n_i
            rleaves, rtime = get_clade_point(
                tree, recomb.name, rtime, rpos+1)
            cleaves, ctime = get_clade_point(
                tree, child.name, coal.age, rpos+1)

            node1 = walk_up(arg2, rleaves, rtime, rpos+1)
            node2 = walk_up(arg2, cleaves, ctime, rpos+1, node1.name)

    
        else:
            # find re-coal for new recomb
            
            assert rtime <= thread[rpos][1], (rtime, thread[rpos][1])
            
            if rleaves == [new_name]:
                # recomb on new branch, coal given thread
                cleaves, ctime = get_clade_point(
                    arg, thread[rpos+1][0], thread[rpos+1][1], rpos+.5)
                assert ctime >= rtime, (rtime, ctime)
                
                node1 = walk_up(arg2, rleaves, rtime, rpos+1)
                node2 = walk_up(arg2, cleaves, ctime, rpos+1, new_name)
                
            else:
                # recomb in ARG, coal on new branch
                cleaves = [new_name]
                ctime = thread[rpos+1][1]
                assert ctime >= rtime, (rtime, ctime)

                # NOTE: new_name is not ignored for walk_up on rleaves
                # because I do not want the recombination to be higher
                # than the coal point, which could happen if the recomb time
                # is the same as the current coal time.                
                node1 = walk_up(arg2, rleaves, rtime, rpos+1)
                node2 = walk_up(arg2, cleaves, ctime, rpos+1, node1.name)


        assert node1.parents
        assert rtime <= ctime

        recomb = add_node(arg2, node1, rtime, rpos, "recomb")
        if node1 == node2:
            node2 = recomb
        coal = add_node(arg2, node2, ctime, rpos, "coal")

        recomb.parents.append(coal)
        coal.children.append(recomb)

        node, time = get_coal_point(arg2, arg2[new_name], rpos+1)
        assert time == thread[rpos+1][1], (time, thread[rpos+1][1])

    
    
    return arg2
    


def arg_lca(arg, leaves, time, pos, ignore=None):
    """Returns Least Common Ancestor for leaves in an ARG at position 'pos'"""

    def is_local_coal(arg, node, pos, local):
        return (len(node.children) == 2 and
                node.children[0] in local and
                arg.get_local_parent(node.children[0], pos-.5) == node and
                node.children[1] in local and
                arg.get_local_parent(node.children[1], pos-.5) == node and
                node.children[0] != node.children[1])


    order = dict((node, i) for i, node in enumerate(
        arg.postorder_marginal_tree(pos-.5)))
    local = set(order.keys())
    if ignore is not None and ignore in arg:
        ptr = arg[ignore]
        local.remove(ptr)
        ptr = arg.get_local_parent(ptr, pos-.5)

        while ptr and ptr in local:
            if (len(ptr.children) == 2 and
                ((ptr.children[0] in local and
                  arg.get_local_parent(ptr.children[0], pos-.5) == ptr) or
                 (ptr.children[1] in local and
                  arg.get_local_parent(ptr.children[1], pos-.5) == ptr))):
                break
            local.remove(ptr)
            ptr = arg.get_local_parent(ptr, pos-.5)

    queue = [(order[arg[x]], arg[x]) for x in leaves]
    seen = set(x[1] for x in queue)
    heapq.heapify(queue)

    while len(queue) > 1:
        i, node = heapq.heappop(queue)
        parent = arg.get_local_parent(node, pos-.5)
        if parent and parent not in seen:
            seen.add(parent)
            heapq.heappush(queue, (order[parent], parent))
    node = queue[0][1]
    parent = arg.get_local_parent(node, pos-.5)

    # walk up appropriate time if given
    if time is not None:
        while parent and parent.age <= time:
            if is_local_coal(arg, parent, pos, local):
                break            
            node = parent
            parent = arg.get_local_parent(node, pos-.5)

        if parent:
            if parent.age < time:
                print (leaves, parent.age, time)
                tree = arg.get_marginal_tree(pos-.5).get_tree()
                tree.write()
                treelib.draw_tree_names(tree, maxlen=8, minlen=8)
                assert False

    return node


def find_recomb_coal(tree, last_tree, recomb_name=None, pos=None):
    """
    Returns the recomb and coal points for the SPR between two trees
    """

    if recomb_name is None:
        recomb = find_tree_next_recomb(last_tree, pos-1, tree=True)
        recomb_name = recomb.name
    
    # find recomb node
    recomb_node = tree[recomb_name]
    recomb_time = recomb_node.age

    # find re-coal point
    coal = recomb_node.parents[0]
    while coal.name not in last_tree and coal.parents:
        coal = coal.parents[0]
    coal_time = coal.age

    # find coal branch in last_tree
    if coal.name not in last_tree:
        # coal above root
        coal_branch = last_tree.root.name
    else:
        ptr = last_tree[coal.name]
        while len(ptr.children) == 1:
            ptr = ptr.children[0]
        coal_branch = ptr.name

    # find recomb branch in tree
    recomb = tree[recomb_name]
    while len(recomb.children) == 1:
        recomb = recomb.children[0]
    recomb_branch = recomb.name

    return (recomb_branch, recomb_time), (coal_branch, coal_time)





#=============================================================================
# probabilities


def calc_C(time_steps, nbranches, popsizes):
    ntimes = len(time_steps)
    C = [0.0]
    for k in xrange(1, ntimes):
        l = k - 1
        C.append(C[-1] + time_steps[l] * nbranches[l] / (2.0 * popsizes[l]))
    return C


def calc_A_matrix(time_steps, nbranches, popsizes):

    ntimes = len(time_steps)
    
    # A_{k,j} =& s'_{j-2} k_{j-2} / (2N) + \sum_{m=k}^{j-3} s'_m k_m / (2N) \\
    #         =& s'_{j-2} k_{j-2} / (2N) + A_{k,j-1}.
    
    A = util.make_matrix(ntimes, ntimes, 0.0)
    for k in xrange(ntimes):
        # A[k][k] = 0
        for j in xrange(k+1, ntimes):
            l = j - 1
            A[k][j] = A[k][j-1] + time_steps[l] * nbranches[l] / (2.0 * popsizes[l])
    return A




def calc_transition_probs(tree, states, nlineages, times,
                          time_steps, popsizes, rho):

    nstates = len(states)
    ntimes = len(time_steps)
    minlen = time_steps[0]
    treelen = sum(x.get_dist() for x in tree)
    nbranches, nrecombs, ncoals = nlineages

    # calculate base case (time=0)
    root_age_index = times.index(tree.root.age)
    treelen_b = treelen + time_steps[root_age_index];
    C = [0.0]
    B = [(nbranches[0] + 1) * time_steps[0] / (nrecombs[0] + 1.0)]
    D = [(1.0 - exp(-max(rho * treelen, rho))) / treelen_b]
    E = [(1.0 - exp(-time_steps[0] * nbranches[0]
                    / (2.0 * popsizes[0]))) / ncoals[0]]
    G = [time_steps[0] / (nrecombs[0] + 1.0)]
    norecombs = [exp(-max(rho * treelen, rho))]

    # calculate all other time points (time>0)
    for b in range(1, ntimes-1):
        # get tree length
        treelen2 = treelen + times[b]
        if b > root_age_index:
            # add wrapped branch
            treelen2 += times[b] - tree.root.age

            # add basal branch
            treelen2_b = treelen2 + time_steps[b]
        else:
            # add basal branch
            treelen2_b = treelen2 + time_steps[root_age_index]

        # due to normalization we do not need exp(-rho * treelen)
        l = b - 1;
        C.append(C[l] + time_steps[l] * nbranches[l] / (2.0 * popsizes[l]))
        eC = exp(C[b])

        B.append(B[b-1] + (nbranches[b] + 1.0) * time_steps[b] / 
                 (nrecombs[b] + 1.0) * eC)
        D.append((1.0 - exp(-max(rho * treelen2, rho))) / treelen2_b)
        E.append((1.0 - exp(-time_steps[b] * nbranches[b] / 
                            (2.0 * popsizes[b]))) / eC / ncoals[b])
        G.append(eC * time_steps[b] / (nrecombs[b] + 1.0))
        norecombs.append(exp(-max(rho * treelen2, rho)))
    E[ntimes-2] = exp(-C[ntimes-2]) / ncoals[ntimes-2]
    

    # calculate full state transition matrix
    transprob = util.make_matrix(nstates, nstates, 0.0)
    for i in range(nstates):
        node1, a = states[i]
        c = times.index(tree[node1].age)
        
        for j in range(nstates):
            node2, b = states[j]
            I = float(a <= b)
            
            if node1 != node2:
                transprob[i][j] = D[a] * E[b] * (B[min(a,b)] - I * G[a])
            else:
                #print "t", a, b, D[a], E[b], B[min(a,b)], norecombs[a], time_steps[:a]
                Bc = B[c-1] if c > 0 else 0.0
                transprob[i][j] = D[a] * E[b] * \
                    (2 * B[min(a,b)] - 2 * I * G[a] - Bc)
                if a == b:
                    transprob[i][j] += norecombs[a]

        # normalize and convert to log scale
        s = sum(transprob[i])
        for j in range(nstates):
            transprob[i][j] = log(transprob[i][j] / s)

    return transprob


def calc_no_recomb_cond_self(tree, states, nlineages, times,
                             time_steps, popsizes, rho):

    nstates = len(states)
    ntimes = len(time_steps)
    minlen = time_steps[0]
    treelen = sum(x.get_dist() for x in tree)
    nbranches, nrecombs, ncoals = nlineages

    # calculate base case (time=0)
    root_age_index = times.index(tree.root.age)
    treelen_b = treelen + time_steps[root_age_index];
    C = [0.0]
    B = [(nbranches[0] + 1) * time_steps[0] / (nrecombs[0] + 1.0)]
    D = [(1.0 - exp(-max(rho * treelen, rho))) / treelen_b]
    E = [(1.0 - exp(-time_steps[0] * nbranches[0]
                    / (2.0 * popsizes[0]))) / ncoals[0]]
    G = [time_steps[0] / (nrecombs[0] + 1.0)]
    norecombs = [exp(-max(rho * treelen, rho))]

    # calculate all other time points (time>0)
    for b in range(1, ntimes-1):
        # get tree length
        treelen2 = treelen + times[b]
        if b > root_age_index:
            # add wrapped branch
            treelen2 += times[b] - tree.root.age

            # add basal branch
            treelen2_b = treelen2 + time_steps[b]
        else:
            # add basal branch
            treelen2_b = treelen2 + time_steps[root_age_index]

        # due to normalization we do not need exp(-rho * treelen)
        l = b - 1;
        C.append(C[l] + time_steps[l] * nbranches[l] / (2.0 * popsizes[l]))
        eC = exp(C[b])

        B.append(B[b-1] + (nbranches[b] + 1.0) * time_steps[b] / 
                 (nrecombs[b] + 1.0) * eC)
        D.append((1.0 - exp(-max(rho * treelen2, rho))) / treelen2_b)
        E.append((1.0 - exp(-time_steps[b] * nbranches[b] / 
                            (2.0 * popsizes[b]))) / eC / ncoals[b])
        G.append(eC * time_steps[b] / (nrecombs[b] + 1.0))
        norecombs.append(exp(-max(rho * treelen2, rho)))
    E[ntimes-2] = exp(-C[ntimes-2]) / ncoals[ntimes-2]
    

    # calculate full state transition matrix
    transprob = util.make_matrix(nstates, nstates, 0.0)
    vec = []
    for i in range(nstates):
        node1, a = states[i]
        c = times.index(tree[node1].age)
        
        for j in range(nstates):
            node2, b = states[j]
            I = float(a <= b)
            
            if node1 != node2:
                transprob[i][j] = D[a] * E[b] * (B[min(a,b)] - I * G[a])
            else:
                #print "t", a, b, D[a], E[b], B[min(a,b)], norecombs[a], time_steps[:a]
                Bc = B[c-1] if c > 0 else 0.0
                transprob[i][j] = D[a] * E[b] * \
                    (2 * B[min(a,b)] - 2 * I * G[a] - Bc)
                if a == b:
                    transprob[i][j] += norecombs[a]

        # normalize and convert to log scale
        #s = sum(transprob[i])
        #for j in range(nstates):
        #    transprob[i][j] = log(transprob[i][j] / s)

        vec.append(norecombs[a] / transprob[i][i])

    return vec



def get_recomb_transition_switch(tree, last_tree, spr, states1, states2,
                                 times):
    
    # SPR subtree moves out from underneath us
    # therefore therefore the new chromosome coalesces with
    # the branch above the subtree

    (recomb_branch, recomb_time), (coal_branch, coal_time) = spr

    # search up for parent
    recomb = last_tree[recomb_branch]
    parent = recomb.parents[0]
    b = times.index(parent.age)

    # find other child
    c = parent.children
    other = (c[0] if c[1] == recomb else c[1])

    # find new state in tree
    if other.name == coal_branch:
        next_state = (tree[other.name].parents[0].name, b)
    else:
        next_state = (other.name, b)
    
    a = states2.index((recomb_branch, recomb_time))
    b = states2.index(next_state)
    return (a, b)

                


def calc_transition_probs_switch(tree, last_tree, recomb_name,
                                 states1, states2,
                                 nlineages, times,
                                 time_steps, popsizes, rho):

    treelen = get_treelen(last_tree, times)
    nbranches, nrecombs, ncoals = nlineages
    (recomb_branch, recomb_time), (coal_branch, coal_time) = \
        find_recomb_coal(tree, last_tree, recomb_name=recomb_name)

    k = times.index(recomb_time)
    coal_time = times.index(coal_time)
    
    last_tree2 = last_tree.copy()
    arglib.remove_single_lineages(last_tree2)
    tree2 = tree.copy()
    arglib.remove_single_lineages(tree2)
    
    # compute transition probability matrix
    transprob = util.make_matrix(len(states1), len(states2), -util.INF)
    
    determ = get_deterministic_transitions(states1, states2, times,
                                           tree2, last_tree2,
                                           recomb_branch, k,
                                           coal_branch, coal_time)


    for i, (node1, a) in enumerate(states1):        
        if (node1, a) == (recomb_branch, k):
            # probabilistic transition case (recomb case)
            spr = (recomb_branch, k), (coal_branch, coal_time)
            recomb_next_states = get_recomb_transition_switch(
                tree2, last_tree2, spr, states1, states2, times)

            # placeholders
            transprob[i][recomb_next_states[0]] = log(.5)
            transprob[i][recomb_next_states[1]] = log(.5)

    
        elif (node1, a) == (coal_branch, coal_time):
            # probabilistic transition case (re-coal case)

            # determine if node1 is still here or not
            last_recomb = last_tree2[recomb_branch]
            last_parent = last_recomb.parents[0]
            if last_parent.name == node1:
                # recomb breaks node1 branch, we need to use the other child
                c = last_parent.children
                node3 = c[0].name if c[1] == last_recomb else c[1].name
            else:
                node3 = node1
            

            # find parent of recomb_branch and node1
            last_parent_age = times.index(last_parent.age)
            parent = tree2[recomb_branch].parents[0]
            assert parent == tree2[node3].parents[0]

            # treelen of T^n_{i-1}
            blen = times[a]
            treelen2 = treelen + blen
            if node1 == last_tree2.root.name:
                treelen2 += blen - last_tree2.root.age
                treelen2 += time_steps[a]
            else:
                treelen2 += time_steps[times.index(last_tree2.root.age)]


            for j, (node2, b) in enumerate(states2):
                transprob[i][j] = 0.0
                if not ((node2 == recomb_branch and b >= k) or
                        (node2 == node3 and b == a) or
                        (node2 == parent.name and b == a)):
                    continue

                # get lineage counts
                # remove recombination branch and add new branch
                kbn = nbranches[b]
                kcn = ncoals[b] + 1
                if times[b] < parent.age:
                    kbn -= 1
                    kcn -= 1
                if b < a:
                    kbn += 1
                
                twon = 2.0 * popsizes[b]

                transprob[i][j] = (
                    (1.0 - exp(- time_steps[b] * kbn / twon)) / kcn *
                    exp(- sum(time_steps[m] * (nbranches[m] + 1
                              - (1 if m < last_parent_age else 0))
                              / (2.0 * popsizes[m])
                              for m in xrange(k, b))))

            # normalize row to ensure they add up to one
            tot = sum(transprob[i])
            for j in xrange(len(states2)):
                x = transprob[i][j]
                if tot > 0.0 and x > 0.0:
                    transprob[i][j] = log(x / tot)
                else:
                    transprob[i][j] = -1e1000

        else:
            # deterministic transition
            assert determ[i] != -1, determ
            transprob[i][determ[i]] = 0.0


    return transprob




def get_deterministic_transitions(states1, states2, times,
                                  tree, last_tree,
                                  recomb_branch, recomb_time,
                                  coal_branch, coal_time):

    # recomb_branch in tree and last_tree
    # coal_branch in last_tree    
    
    state2_lookup = util.list2lookup(states2)
    
    next_states = []
    for i, state1 in enumerate(states1):
        node1, a = state1
        
        if (node1, a) == (coal_branch, coal_time):
            # not a deterministic case
            next_states.append(-1)
        
        elif node1 != recomb_branch:
            # SPR only removes a subset of descendents, if any
            # trace up from remaining leaf to find correct new state

            node = last_tree.nodes.get(node1, None)
            if node is None:
                print node1
                treelib.draw_tree_names(last_tree.get_tree(),
                                        minlen=8, maxlen=8)
                raise Exception("unknown node name '%s'" % node1)

            
            if node.is_leaf():
                # SPR can't disrupt leaf branch
                node2 = node1

            else:
                child1 = node.children[0]
                child2 = node.children[1]
                
                if recomb_branch == child1.name:
                    # right child is not disrupted
                    node2 = child2.name

                elif recomb_branch == child2.name:
                    # left child is not disrupted
                    node2 = child1.name

                else:
                    # node is not disrupted
                    node2 = node1

            # optionally walk up
            if ((coal_branch == node1 or coal_branch == node2) and
                coal_time <= a):
                # coal occurs under us
                node2 = tree[node2].parents[0].name
            next_states.append(state2_lookup[(node2, a)])
                
        else:
            # SPR is on same branch as new chromosome
            if recomb_time >= a:
                # we move with SPR subtree
                # TODO: we could probabilistically have subtree move
                # out from underneath.
                next_states.append(state2_lookup[(recomb_branch, a)])

            else:
                # SPR should not be able to coal back onto same branch
                # this would be a self cycle
                assert coal_branch != node1
                
                # SPR subtree moves out from underneath us
                # therefore therefore the new chromosome coalesces with
                # the branch above the subtree

                # search up for parent
                recomb = last_tree[recomb_branch]
                parent = recomb.parents[0]
                b = times.index(parent.age)

                # find other child
                c = parent.children
                other = (c[0] if c[1] == recomb else c[1])

                # find new state in tree
                if other.name == coal_branch:
                    next_state = (tree[other.name].parents[0].name, b)
                else:
                    next_state = (other.name, b)
                
                next_states.append(state2_lookup[next_state])

    return next_states


def calc_state_priors(tree, states, nlineages,
                      times, time_steps, popsizes, rho):
    """Calculate state priors"""

    priormat = [
        log((1 - exp(- time_steps[b] * nlineages[0][b] /
                 (2.0 * popsizes[b]))) / nlineages[2][b] *
             exp(-sum(time_steps[m] * nlineages[0][m] /
                      (2.0 * popsizes[m])
                      for m in range(0, b))))
            for node, b in states]
    
    return priormat







#=============================================================================
# ArgHmm model

"""
        bases     0     1     2     3     4     5
               |-----|-----|-----|-----|-----|-----|
recomb points        0     1     2     3     4     5
                                 *           *
local blocks   (0, 3) (3, 5) (5, 6)

"""


class ArgHmm (hmm.HMM):

    def __init__(self, arg, seqs, new_name=None,
                 popsize=1e4, rho=1.5e-8, mu=2.5e-8,
                 times=None,
                 ntimes=30, maxtime=100000.0, delta=.01):

        assert arg.start == 0

        # setup model
        self.new_name = new_name
        if times is None:
            self.times = get_time_points(ntimes, maxtime, delta)
        else:
            self.times = times
            ntimes = len(self.times) - 1
        self.time_steps = [self.times[i] -  self.times[i-1]
                           for i in range(1, ntimes+1)]
        self.time_steps.append(maxtime*10000.0)
        self.ntimes = ntimes
        self.arg = arg
        self.seqs = seqs
        self.rho = rho
        self.mu = mu
        if hasattr(popsize, "__len__"):
            self.popsizes = popsize
        else:
            self.popsizes = [popsize] * len(self.time_steps)

        # assert none of the recombination are between the same sites
        recombs = [int(x.pos) for x in arg if x.event == "recomb"]
        assert len(recombs) == len(set(recombs))

        # determine states
        self.recomb_pos = [-1] + list(
            x.pos for x in iter_visible_recombs(arg))
        self.recomb_pos.append(arg.end - 1)
        self.states = []
        self.state_spaces = []
        j = 0
        last_recomb = None
        for i in xrange(arg.end):
            while j < len(self.recomb_pos) and i > self.recomb_pos[j]:
                j += 1
                last_recomb = None
            if j != last_recomb:
                last_recomb = j
                self.states.append(
                    list(iter_coal_states(arg.get_marginal_tree(i-.5),
                                          self.times)))
                self.state_spaces.append(j-1)
            else:
                self.states.append(self.states[-1])
                self.state_spaces.append(j-1)

        #print self.state_spaces[:30]

        # current local tree
        self.local_block = [-1, self.recomb_pos[1]]
        self.local_tree = None
        self.emit_col = None
        self.last_tree = None
        self.last_pos = None
        self.transmat = None
        self.transmat_switch = None
        
        self.check_local_tree(0, force=True)


    def get_state_space(self, pos):
        """Returns the state_space ids for (pos-1, pos)"""
        return self.state_spaces[pos]


    def get_local_block(self, space):
        return (self.recomb_pos[space]+1, self.recomb_pos[space+1]+1)
        


    def check_local_tree(self, pos, force=False):

        # update local block
        if force or not (self.local_block[0] <=  pos < self.local_block[1]):

            # get new local information
            self.local_tree = self.arg.get_marginal_tree(pos-.5)
            self.local_block = self.get_local_block(self.state_spaces[pos])
            self.nlineages = get_nlineages_recomb_coal(
                self.local_tree, self.times)
            
            # get new transition matrices
            self.transmat = calc_transition_probs(
                self.local_tree, self.states[pos], self.nlineages,
                self.times, self.time_steps, self.popsizes, self.rho)

            assert len(self.transmat) == len(self.states[pos])
            assert len(self.transmat[0]) == len(self.states[pos])

            # get switch matrix for beginning of block
            start = self.local_block[0]
            recomb = find_tree_next_recomb(self.arg, start - 1)
            if start > 0 and recomb is not None:
                last_tree = self.arg.get_marginal_tree(start-1-.5)
                self.transmat_switch = calc_transition_probs_switch(
                    self.local_tree, last_tree, recomb.name,
                    self.states[start-1], self.states[start],
                    self.nlineages, self.times,
                    self.time_steps, self.popsizes, self.rho)

                assert len(self.transmat_switch) == len(self.states[start-1])
                assert len(self.transmat_switch[0]) == len(self.states[start])
            else:
                self.transmat_switch = None

            # get prior matrix if needed
            self.priormat = calc_state_priors(
                self.local_tree, self.states[pos], self.nlineages,
                self.times, self.time_steps, self.popsizes, self.rho)

            # makes computing emissions easier
            arglib.remove_single_lineages(self.local_tree)


        # update local site
        if force or pos != self.last_pos:
            self.emit_col = emit.calc_emission(self.local_tree, self, pos,
                                               self.new_name)
            

        self.last_pos = pos


    def get_num_states(self, pos):
        return len(self.states[pos])


    def prob_prior(self, pos, state):

        self.check_local_tree(pos)
        return self.priormat[state]
    
        
    def prob_transition(self, pos1, state1, pos2, state2):

        assert pos1 == pos2 - 1
        self.check_local_tree(pos2)
        
        if pos2 == self.local_block[0] and self.transmat_switch:
            return self.transmat_switch[state1][state2]
        else:
            return self.transmat[state1][state2]
        

    def prob_emission(self, pos, state):

        self.check_local_tree(pos)
        return self.emit_col[state]
    


    
#=============================================================================
# HMM methods



def iter_trans_emit_matrices(model, n):

    last_tree = None
    last_nlineages = None
    
    # get transition matrices and emissions
    for rpos in model.recomb_pos[:-1]:
        pos = rpos + 1
        
        # get new local information
        tree = model.arg.get_marginal_tree(pos-.5)
        block = model.get_local_block(model.state_spaces[pos])
        nlineages = get_nlineages_recomb_coal(tree, model.times)
        nbranches, nrecombs, ncoals = nlineages
        times_lookup = dict((t, i) for i, t in enumerate(model.times))
        tree2 = tree.get_tree()
        ptree, nodes, nodelookup = make_ptree(tree2)
        int_states = [[nodelookup[tree2[node]], timei]
                      for node, timei in model.states[pos]]
        nstates = len(int_states)
        ages = [tree[node.name].age for node in nodes]
        ages_index = [times_lookup[tree[node.name].age]
                      for node in nodes]
        #treelen = sum(x.dist for x in tree2)
        treelen = get_treelen(tree, model.times)

            
        # get new transition matrices
        transmat = new_transition_probs(
            len(nodes), ptree, ages_index, treelen, 
            ((c_int * 2) * nstates)
            (* ((c_int * 2)(n, t) for n, t in int_states)), nstates,
            len(model.time_steps), model.times, model.time_steps,
            nbranches, nrecombs, ncoals, 
            model.popsizes, model.rho)

        
        # get switch matrix for beginning of block
        start = block[0]
        recomb = find_tree_next_recomb(model.arg, start - 1)
        if start > 0 and recomb is not None:
            assert last_tree is not None

            # debug:
            states1 = list(iter_coal_states(last_tree, model.times))
            if states1 != model.states[start-1]:
                print "start", start
                print "recomb_pos", model.recomb_pos
                print "states1", states1
                print "model.states-1", model.states[start-1]
                print "model.states", model.states[start]
                
                assert states1 == model.states[start-1]
            
            #transmat_switch = calc_transition_probs_switch(
            #    tree, last_tree, recomb.name,
            #    model.states[start-1], model.states[start],
            #    last_nlineages, model.times,
            #    model.time_steps, model.popsizes, model.rho)

            transmat_switch = calc_transition_probs_switch_c(
                tree, last_tree, recomb.name,
                model.states[start-1], model.states[start],
                last_nlineages, model.times,
                model.time_steps, model.popsizes, model.rho, raw=False)
            
        else:
            transmat_switch = None

        # get emission matrix
        seqs = [model.seqs[node.name][block[0]:block[1]]
                for node in nodes if node.is_leaf()]
        seqs.append(model.seqs[model.new_name][block[0]:block[1]])

        emit = new_emissions(
            ((c_int * 2) * nstates)
            (* ((c_int * 2)(n, t) for n, t in int_states)), nstates, 
            ptree, len(ptree), ages_index,
            (c_char_p * len(seqs))(*seqs), len(seqs), len(seqs[0]),
            model.times, len(model.times), model.mu)

        last_tree = tree
        last_nlineages = nlineages

        yield block, nstates, transmat, transmat_switch, emit

    

def forward_algorithm(model, n, verbose=False, matrices=None):

    probs = []

    if verbose:
        util.tic("forward")


    (ptrees, ages, sprs, blocks), all_nodes = get_treeset(
        model.arg, model.times)
    blocklens = [x[1] - x[0] for x in blocks]
    seqlen = sum(blocklens)

    seqs = [model.seqs[node] for node in all_nodes[0]
            if model.arg[node].is_leaf()]
    seqs.append(model.seqs[model.new_name])
    nnodes = len(ptrees[0])
    fw = arghmm_forward_alg(ptrees, ages, sprs, blocklens,
                            len(ptrees), nnodes, 
                            model.times, len(model.times),
                            model.popsizes, model.rho, model.mu,
                            (c_char_p * len(seqs))(*seqs), len(seqs),
                            seqlen)

    # map states to python state space
    all_states = get_state_spaces(ptrees, ages, sprs, blocklens,
                                  len(ptrees), nnodes, len(model.times))

    probs = []
    for k, (start, end) in enumerate(blocks):
        states = model.states[start]
        nstates = len(states)
        istates = all_states[k]
        nodes = all_nodes[k]
        lookup = util.list2lookup(states)
        
        mapping = [0] * nstates
        for j, (inode, itime) in enumerate(istates[:nstates]):
            s = (nodes[inode], itime)
            mapping[lookup[s]] = j
        
        for i in range(start, end):
            col = []
            probs.append(col)
            for j in xrange(nstates):
                col.append(fw[i][mapping[j]])


    delete_state_spaces(all_states, len(ptrees))
    delete_double_matrix(fw, seqlen)

    if verbose:
        util.toc()
            
    return probs



def backward_algorithm(model, n, verbose=False, matrices=None):

    probs = []

    if verbose:
        util.tic("backward")

    # get prior matrix
    local_tree = model.arg.get_marginal_tree(n-.5)
    nlineages = get_nlineages_recomb_coal(local_tree, model.times)
    priors = calc_state_priors(
        local_tree, model.states[n-1], nlineages,
        model.times, model.time_steps, model.popsizes, model.rho)
    
    for i in xrange(n):
        probs.append(None)
    nstates = model.get_num_states(n-1)
    probs[n-1] = [priors[j] + model.prob_emission(n-1, j)
                  for j in xrange(nstates)]

    # get transition matrices
    matrices_given = matrices is not None
    if matrices is None:
        matrices = list(iter_trans_emit_matrices(model, n))

    # iterate over blocks
    for block, nstates, transmat, transmat_switch, emit in reversed(matrices):
        if verbose:
            util.logger(" pos %d" % block[0])

        blocklen = block[1] - block[0]

        # use transmat for rest of block
        # make forward table for block
        bw = []
        for pos in xrange(block[0], block[1]-1):
            bw.append([0.0 for k in xrange(nstates)])
        assert len(probs[block[1]-1]) == nstates
        bw.append(probs[block[1]-1])

        assert len(bw) == blocklen
        
        backward_alg(blocklen, nstates, transmat, emit, bw)

        i = block[0]
        for j in range(blocklen-1):
            probs[i+j] = bw[j][:nstates]

        # use switch matrix for first col
        if block[0] > 0:
            nstates1 = len(transmat_switch)
            nstates2 = len(transmat_switch[0])
            i = block[0] - 1
            e = emit[0]
            col2 = probs[i+1]
            assert len(col2) == nstates2 == nstates
            col1 = []
            for j in xrange(nstates1):
                col1.append(stats.logsum([
                    col2[k] + e[k] + transmat_switch[j][k]
                    for k in xrange(nstates2)]))
            probs[i] = col1

        if not matrices_given:
            delete_emissions(emit, blocklen)
            delete_transition_probs(transmat, nstates)


    if verbose:
        util.toc()
            
    return probs


def get_posterior_probs(model, n, verbose=False,
                        probs_forward=None, probs_backward=None):
    """
    Returns posterior decoding of thread state
    """
    
    matrices = None

    if probs_forward is None:
        if verbose:
            util.tic("calculate transition matrices")
        matrices = list(iter_trans_emit_matrices(model, n))
        if verbose:
            util.toc()
        probs_forward = forward_algorithm(model, n, matrices=matrices,
                                          verbose=verbose)
    if probs_backward is None:
        probs_backward = backward_algorithm(model, n, matrices=matrices,
                                            verbose=verbose)

    if matrices:
        delete_trans_emit_matrices(matrices)
            

    total_prob = -util.INF
    for j in xrange(model.get_num_states(0)):
        total_prob = logadd(total_prob,
                            model.prob_prior(0, j) +
                            model.prob_emission(0, j) +
                            probs_backward[0][j])

    probs_post = [
        [probs_forward[i][j] + probs_backward[i][j] - total_prob
         for j in xrange(model.get_num_states(i))]
        for i in xrange(n)]

    return probs_post




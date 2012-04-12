#
# Ancestral Recombination Graph Hidden Markov Model (ArgHmm)
#

from math import exp, log
import random
from itertools import chain, izip
import heapq

from rasmus import hmm, util, stats, treelib
from rasmus.stats import logadd
from compbio import arglib, fasta, phylo


# import arghmm C lib
from arghmm.ctypes_export import *
arghmmc = load_library(["..", "lib"], "libarghmm.so")

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
# export c functions

ex = Exporter(globals())
export = ex.export


if arghmmc:
    # replace python function with c
    
    export(arghmmc, "forward_alg", c_int,
           [c_int, "n", c_int, "nstates",
            c_double_p_p, "trans", c_double_p_p, "emit",
            c_double_matrix, "fw"])

    export(arghmmc, "backward_alg", c_int,
           [c_int, "n", c_int, "nstates",
            c_double_p_p, "trans", c_double_p_p, "emit",
            c_double_matrix, "bw"])

    export(arghmmc, "sample_hmm_posterior", c_int,
           [c_int, "n", c_int, "nstates",
            c_double_p_p, "trans", c_double_p_p, "emit", 
            c_double_matrix, "fw", c_int_list, "path"])


    export(arghmmc, "new_transition_probs", c_double_p_p,
           [c_int, "nnodes", c_int_list, "ptree",
            c_int_list, "ages_index", c_double, "treelen",
            POINTER(c_int * 2), "states", c_int, "nstates",
            c_int, "ntimes", c_double_list, "times",
            c_double_list, "time_steps",
            c_int_list, "nbranches", c_int_list, "nrecombs",
            c_int_list, "ncoals", 
            c_double_list, "popsizes", c_double, "rho"])

    export(arghmmc, "new_transition_probs_switch", c_double_p_p,
           [c_int_list, "ptree", c_int_list, "last_ptree", c_int, "nnodes",
            c_int, "recomb_name", c_int, "recomb_time",
            c_int, "coal_name", c_int, "coal_time",
            c_int_list, "ages_index", c_int_list, "last_ages_index",
            c_double, "treelen", c_double, "last_treelen",
            POINTER(c_int * 2), "states1", c_int, "nstates1",
            POINTER(c_int * 2), "states2", c_int, "nstates2",
            c_int, "ntimes", c_double_list, "times",
            c_double_list, "time_steps",
            c_int_list, "nbranches", c_int_list, "nrecombs",
            c_int_list, "ncoals", 
            c_double_list, "popsizes", c_double, "rho"])

    export(arghmmc, "delete_transition_probs", c_int,
           [c_double_p_p, "transition_probs", c_int, "nstates"])

    export(arghmmc, "new_emissions", c_double_p_p,
           [POINTER(c_int * 2), "states",
            c_int, "nstates", 
            c_int_list, "ptree", c_int, "nnodes", c_int_list, "ages",
            c_char_p_p, "seqs", c_int, "nseqs", c_int, "seqlen",
            c_double_list, "times", c_int, "ntimes",
            c_double, "mu"])

    export(arghmmc, "delete_emissions", c_int,
           [c_double_p_p, "emit", c_int, "seqlen"])


    export(arghmmc, "arghmm_forward_alg", c_double_p_p,
           [c_int_matrix, "ptrees", c_int_matrix, "ages",
            c_int_matrix, "sprs", c_int_list, "blocklens",
            c_int, "ntrees", c_int, "nnodes", 
            c_double_list, "times", c_int, "ntimes",
            c_double_list, "popsizes", c_double, "rho", c_double, "mu",
            c_char_p_p, "seqs", c_int, "nseqs", c_int, "seqlen",
            c_double_p_p, "fw"])

    export(arghmmc, "delete_double_matrix", c_int,
           [c_double_p_p, "mat", c_int, "nrows"])

    export(arghmmc, "arghmm_sample_posterior", POINTER(c_int *2),
           [c_int_matrix, "ptrees", c_int_matrix, "ages",
            c_int_matrix, "sprs", c_int_list, "blocklens",
            c_int, "ntrees", c_int, "nnodes", 
            c_double_list, "times", c_int, "ntimes",
            c_double_list, "popsizes", c_double, "rho", c_double, "mu",
            c_char_p_p, "seqs", c_int, "nseqs", c_int, "seqlen",
            POINTER(POINTER(c_int *2)), "path"])

    export(arghmmc, "delete_path", c_int,
           [POINTER(c_int * 2), "path"])


    export(arghmmc, "get_state_spaces", POINTER(POINTER(c_int * 2)),
           [c_int_matrix, "ptrees", c_int_matrix, "ages",
            c_int_matrix, "sprs", c_int_list, "blocklens",
            c_int, "ntrees", c_int, "nnodes", c_int, "ntimes"])

    export(arghmmc, "delete_state_spaces", c_int,
           [POINTER(POINTER(c_int * 2)), "all_states", c_int, "ntrees"])


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


def discretize_arg(arg, times, ignore_top=True):
    """
    Round node ages to the nearest time point

    If 'ignore_top' is True, then do not use last time point for rounding
    """

    if ignore_top:
        times = times[:-1]
    
    for node in arg:
        i, j = util.binsearch(times, node.age)
        if j is None: j = len(times) - 1
        node.age = times[j]

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
            

def get_treelen(tree, times):
    """Calculate tree length"""
    treelen = sum(x.get_dist() for x in tree)
    rooti = times.index(tree.root.age)
    root_time = times[rooti+1] - times[rooti]
    return treelen + root_time


def get_treelen_branch(tree, times, node, time, treelen=None):

    if treelen is None:
        treelen = sum(x.get_dist() for x in tree)
    else:
        rooti = times.index(tree.root.age)
        root_time = times[rooti+1] - times[rooti]
        treelen -= root_time

    blen = time
    treelen2 = treelen + blen
    if node == tree.root.name:
        treelen2 += blen - tree.root.age
        rooti = times.index(time)
        root_time = times[rooti+1] - times[rooti]
    else:
        rooti = times.index(tree.root.age)
        root_time = times[rooti+1] - times[rooti]
    
    return treelen2 + root_time



#=============================================================================
# helper functions


def parsimony_ancestral_seq(tree, seqs, pos):
    """Calculates ancestral sequence for a local tree using parsimony"""

    ancestral = {}
    sets = {}

    # do unweight parsimony
    for node in tree.postorder():
        if node.is_leaf():
            sets[node] = set([seqs[node.name][pos]])
        else:
            lset = sets[node.children[0]]
            rset = sets[node.children[1]]
            intersect = lset & rset
            if len(intersect) > 0:
                sets[node] = intersect
            else:
                sets[node] = lset | rset

    # traceback
    for node in tree.preorder():
        s = sets[node]
        if len(s) == 1 or not node.parents:
            # NOTE: this technique is used to make assignment deterministic
            ancestral[node.name] = ("A" if "A" in s else
                                    "C" if "C" in s else
                                    "G" if "G" in s else
                                    "T")
        else:
            pchar = ancestral[node.parents[0].name]
            if pchar in s:
                ancestral[node.name] = pchar
            else:
                ancestral[node.name] = ("A" if "A" in s else
                                        "C" if "C" in s else
                                        "G" if "G" in s else
                                        "T")

    return ancestral


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



def sample_recombinations_thread(model, thread, use_times=True):
    """Samples new recombination for a thread"""
    
    r = 0
    
    # assumes that recomb_pos starts with -1 and ends with arg.end
    arg_recomb = model.recomb_pos
    time_lookup = util.list2lookup(model.times)
    minlen = model.time_steps[0]
    
    tree = model.arg.get_marginal_tree(-.5)
    treelen = get_treelen(tree, model.times)
    new_node = model.new_name
    transmat = None
    nstates = 0
    selftrans = None

    next_recomb = -1
    
    for pos, state in enumerate(thread):
        node, node_time = state
        timei = time_lookup[node_time]
        
        # update local tree if needed
        while r < len(arg_recomb) and arg_recomb[r] < pos:
            r += 1
            tree = model.arg.get_marginal_tree(pos-.5)
            treelen = get_treelen(tree, model.times)
            nlineages = get_nlineages_recomb_coal(tree, model.times)
            nbranches, nrecombs, ncoals = nlineages

            if transmat is not None:
                delete_transition_probs(transmat, nstates)
            transmat = calc_transition_probs_c(
                tree, model.states[pos], nlineages,
                model.times, model.time_steps, model.popsizes, model.rho)
            nstates = len(model.states[pos])
            statei = model.states[pos].index((node, timei))
            selftrans = transmat[statei][statei]
            

        if pos == 0 or arg_recomb[r-1] == pos - 1:
            # previous arg recomb is right behind us, sample no recomb
            next_recomb = -1
            continue

        # get information about pos-1
        # since their no recomb in G_{n-1}, last_tree == tree
        last_state = thread[pos-1]
        last_node, last_time = last_state
        last_timei = time_lookup[last_time]
        last_tree = tree
        last_treelen = treelen
        
        if state == last_state:
            if pos > next_recomb:
                # sample the next recomb pos
                last_treelen2 = get_treelen_branch(
                    last_tree, model.times, last_node, last_time)
                rate = max(1.0 - exp(-model.rho * (last_treelen2 - last_treelen)
                                     - selftrans), model.rho)
                next_recomb = pos + int(random.expovariate(rate))
                
            if pos < next_recomb:
                continue


        next_recomb = -1
        last_treelen2 = get_treelen_branch(
            last_tree, model.times, last_node, last_time)
        statei = model.states[pos].index((node, timei))
        selftrans = transmat[statei][statei]

        # there must be a recombination
        # either because state changed or we choose to recombine
        if node == last_node:
            if timei == last_timei:
                # y = v, k in [0, min(timei, last_timei)]
                # y = node, k in Sr(node)
                # if node.parent.age == model.times[timei],
                #   y = sis(last_tree, node.name), k in Sr(y)
                node_timei = time_lookup[tree[node].age]
                recombs = [(new_node, k) for k in
                           range(0, min(timei, last_timei)+1)] + \
                          [(node, k) for k in
                           range(node_timei, min(timei, last_timei)+1)]
            else:
                # y = v, k in [0, min(timei, last_timei)]
                # y = node, k in Sr(node)
                node_timei = time_lookup[tree[node].age]
                recombs = [(new_node, k) for k in
                           range(0, min(timei, last_timei)+1)] + \
                          [(node, k) for k in
                           range(node_timei, min(timei, last_timei)+1)]
        else:
            # y = v, k in [0, min(timei, last_timei)]
            recombs = [(new_node, k)
                       for k in range(0, min(timei, last_timei)+1)]

        if len(recombs) == 0:
            print ((last_node, last_timei), (node, timei))
            raise Exception("recomb not sampled!")

        C = calc_C(model.time_steps, nbranches, model.popsizes)
        j = timei
        probs = []
        for recomb in recombs:
            k = recomb[1]
            probs.append((nbranches[k] + 1) * model.time_steps[k] /
                         (ncoals[j] * (nrecombs[k] + 1.0) * last_treelen2) *
                         (1.0 - exp(-model.time_steps[j] * nbranches[j] /
                                    (2.0 * model.popsizes[j-1]))) *
                         (1.0 - exp(-model.rho * last_treelen2)) *
                         exp(-C[k] + C[k]))
        recomb_node, recomb_time = recombs[stats.sample(probs)]
        
        if use_times:
            recomb_time = model.times[recomb_time]
        yield (pos, recomb_node, recomb_time)

        

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



def add_arg_thread(arg, new_name, thread, recombs, arg3=None):


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
    




def get_clade_point(arg, node_name, time, pos):

    if node_name in arg:
        tree = arg.get_marginal_tree(pos - .5)
        if (time > tree.root.age or
            (time == tree.root.age and node_name not in tree)):
            return (list(tree.leaf_names()), time)
        return (list(tree.leaf_names(tree[node_name])), time)
    else:
        return ([node_name], time)




def sample_arg(seqs, ntimes=20, rho=1.5e-8, mu=2.5e-8, popsize=1e4,
               verbose=False, force=False):
    """
    Sample ARG for sequences
    """
    
    def add_chrom(arg, new_name):
        if verbose:
            util.tic("adding %s..." % new_name)

        model = ArgHmm(arg, seqs, new_name=new_name,
                       times=times, rho=rho, mu=mu, popsize=popsize)

        if verbose:
            util.tic("sample thread")
        path = sample_posterior(model, arg.end)
        if verbose:
            util.toc()

        if verbose:
            util.tic("sample recombs")
        thread = list(iter_thread_from_path(model, path))
        recombs = list(sample_recombinations_thread(
            model, thread))
        if verbose:
            util.toc()

        if verbose:
            util.tic("add thread")
        arg = add_arg_thread(arg, new_name, thread, recombs)
        if verbose:
            util.toc()

        if verbose:
            util.toc()
        return arg, thread

    names = seqs.keys()
    length = len(seqs[names[0]])
    arg = make_trunk_arg(0, length, name=names[0])
    times = get_time_points(ntimes=ntimes, maxtime=80000, delta=.01)

    util.tic("sample ARG of %d sequences" % len(seqs))
    for j in xrange(1, len(names)):
        while True:
            try:
                arg, thread = add_chrom(arg, names[j])
                break
            except:
                if not force:
                    raise
        
    util.toc()
    
    return arg


def resample_arg(arg, seqs, new_name,
                 ntimes=20, rho=1.5e-8, mu=2.5e-8, popsize=1e4,
                 times=None,
                 verbose=False):

    keep = [x for x in seqs.keys() if x != new_name]
    arglib.subarg_by_leaf_names(arg, keep)
    arg = arglib.smcify_arg(arg)
    if times is None:
        times = get_time_points(ntimes=ntimes, maxtime=80000, delta=.01)

    if verbose:
        util.tic("adding %s..." % new_name)
    model = ArgHmm(arg, seqs, new_name=new_name,
                   times=times, rho=rho, mu=mu, popsize=popsize)

    if verbose:
        util.tic("sample thread")
    path = sample_posterior(model, arg.end)
    if verbose:
        util.toc()

    if verbose:
        util.tic("sample recombs")
    thread = list(iter_thread_from_path(model, path))
    recombs = list(sample_recombinations_thread(
        model, thread))
    if verbose:
        util.toc()

    if verbose:
        util.tic("add thread")
    arg = add_arg_thread(arg, new_name, thread, recombs)
    if verbose:
        util.toc()

    if verbose:
        util.toc()

    return arg


def resample_arg_all(arg, seqs, ntimes=20, rho=1.5e-8, mu=2.5e-8, popsize=1e4,
                     times=None, verbose=False):
    if verbose:
        util.tic("resample all chromosomes")
    for new_name in sorted(seqs.keys()):
        arg = resample_arg(arg, seqs, new_name,
                           ntimes=ntimes, rho=rho, mu=mu, popsize=popsize,
                           times=times, verbose=verbose)
    if verbose:
        util.toc()
    return arg
    

def resample_arg_max(arg, seqs, new_name,
                     ntimes=20, rho=1.5e-8, mu=2.5e-8, popsize=1e4,
                     times=None,
                     verbose=False):

    keep = [x for x in seqs.keys() if x != new_name]
    arglib.subarg_by_leaf_names(arg, keep)
    arg = arglib.smcify_arg(arg)
    if times is None:
        times = get_time_points(ntimes=ntimes, maxtime=80000, delta=.01)

    if verbose:
        util.tic("adding %s..." % new_name)
    model = ArgHmm(arg, seqs, new_name=new_name,
                   times=times, rho=rho, mu=mu, popsize=popsize)

    if verbose:
        util.tic("sample thread")
    path = hmm.viterbi(model, arg.end) #sample_posterior(model, arg.end)
    if verbose:
        util.toc()

    if verbose:
        util.tic("sample recombs")
    thread = list(iter_thread_from_path(model, path))
    recombs = list(sample_recombinations_thread(
        model, thread))
    if verbose:
        util.toc()

    if verbose:
        util.tic("add thread")
    arg = add_arg_thread(arg, new_name, thread, recombs)
    if verbose:
        util.toc()

    if verbose:
        util.toc()

    return arg


    
    



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

    ntimes = len(time_steps)
    minlen = time_steps[0]
    treelen = get_treelen(tree, times)
    #treelen = max(sum(x.get_dist() for x in tree), minlen)
    nbranches, nrecombs, ncoals = nlineages
    
    # A_{k,j} =& s'_{j-2} k_{j-1} / (2N) + \sum_{m=k}^{j-2} s'_m k_m / (2N) \\
    #         =& s'_{j-2} k_{j-1} / (2N) + A_{k,j-1}.
    
    A = util.make_matrix(ntimes, ntimes, 0.0)
    for k in xrange(ntimes):
        # A[k][k] = 0
        for j in xrange(k+1, ntimes):
            l = j - 1
            A[k][j] = A[k][j-1] + time_steps[l] * nbranches[l] / (2.0 * popsizes[l])

    # B_{c,a} =& \sum_{k=0}^{c} \exp(- A_{k,a}) \\
    #         =& B_{c-1,a} + \exp(- A_{c,a}).

    B = util.make_matrix(ntimes, ntimes, 0.0)
    for b in xrange(ntimes):
        B[0][b] = nbranches[0] * time_steps[0] / (nrecombs[0] + 1.0) * exp(-A[0][b])
        for c in xrange(1, min(b+1, ntimes-1)):
            B[c][b] = (B[c-1][b] + nbranches[c] * time_steps[c] / (nrecombs[c] + 1.0)
                       * exp(-A[c][b]))

    # S_{a,b} &= B_{min(a,b),b}
    S = util.make_matrix(ntimes, ntimes, 0.0)
    for a in xrange(ntimes):
        for b in xrange(ntimes):
            S[a][b] = B[min(a, b)][b]

    # f =\frac{[1 - \exp(- \rho (|T^{n-1}_{i-1}| + s_a))] 
    #       [1 - \exp(- s'_b k_b / (2N))]}
    #      {\exp(-\rho |T^{n-1}_{i-1}|) (|T^{n-1}_{i-1}| + s_a) k^C_b}
    # |T^{n-1}_{i-1}| = treelen
    
    time_lookup = util.list2lookup(times)
    transprob = util.make_matrix(len(states), len(states), 0.0)
    for i, (node1, a) in enumerate(states):
        c = time_lookup[tree[node1].age]

        blen = times[a]
        treelen2 = treelen + blen
        if node1 == tree.root.name:
            treelen2 += blen - tree.root.age
            treelen2 += time_steps[a]
        else:
            treelen2 += time_steps[time_lookup[tree.root.age]]
        
        for j, (node2, b) in enumerate(states):
            f = ((1.0 - exp(-rho * treelen2)) *
                 (1.0 - exp(-time_steps[b] * nbranches[b]
                            / (2.0 * popsizes[b]))) /
                 (exp(-rho * treelen) * treelen2 * ncoals[b]))
            
            if node1 != node2:
                transprob[i][j] = f * S[a][b]
            else:
                transprob[i][j] = f * (2*S[a][b] - S[c][b])
                if a == b:
                    transprob[i][j] += exp(-rho * (treelen2 - treelen))

        # normalize row to sum to one
        tot = sum(transprob[i])
        for j in xrange(len(states)):
            transprob[i][j] = util.safelog(transprob[i][j] / tot)
        #for j in xrange(len(states)):
        #    transprob[i][j] = util.safelog(transprob[i][j])

    return transprob



def calc_transition_probs_c(tree, states, nlineages, times,
                            time_steps, popsizes, rho, raw=True):
    
    nbranches, nrecombs, ncoals = nlineages

    times_lookup = dict((t, i) for i, t in enumerate(times))
    tree2 = tree.get_tree()
    ptree, nodes, nodelookup = make_ptree(tree2)
    int_states = [[nodelookup[tree2[node]], timei]
                  for node, timei in states]
    nstates = len(int_states)
    ages_index = [times_lookup[tree[node.name].age]
                  for node in nodes]
    #treelen = sum(x.dist for x in tree2)
    treelen = get_treelen(tree, times)
    transmat = new_transition_probs(
        len(nodes), ptree, ages_index, treelen, 
        ((c_int * 2) * nstates)
        (* ((c_int * 2)(n, t) for n, t in int_states)), nstates,
        len(time_steps), times, time_steps,
        nbranches, nrecombs, ncoals, 
        popsizes, rho)

    if raw:
        return transmat
    else:
        transmat2 = [transmat[i][:nstates]
            for i in range(nstates)]
        delete_transition_probs(transmat, nstates)
        return transmat2


def calc_transition_probs_switch_c(tree, last_tree, recomb_name,
                                   states1, states2,
                                   nlineages, times,
                                   time_steps, popsizes, rho, raw=True):

    times_lookup = dict((t, i) for i, t in enumerate(times))
    nbranches, nrecombs, ncoals = nlineages
    (recomb_branch, recomb_time), (coal_branch, coal_time) = \
        find_recomb_coal(tree, last_tree, recomb_name=recomb_name)
    
    recomb_time = times.index(recomb_time)
    coal_time = times.index(coal_time)

    last_tree2 = last_tree.copy()
    arglib.remove_single_lineages(last_tree2)
    tree2 = tree.copy()
    arglib.remove_single_lineages(tree2)

    # get last ptree
    last_tree2 = last_tree2.get_tree()
    tree2 = tree2.get_tree()
    last_ptree, last_nodes, last_nodelookup = make_ptree(last_tree2)

    # find old node and new node
    recomb_parent = last_tree2[recomb_branch].parent
    recoal = [x for x in tree2 if x.name not in last_tree2][0]

    # make nodes array consistent
    nodes = [tree2.nodes.get(x.name, None) for x in last_nodes]
    i = last_nodes.index(recomb_parent)
    assert nodes[i] == None
    nodes[i] = recoal

    # get ptree
    ptree, nodes, nodelookup = make_ptree(tree2, nodes=nodes)

    # get recomb and coal branches
    recomb_name = last_nodelookup[last_tree2[recomb_branch]]
    coal_name = last_nodelookup[last_tree2[coal_branch]]
    
    int_states1 = [[last_nodelookup[last_tree2[node]], timei]
                  for node, timei in states1]
    nstates1 = len(int_states1)
    int_states2 = [[nodelookup[tree2[node]], timei]
                  for node, timei in states2]
    nstates2 = len(int_states2)
    
    last_ages_index = [times_lookup[last_tree[node.name].age]
                       for node in last_nodes]
    ages_index = [times_lookup[tree[node.name].age]
                  for node in nodes]

    last_treelen = sum(x.dist for x in last_tree2)
    treelen = sum(x.dist for x in tree2)
    
    transmat = new_transition_probs_switch(
        ptree, last_ptree, len(nodes),
        recomb_name, recomb_time, coal_name, coal_time,

        ages_index, last_ages_index,
        treelen, last_treelen,
        ((c_int * 2) * nstates1)
        (* ((c_int * 2)(n, t) for n, t in int_states1)), nstates1, 
        ((c_int * 2) * nstates2)
        (* ((c_int * 2)(n, t) for n, t in int_states2)), nstates2,
        
        len(time_steps), times, time_steps,
        nbranches, nrecombs, ncoals, 
        popsizes, rho)

    if raw:
        return transmat
    else:
        transmat2 = [transmat[i][:nstates2]
            for i in range(nstates1)]
        delete_transition_probs(transmat, nstates1)
        return transmat2


def find_recomb_coal(tree, last_tree, recomb_name=None, pos=None):

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
        if (node1, a) != (coal_branch, coal_time):
            # deterministic transition
            assert determ[i] != -1, determ
            transprob[i][determ[i]] = 0.0

        else:
            # probabilistic transition case

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
                #node2 = walk_up(node1, node1, a)
                #next_states.append(state2_lookup[(node2, a)])
                start = node1
                #ignore = None

            else:
                child1 = node.children[0]
                child2 = node.children[1]
                
                if recomb_branch == child1.name:
                    # right child is not disrupted
                    #node2 = walk_up(node1, child2.name, a, node1)
                    #next_states.append(state2_lookup[(node2, a)])
                    start = child2.name
                    #ignore = node1

                elif recomb_branch == child2.name:
                    # left child is not disrupted
                    #node2 = walk_up(node1, child1.name, a, node1)
                    #next_states.append(state2_lookup[(node2, a)])
                    start = child1.name
                    #ignore = node1

                else:
                    # node is not disrupted
                    #node2 = walk_up(node1, node1, a)
                    #next_states.append(state2_lookup[(node2, a)])
                    start = node1
                    #ignore = None

            #assert ignore not in tree

            # optionally walk up
            if ((coal_branch == node1 or coal_branch == start) and
                coal_time < a):
                # coal occurs under us
                # TODO: make this probabilistic
                ptr = tree[start].parents[0]
                node2 = ptr.name
            else:
                node2 = start
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

    priormat = [
        log((1 - exp(- time_steps[b] * nlineages[0][b] /
                 (2.0 * popsizes[b]))) / nlineages[2][b] *
             exp(-sum(time_steps[m] * nlineages[0][m] /
                      (2.0 * popsizes[m])
                      for m in range(0, b))))
            for node, b in states]
    
    return priormat



def arg_lca(arg, leaves, time, pos, ignore=None):

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



def iter_arg_sprs(arg, start=None, end=None):

    if start is None:
        start = arg.start
    if end is None:
        end = arg.end

    last_tree_full = None
    last_tree = None
    for block, tree_full in arglib.iter_tree_tracks(arg, start, end):
        if last_tree_full:
            recomb = (x for x in tree_full if x.pos == block[0]).next()
            spr = find_recomb_coal(tree_full, last_tree_full,
                                   recomb_name=recomb.name)
        else:
            spr = None
        
        tree = tree_full.copy()
        tree = arglib.remove_single_lineages(tree)

        # convert block to our system
        a, b = block
        if a == 0:
            a = -1
        if b == end:
            b -= 1
        block = [a+1, b+1]

        #print "-----------"
        #print spr
        #if last_tree:
        #    last_tree.get_tree().write()
        #tree.get_tree().write()
        
        yield block, tree, last_tree, spr

        last_tree_full = tree_full
        last_tree = tree


def get_treeset(arg, times, start=None, end=None):

    times_lookup = dict((t, i) for i, t in enumerate(times))

    ptrees  = []
    ages = []
    sprs = []
    blocks = []
    all_nodes = []

    for block, tree, last_tree, spr in iter_arg_sprs(arg, start, end):
        
        if last_tree is None:
            # get frist ptree
            tree2 = tree.get_tree()
            ptree, nodes, nodelookup = make_ptree(tree2)
            ispr = [-1, -1, -1, -1]
            age = [times_lookup[tree[x.name].age] for x in nodes]

        else:
            tree2 = tree.get_tree()
            (rname, rtime), (cname, ctime) = spr
            
            # find old node and new node
            recomb_parent = last_tree2[rname].parent
            recoal = [x for x in tree2 if x.name not in last_tree2][0]

            # make nodes array consistent
            nodes = [tree2.nodes.get(x.name, None) for x in last_nodes]
            i = last_nodes.index(recomb_parent)
            assert nodes[i] is None
            nodes[i] = recoal

            # get ptree
            ptree, nodes, nodelookup = make_ptree(tree2, nodes=nodes)
            age = [times_lookup[tree[x.name].age] for x in nodes]

            # get integer-based spr
            recomb_name = last_nodelookup[last_tree2[rname]]
            coal_name = last_nodelookup[last_tree2[cname]]
            ispr = [recomb_name, times_lookup[rtime],
                    coal_name, times_lookup[ctime]]

        # append integer-based data
        ptrees.append(ptree)
        ages.append(age)
        sprs.append(ispr)
        blocks.append(block)
        all_nodes.append([x.name for x in nodes])

        # setup last tree
        last_tree = tree
        last_tree2 = tree2
        last_ptree, last_nodes, last_nodelookup = ptree, nodes, nodelookup

    return (ptrees, ages, sprs, blocks), all_nodes

        

#=============================================================================
# trunk ARG

def make_trunk_arg(start, end, name="ind1"):
    
    arg = arglib.ARG(start=start, end=end)
    node = arg.new_node(name, event="gene", age=0)
    return arg


def make_single_seq(length, name="ind1"):

    dna = "ACGT"
    seqs = fasta.FastaDict()
    seqs[name] = "".join(dna[random.randint(0, 3)]
                         for i in xrange(length))
    return seqs






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
        self.local_site = None
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
            self.priormat = [
                log((1 - exp(- self.time_steps[b-1] * self.nlineages[0][b-1] /
                         (2.0 * self.popsizes[b-1]))) / self.nlineages[2][b] *
                     exp(-sum(self.time_steps[m] * self.nlineages[0][m] /
                              (2.0 * self.popsizes[m])
                              for m in range(0, b-1))))
                for node, b in self.states[pos]]

            # makes computing emissions easier
            arglib.remove_single_lineages(self.local_tree)


        # update local site
        if force or pos != self.last_pos:
            self.local_site = parsimony_ancestral_seq(
                self.local_tree, self.seqs, pos)
            

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
        node_name, timei = self.states[pos][state]
        node = self.local_tree[node_name]
        time = self.times[timei]
        mu = self.mu

        mintime = self.time_steps[0]

        # v = new chromosome
        # x = current branch
        # p = parent of current branch

        #if state == 0:
        #    ptree, nodes, node_lookup = make_ptree(self.local_tree.get_tree())
        #    print pos, "".join(self.local_site[x.name] for x in nodes), "'"

        if node.parents:
            parent = node.parents[0]
            parent_age = parent.age

            #print "parent.parents", parent.parents

            if not parent.parents:
                # unwrap top branch
                c = parent.children
                sib = (c[1] if node == c[0] else c[0])
                
                v = self.seqs[self.new_name][pos]
                x = self.local_site[node.name]
                p = self.local_site[sib.name]

                # modify (x,p) length to (x,p) + (sib,p)
                parent_age = 2 * parent_age - sib.age

            else:
                v = self.seqs[self.new_name][pos]
                x = self.local_site[node.name]
                p = self.local_site[parent.name]

        else:
            #print "parent", None
            parent = None
            parent_age = None

            # adjust time by unwrapping branch
            time = 2 * time - node.age

            v = self.seqs[self.new_name][pos]
            x = self.local_site[node.name]
            p = x

        time = max(time, mintime)

        #print " ", pos, state, v, x, p, "'"

        if v == x == p:
            # no mutation
            return - self.mu * time

        elif v != p == x:
            # mutation on v
            return log(.33 - .33 * exp(-mu * time))

        elif v == p != x:
            # mutation on x
            t1 = max(parent_age - node.age, mintime)
            t2 = max(time - node.age, mintime)

            return log((1 - exp(-mu *t2)) / (1 - exp(-mu * t1))
                       * exp(-mu * (time + t2 - t1)))

        elif v == x != p:
            # mutation on (y,p)
            t1 = max(parent_age - node.age, mintime)
            t2 = max(parent_age - time, mintime)

            return log((1 - exp(-mu * t2)) / (1 - exp(-mu * t1))
                       * exp(-mu * (time + t2 - t1)))

        else:
            # two mutations (v,x)
            # mutation on x
            if parent:
                t1 = max(parent_age - node.age, mintime)
                t2a = max(parent_age - time, mintime)
            else:
                t1 = max(self.times[-1] - node.age, mintime)
                t2a = max(self.times[-1] - time, mintime)
            t2b = max(time - node.age, mintime)
            t2 = max(t2a, t2b)
            t3 = time

            return log((1 - exp(-mu *t2)) * (1 - exp(-mu *t3))
                       / (1 - exp(-mu * t1))
                       * exp(-mu * (time + t2 + t3 - t1)))


    def emit(self, pos, state):

        self.check_local_tree(pos)
        node_name, timei = self.states[pos][state]
        time = self.times[timei]
        base = self.local_site[node_name]

        # sample whether to mutation from an exponential distrib
        if random.expovariate(self.mu) < time:
            while True:
                x = "ACGT"[random.randint(0, 3)]
                if x != base:
                    return x
        else:
            return base
        
            
            



def arghmm_sim(arg, seqs, name=None, times=None,
               ntimes=30, maxtime=45000.0, delta=.01):

    model = ArgHmm(arg, seqs, name=name, times=times,
                   ntimes=ntimes, maxtime=maxtime, delta=delta)
    


    
#=============================================================================
# custom HMM methods



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
            
            transmat_switch = calc_transition_probs_switch(
                tree, last_tree, recomb.name,
                model.states[start-1], model.states[start],
                last_nlineages, model.times,
                model.time_steps, model.popsizes, model.rho)

            #transmat_switch = calc_transition_probs_switch_c(
            #    tree, last_tree, recomb.name,
            #    model.states[start-1], model.states[start],
            #    last_nlineages, model.times,
            #    model.time_steps, model.popsizes, model.rho, raw=False)
            
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


def delete_trans_emit_matrices(matrices):
    """Delete matrices"""
    for block, nstates, transmat, transmat_switch, emit in matrices:
        delete_emissions(emit, block[1] - block[0])
        delete_transition_probs(transmat, nstates)
        

    

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
    fw = arghmm_forward_alg(ptrees, ages, sprs, blocklens,
                            len(ptrees), len(ptrees[0]), 
                            model.times, len(model.times),
                            model.popsizes, model.rho, model.mu,
                            (c_char_p * len(seqs))(*seqs), len(seqs),
                            seqlen, None)

    # map states to python state space
    all_states = get_state_spaces(ptrees, ages, sprs, blocklens,
                                  len(ptrees), len(ptrees[0]), len(model.times))

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


def sample_posterior(model, n, probs_forward=None,
                     verbose=False, matrices=None):

    if verbose:
        util.tic("sample thread")

    (ptrees, ages, sprs, blocks), all_nodes = get_treeset(
        model.arg, model.times)
    blocklens = [x[1] - x[0] for x in blocks]
    seqlen = sum(blocklens)

    seqs = [model.seqs[node] for node in all_nodes[0]
            if model.arg[node].is_leaf()]
    seqs.append(model.seqs[model.new_name])
    path = arghmm_sample_posterior(
        ptrees, ages, sprs, blocklens,
        len(ptrees), len(ptrees[0]), 
        model.times, len(model.times),
        model.popsizes, model.rho, model.mu,
        (c_char_p * len(seqs))(*seqs), len(seqs),
        seqlen, None)
    
    path2 = []
    for k, (start, end) in enumerate(blocks):
        states = model.states[start]
        lookup = util.list2lookup(states)
        nodes = all_nodes[k]
        
        for i in range(start, end):
            s = (nodes[path[i][0]], path[i][1])
            path2.append(lookup[s])

    delete_path(path)

    if verbose:
        util.toc()

    return path2



def sample_posterior2(model, n, probs_forward=None, matrices=None,
                     verbose=False):

    # get transition matrices
    matrices_given = matrices is not None
    if matrices is None:
        matrices = list(iter_trans_emit_matrices(model, n))

    # get forward table
    if probs_forward:
        probs = probs_forward
    else:
        probs = forward_algorithm(model, n, matrices=matrices, verbose=verbose)

    if verbose:
        util.tic("sample thread")

    # choose last column first
    path = range(n)
    i = n-1
    total = stats.logsum(probs[-1])
    path[i] = stats.sample([exp(x - total) for x in probs[-1]])
    
    for block, nstates, transmat, transmat_switch, emit in reversed(matrices):
        #if verbose:
        #    util.logger(" pos %d" % block[0])
        blocklen = block[1] - block[0]
        # use transmat and sample path for block

        # make fw and path inputs
        fw = []
        for i in xrange(block[0], block[1]):
            fw.append(probs[i])
        path2 = range(block[0], block[1])
        path2[-1] = path[block[1]-1]
        
        sample_hmm_posterior(blocklen, nstates, transmat, emit, fw, path2)
        
        # get path output
        for i in xrange(block[0], block[1]-1):
            path[i] = path2[i-block[0]]


        # use switch matrix for last col of next block
        if block[0] > 0:
            nstates1 = len(transmat_switch)
            nstates2 = len(transmat_switch[0])
            assert path[block[0]] < nstates2
            
            i = block[0] - 1
            A = []
            k = path[i+1]
            for j in range(nstates1):
                A.append(probs[i][j] + transmat_switch[j][k])
            tot = stats.logsum(A)
            path[i] = stats.sample([exp(x - tot) for x in A])

        if not matrices_given:
            delete_emissions(emit, blocklen)
            delete_transition_probs(transmat, nstates)

    if verbose:
        util.toc()
            
    return path




def iter_forward_algorithm(model, n, verbose=False):

    # NOTE: not currently used
    
    # calc first position
    nstates = model.get_num_states(0)
    col1 = [model.prob_prior(0, j) + model.prob_emission(0, j)
            for j in xrange(nstates)]
    
    if n > 20:
        step = (n // 20)
    else:
        step = 1
    
    # loop through positions
    nstates1 = nstates
    i = 1
    next_print = step
    while i < n:
        while verbose and i > next_print:
            next_print += step
            print " forward iter=%d/%d" % (i+1, n)
        
        nstates2 = model.get_num_states(i)
        col2 = [0] * nstates2
        emit = [model.prob_emission(i, k) for k in xrange(nstates2)]
        trans = [[model.prob_transition(i-1, j, i, k)
                  for j in xrange(nstates1)]
                 for k in xrange(nstates2)]
        forward_step(i, col1, col2, nstates1, nstates2, trans, emit)
        
        yield col2
        col1 = col2
        nstates1 = nstates2
        i += 1



#=============================================================================
# python equivalent functions


def py_forward_algorithm(model, n, verbose=False):

    probs = []

    # calc first position
    nstates = model.get_num_states(0)
    probs.append([model.prob_prior(0, j) + model.prob_emission(0, j)
                  for j in xrange(nstates)])
    
    if n > 20:
        step = (n // 20)
    else:
        step = 1
    
    # loop through positions
    nstates1 = nstates
    i = 1
    next_print = step    
    while i < n:
        while verbose and i > next_print:
            next_print += step
            print " forward iter=%d/%d" % (i+1, n)

        # do first position manually
        nstates2 = model.get_num_states(i)
        model.check_local_tree(i)
        if i == model.local_block[0] and model.transmat_switch:
            trans = model.transmat_switch
        else:
            trans = model.transmat
        
        col1 = probs[i-1]

        # find total transition and emission
        col2 = []
        for k in xrange(nstates2):
            tot = -util.INF
            emit = model.prob_emission(i, k)
            for j in xrange(nstates1):
                p = col1[j] + trans[j][k] + emit
                tot = logadd(tot, p)
            col2.append(tot)
                
        probs.append(col2)
        nstates1 = nstates2
        i += 1
        if i >= n:
            break

        # do rest of block quickly
        space = model.get_state_space(i)
        block = model.get_local_block(space)
        blocklen = block[1] - i

        if i > block[0] and blocklen > 4:
            nstates = model.get_num_states(i)

            # setup tree and states
            tree = model.arg.get_marginal_tree(i-.5)
            tree2 = tree.get_tree()
            ptree, nodes, nodelookup = make_ptree(tree2)
            int_states = [[nodelookup[tree2[node]], timei]
                          for node, timei in model.states[i]]
            ages = [model.times.index(tree[node.name].age) for node in nodes]
            seqs = [model.seqs[node.name][i-1:block[1]]
                    for node in nodes if node.is_leaf()]
            seqs.append(model.seqs[model.new_name][i-1:block[1]])
            seqlen = blocklen + 1
            
            emit = new_emissions(
                ((c_int * 2) * nstates)
                (* ((c_int * 2)(n, t) for n, t in int_states)), nstates, 
                ptree, len(ptree), ages,
                (c_char_p * len(seqs))(*seqs), len(seqs), seqlen,
                model.times, len(model.times), model.mu)

            trans = c_matrix(
                c_double,
                [[model.prob_transition(i-1, j, i, k)
                  for k in xrange(nstates)] for j in xrange(nstates)])
            
            fw = [probs[-1]]
            for pos in xrange(i, block[1]):
                fw.append([0.0 for k in xrange(nstates)])
                
            forward_alg(blocklen+1, nstates, trans, emit, fw)

            delete_emissions(emit, blocklen)

            for col in fw[1:]:
                probs.append(col[:nstates])
            nstates1 = nstates
            i = block[1]
            
    return probs


def py_forward_algorithm2(model, n, verbose=False, matrices=None):

    probs = []

    if verbose:
        util.tic("forward")

    # get prior matrix
    local_tree = model.arg.get_marginal_tree(-.5)
    nlineages = get_nlineages_recomb_coal(local_tree, model.times)
    priors = calc_state_priors(
        local_tree, model.states[0], nlineages,
        model.times, model.time_steps, model.popsizes, model.rho)
    probs.append(priors)

    matrices_given = matrices is not None
    if matrices is None:
        matrices = iter_trans_emit_matrices(model, n)

    # iterate over blocks
    for block, nstates, transmat, transmat_switch, emit in matrices:
        if verbose:
            util.logger(" pos %d %d" % (block[0], block[1] - block[0]))

        blocklen = block[1] - block[0]

        # use switch matrix for first col
        if block[0] > 0:
            nstates1 = len(transmat_switch)
            nstates2 = len(transmat_switch[0])
            
            col1 = probs[-1]
            col2 = []
            for k in xrange(nstates2):
                e = emit[0][k]
                col2.append(stats.logsum([col1[j] + transmat_switch[j][k] + e
                                          for j in xrange(nstates1)]))
            probs.append(col2)

        # use transmat for rest of block
        # make forward table for block
        fw = [probs[-1]]
        for pos in xrange(block[0]+1, block[1]):
            fw.append([0.0 for k in xrange(nstates)])
        
        forward_alg(blocklen, nstates, transmat, emit, fw)

        if not matrices_given:
            delete_emissions(emit, blocklen)
            delete_transition_probs(transmat, nstates)
        
        for col in fw[1:]:
            probs.append(col[:nstates])


    if verbose:
        util.toc()
            
    return probs



def py_sample_posterior(model, n, probs_forward=None, verbose=False):

    # NOTE: logsum is used for numerical stability

    path = range(n)

    # get forward probabilities
    if probs_forward is None:
        probs_forward = py_forward_algorithm(model, n, verbose=verbose)

    # base case i=n-1
    i = n-1
    A = [probs_forward[i][j] for j in range(model.get_num_states(i))]
    tot = stats.logsum(A)
    path[i] = stats.sample([exp(x - tot) for x in A])
    #path[i] = stats.sample(map(exp, A))
  
    # recurse
    for i in xrange(n-2, -1, -1):
        C = []
        A = []
        for j in range(model.get_num_states(i)):
            # C_{i,j} = trans(j, Y[i+1]) * emit(X[i+1], Y[i+1])
            # !$A_{j,i} = F_{i,j} C_{i,j}$!
            C.append(
                model.prob_transition(i, j, i+1, path[i+1]) +
                model.prob_emission(i+1, path[i+1]))
            A.append(probs_forward[i][j] + C[j])
        tot = stats.logsum(A)
        path[i] = j = stats.sample([exp(x - tot) for x in A])
        #path[i] = j = stats.sample(map(exp, A))
    
    return path



def sample_posterior_old(model, n, probs_forward=None, verbose=False):

    # NOTE: logsum is used for numerical stability

    path = range(n)

    # get forward probabilities
    if probs_forward is None:
        probs_forward = forward_algorithm(model, n, verbose=verbose)

    # base case i=n-1
    B = 0.0
    i = n-1
    A = [probs_forward[i][j] for j in range(model.get_num_states(i))]
    tot = stats.logsum(A)
    path[i] = stats.sample([exp(x - tot) for x in A])
    #path[i] = stats.sample(map(exp, A))
  
    # recurse
    for i in xrange(n-2, -1, -1):
        C = []
        A = []
        for j in range(model.get_num_states(i)):
            # C_{i,j} = trans(j, Y[i+1]) * emit(X[i+1], Y[i+1])
            # !$A_{j,i} = F_{i,j} C_{i,j} B_{i+1,l}$!
            C.append(
                model.prob_transition(i, j, i+1, path[i+1]) +
                model.prob_emission(i+1, path[i+1]))
            A.append(probs_forward[i][j] + C[j] + B)
        #tot = stats.logsum(A)
        path[i] = j = stats.sample([exp(x - tot) for x in A])
        #path[i] = j = stats.sample(map(exp, A))
        # !$B_{i,j} = C_{i,j} B_{i+1,l}$!
        B += C[j]
    
    return path


#=============================================================================
# C interface functions

def make_ptree(tree, skip_single=True, nodes=None):
    """Make parent tree array from tree"""

    ptree = []

    if nodes is None:
        nodes = []
        if skip_single:
            nodes = list(x for x in tree.postorder() if len(x.children) != 1)
        else:
            nodes = list(tree.postorder())
        assert nodes[-1] == tree.root
    
    # ensure sort is stable
    def leafsort(a, b):
        if a.is_leaf():
            if b.is_leaf():
                return 0
            else:
                return -1
        else:
            if b.is_leaf():
                return 1
            else:
                return 0
    
    # bring leaves to front
    nodes.sort(cmp=leafsort)
    
    # make lookup
    nodelookup = {}
    for i, n in enumerate(nodes):
        nodelookup[n] = i

    # make ptree
    for node in nodes:
        if node == tree.root:
            ptree.append(-1)
        else:
            parent = node.parent
            if skip_single:
                while len(parent.children) == 1:
                    parent = parent.parent
            ptree.append(nodelookup[parent])
        
    return ptree, nodes, nodelookup





#=============================================================================
# OLD CODE


'''

def get_deterministic_transitions(states1, states2, times, tree, last_tree,
                                  recomb_branch, recomb_time,
                                  coal_branch, coal_time):

    # recomb_branch in tree
    # coal_branch in last_tree
    
    # get leaves under recomb_node
    recomb_leaves = set(last_tree.leaves(last_tree[recomb_branch]))

    def find_state(node, time):
        b = util.INF
        state2 = None
        
        while len(node.children) == 1:
            node = node.children[0]
        
        for j, (n, t) in enumerate(states2):
            if node.name == n and time <= t < b:
                b = t
                state2 = j
        assert state2 is not None, ((node, time), states2)
        return state2
                

    def trace_up(node, time):
        last = node
        while node.age <= times[time]:
            if len(node.children) != 1:
                last = node
            if not node.parents:
                break
            node = node.parents[0]
        return last

    next_states = []
    for i, state1 in enumerate(states1):
        node1, a = state1
        leaves1 = set(last_tree.leaves(last_tree[node1]))
        remain = leaves1 - recomb_leaves

        if (node1, a) == (coal_branch, coal_time):
            # not a deterministic case (just mark i-->i)
            next_states.append(i)
        
        elif len(remain) > 0:
            # SPR only removes a subset of descendents
            # trace up from remaining leaf to find correct new state
            ptr = tree[iter(remain).next().name]
            node = trace_up(ptr, a)
            next_states.append(find_state(node, a))

        else:
            # SPR is on same branch as new chromosome
            if recomb_time >= a:
                # we move with SPR subtree
                ptr = tree[iter(recomb_leaves).next().name]
                node = trace_up(ptr, a)
                next_states.append(find_state(node, a))

            elif coal_time <= a and coal_branch == node1:
                # SPR subtree coals back underneath us
                next_states.append(find_state(tree[node1], a))

            else:
                # SPR subtree moves out from underneath us
                # therefore therefore the new chromosome coalesces with
                # the branch above the subtree

                # search up for parent
                ptr = last_tree[recomb_branch]
                ptr = ptr.parents[0]
                while len(ptr.children) == 1:
                    ptr = ptr.parents[0]
                b = times.index(ptr.age)

                # go over to new tree
                if ptr.name not in tree:
                    # we are above root
                    assert ptr.age >= tree.root.age
                    next_states.append(find_state(tree.root, b))
                else:
                    ptr = tree[ptr.name]
                    next_states.append(find_state(ptr, b))

    return next_states

'''

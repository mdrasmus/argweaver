#
# Ancestral Recombination Graph Hidden Markov Model (ArgHmm)
#

from math import exp, log
import random

from rasmus import hmm, util, stats, treelib
from rasmus.stats import logadd
from compbio import arglib, fasta


# import arghmm C lib
from arghmm.ctypes_export import *
arghmmc = load_library(["..", "lib"], "libarghmm.so")

#=============================================================================
# constants

PROGRAM_NAME = u"ArgHmm"
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
    export(arghmmc, "forward_step", c_int, 
           [c_int, "i", c_double_list, "col1", c_double_list, "col2",
            c_int, "nstates1", c_int, "nstates2",
            c_double_matrix, "trans", c_double_list, "emit"])

    export(arghmmc, "forward_alg", c_int,
           [c_int, "n", c_int, "nstates1", c_int, "nstates2",
            c_double_matrix, "fw", c_double_p_p, "trans",
            c_double_p_p, "emit"])

    export(arghmmc, "backward_alg", c_int,
           [c_int, "n", c_int, "nstates1", c_int, "nstates2",
            c_double_matrix, "fw", c_double_p_p, "trans",
            c_double_p_p, "emit"])

    export(arghmmc, "new_emissions", c_double_p_p,
           [POINTER(c_int * 2), "states",
            c_int, "nstates", 
            c_int_list, "ptree", c_int, "nnodes", c_double_list, "ages",
            c_char_p_p, "seqs", c_int, "nseqs", c_int, "seqlen",
            c_double_list, "times", c_int, "ntimes",
            c_double, "mu"])

    export(arghmmc, "delete_emissions", c_int,
           [c_double_p_p, "emit"])



#=============================================================================
# discretization


def get_time_point(i, ntimes, maxtime, delta=10):
    return (exp(i/float(ntimes) * log(1 + delta * maxtime)) - 1) / delta


def get_time_points(ntimes=30, maxtime=45000, delta=.01):
    return [get_time_point(i, ntimes, maxtime, delta)
            for i in range(ntimes+1)]


def iter_coal_states(tree, times):

    seen = set()
    
    for node in tree.preorder():
        if len(node.children) == 1:
            continue
        i, j = util.binsearch(times, node.age)
        
        if node.parents:
            parent = node.parents[0]
            while parent not in seen:
                parent = parent.parents[0]
            
            # do not offer coalescing at bottom of branch
            # if branch is non-zero len
            #if parent.age > node.age:
            #    i += 1
            while i < len(times) and times[i] <= parent.age:
                yield (node.name, i)
                i += 1
        else:
            # do not coalesce at bottom of root branch, unless root of tree
            # is at top time
            #if i < len(times) - 1:
            #    i += 1
            while i < len(times):
                yield (node.name, i)
                i += 1

        seen.add(node)


def get_nlineages(tree, times):
    """Count the number of lineages in each time segment"""
    nlineages = [0 for i in times]
    for name, i in iter_coal_states(tree, times):
        node = tree[name]
        if node.parents:
            parent = node.parents[0]
            while len(parent.children) == 1:
                parent = parent.parents[0]
        if not node.parents or parent.age > node.age:
            nlineages[i-1] += 1
    nlineages[-1] = 1
    return nlineages


def discretize_arg(arg, times):
    for node in arg:
        i, j = util.binsearch(times, node.age)
        if j is None: j = len(times) - 1
        node.age = times[j]

        if node.event == "recomb":
            node.pos = int(node.pos)



#=============================================================================
# helper functions

def parsimony_ancestral_seq(tree, seqs, pos):

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
            ancestral[node.name] = s.pop()
        else:
            pchar = ancestral[node.parents[0].name]
            if pchar in s:
                ancestral[node.name] = pchar
            else:
                ancestral[node.name] = s.pop()

    return ancestral



def iter_chrom_thread(arg, node, by_block=True):

    for block, tree in arglib.iter_tree_tracks(arg):

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

        if by_block:
            yield (sib.name, parent.age, block)
        else:
            for i in range(int(block[0]+1), int(block[1]+1)):
                yield (sib.name, parent.age)


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



def add_arg_thread(arg, new_name, thread, by_block=False):
    pass



#=============================================================================
# probabilities


def calc_transition_probs(tree, states, nlineages, times,
                          time_steps, popsizes, rho):

    ntimes = len(time_steps)
    treelen = sum(x.get_dist() for x in tree)
    mintime = time_steps[0]
    
    # A_{k,j} =& s'_{j-2} k_{j-2} / (2N) + \sum_{m=k}^{j-3} s'_m k_m / (2N) \\
    #         =& s'_{j-2} k_{j-2} / (2N) + A_{k,j-1}.
    
    A = util.make_matrix(ntimes, ntimes, 0.0)
    for k in xrange(ntimes):
        # A[k][k] = A[k][k+1] = 0
        for j in xrange(k+2, ntimes):
            l = j - 2
            A[k][j] = A[k][j-1] + time_steps[l] * nlineages[l] / (2.0 * popsizes[l])

    # B_{c,a} =& \sum_{k=0}^{c} \exp(- A_{k,a}) \\
    #         =& B_{c-1,a} + \exp(- A_{c,a}).

    B = util.make_matrix(ntimes, ntimes, 0.0)
    for b in xrange(ntimes):
        B[0][b] = time_steps[0] * exp(-A[0][b])
        for c in xrange(1, b):
            B[c][b] = B[c-1][b] + time_steps[c] * exp(-A[c][b])


    # S_{a,b} &= B_{min(a-1,b-1),a}
    
    S = util.make_matrix(ntimes, ntimes, 0.0)
    for a in xrange(1, ntimes):
        for b in xrange(1, ntimes):
            S[a][b] = B[min(a-1, b-1)][b]

    # f =\frac{[1 - \exp(- \rho (|T^{n-1}_{i-1}| + s_a))] 
    #       [1 - \exp(- s'_{b-1} k_{b-1} / (2N))]}
    #      {\exp(-\rho |T^{n-1}_{i-1}|) (|T^{n-1}_{i-1}| + s_a) k_{b-1}}
    # |T^{n-1}_{i-1} = treelen

    transprob = util.make_matrix(len(states), len(states), 0.0)
    for i, (node1, a) in enumerate(states):
        c = times.index(tree[node1].age)
        for j, (node2, b) in enumerate(states):
            treelen2 = treelen + max(times[a], mintime)
            f = ((1.0 - exp(-rho * treelen2)) *
                 (1.0 - exp(-time_steps[b-1] * nlineages[b-1]
                            / (2.0 * popsizes[b-1]))) /
                 (exp(-rho * treelen) * treelen2 * nlineages[b-1]))
            if node1 != node2:
                transprob[i][j] = f * S[a][b]
            elif a != b:
                transprob[i][j] = f * (2*S[a][b] - S[c][b])
            else:
                # compute at the end
                pass

        transprob[i][i] = 1.0 - sum(transprob[i])
        for j in xrange(len(states)):
            transprob[i][j] = util.safelog(transprob[i][j])
        
    return transprob


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

    treelen = sum(x.get_dist() for x in last_tree)

    # find recomb node
    recomb_node = tree[recomb_name]
    k = times.index(recomb_node.age)

    # find re-coal point
    coal = recomb_node.parents[0]
    while coal.name not in last_tree and coal.parents:
        coal = coal.parents[0]
    coal_time = times.index(coal.age)

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
    
    # compute transition probability matrix
    transprob = util.make_matrix(len(states1), len(states2), 0.0)
    determ = get_deterministic_transitions(states1, states2, times,
                                           tree, last_tree,
                                           recomb_branch, k,
                                           coal_branch, coal_time)
    
    for i, (node1, a) in enumerate(states1):        
        if (node1, a) == (coal_branch, coal_time):
            # probabilistic transition case
            #print "complex", i
            
            # \frac{k^{(n-1)}_{j-1}}{k^{(n-1)}_{j-1} + 1}
            # \frac{[1 - \exp(- s'_{j-1} k^{(n)}_{j-1} / (2N))]}
            #      {[1 - \exp(- s'_{j-1} k^{(n-1)}_{j-1} / (2N))]}
            # \frac{|T^{n-1}_{i-1}|
            #       [1 - \exp(- \rho (|T^{n-1}_{i-1}| + t_{i-1}))]}
            #      {[|T^{n-1}_{i-1}| + t_{i-1}]
            #       [1 - \exp(- \rho |T^{n-1}_{i-1}|)]} 
            # \exp(- \sum_{m=k}^{j-2} s'_k / (2N))                 

            #transprob[i][i] = .1
            #continue
            
            if (node1, a) in states2:
                j = states2.index((node1, a))

                # TODO: fix
                # - 1 lineage if recomb br in time segment b-1
                b = a                
                if (recomb_branch, b) in states2:
                    kbn1 = max(nlineages[b-1] - 1, 1.0)
                else:
                    kbn1 = nlineages[b-1]
                kbn  = kbn1 + 1

                transprob[i][j] = (
                    (kbn1/float(kbn)) *
                    ((1.0 - exp(-time_steps[b-1] * kbn /
                                (2.0*popsizes[b-1])))/
                     (1.0 - exp(-time_steps[b-1] * kbn1 /
                                (2.0*popsizes[b-1]))))*
                    (treelen / (treelen + times[a])) *
                    ((1.0 - exp(-rho * (treelen + times[a]))) /
                     (1.0 - exp(-rho * treelen))) *
                    exp(- sum(time_steps[m] / (2.0 * popsizes[m])
                              for m in xrange(k, b-1))))

            for j, (node2, b) in enumerate(states2):
                if node2 != recomb_branch:
                    continue

                # TODO: fix
                # - 1 lineage if recomb br in time segment b-1
                if (recomb_branch, b) in states2:
                    kbn1 = max(nlineages[b-1] - 1, 1.0)
                else:
                    kbn1 = nlineages[b-1]
                kbn  = kbn1 + 1

                transprob[i][j] = (
                    (kbn1/float(kbn)) *
                    ((1.0 - exp(-time_steps[b-1] * kbn /
                                (2.0*popsizes[b-1])))/
                     (1.0 - exp(-time_steps[b-1] * kbn1 /
                                (2.0*popsizes[b-1]))))*
                    (treelen / (treelen + times[a])) *
                    ((1.0 - exp(-rho * (treelen + times[a]))) /
                     (1.0 - exp(-rho * treelen))) *
                    exp(- sum(time_steps[m] / (2.0 * popsizes[m])
                              for m in xrange(k, b-1))))
               

        else:
            # otherwise use deterministic transition
            #state1 = (node1, a)
            #j = get_deterministic_transition(state1, states2, times,
            #                                 tree, last_tree,
            #                                 recomb_branch, k,
            #                                 coal_branch, coal_time)
            j = determ[i]
            transprob[i][j] = 1.0

    for i in xrange(len(states1)):
        # HACK for now:  renormalize row to ensure they add up to one
        tot = sum(transprob[i])
        
        for j in xrange(len(states2)):
            transprob[i][j] = util.safelog(transprob[i][j] / tot)
        
    return transprob



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
    for state1 in states1:
        node1, a = state1
        leaves1 = set(last_tree.leaves(last_tree[node1]))
        remain = leaves1 - recomb_leaves
        
        if len(remain) > 0:
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
        



def get_deterministic_transition(state1, states2, times, tree, last_tree,
                                 recomb_branch, recomb_time,
                                 coal_branch, coal_time):

    # recomb_branch in tree
    # coal_branch in last_tree
    
    node1, a = state1

    # get leaves under node1
    leaves1 = set(last_tree.leaves(last_tree[node1]))

    # get leaves under recomb_node
    recomb_leaves = set(last_tree.leaves(last_tree[recomb_branch]))

    remain = leaves1 - recomb_leaves

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
    
    
    if len(remain) > 0:
        # SPR only removes a subset of descendents
        # trace up from remaining leaf to find correct new state
        ptr = tree[remain.pop().name]
        node = trace_up(ptr, a)
        return find_state(node, a)

    else:
        # SPR is on same branch as new chromosome
        if recomb_time >= a:
            # we move with SPR subtree
            ptr = tree[recomb_leaves.pop().name]
            node = trace_up(ptr, a)
            return find_state(node, a)

        elif coal_time <= a and coal_branch == node1:
            # SPR subtree coals back underneath us
            return find_state(tree[node1], a)
            
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
                return find_state(tree.root, b)
            ptr = tree[ptr.name]
            
            return find_state(ptr, b)
        

def find_tree_next_recomb(arg, pos, tree=False):

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
    """iterate through recombination visible in local trees"""
    
    pos = start if start is not None else 0
    while True:
        recomb = find_tree_next_recomb(arg, pos+1)
        if recomb:
            yield recomb
            pos = recomb.pos
        else:
            break



'''
def get_deterministic_transition2(state1, states2, times, tree, last_tree,
                                 recomb_branch, recomb_time,
                                 coal_branch, coal_time):

    # recomb_branch in tree
    # coal_branch in last_tree

    node1, a = state1

    if recomb_branch == node1:
        if recomb_time >= a:
            # state unaffected
            return states2.index(state1)
        else:
            # search up for parent
            ptr = last_tree[recomb_branch]
            last_ptr = ptr
            ptr = ptr.parents[0]
            while len(ptr.children) == 1:
                last_ptr = ptr
                ptr = ptr.parents[0]
            b = times.index(ptr.age)

            # go over to new tree
            if ptr.name not in tree:
                # we are above root
                assert ptr.age >= tree.root.age
                print tree.root.name, b
                # may have to go up coal point is not available
                #while (tree.root.name, b) not in states2:
                #    b += 1
                return states2.index((tree.root.name, b))
            ptr = tree[ptr.name]

            # go down to beginning of branch
            while len(ptr.children) == 1:
                ptr = ptr.children[0]

            # may have to go up coal point is not available
            while (ptr.name, b) not in states2:
                b += 1
            print ptr.name, b
            return states2.index((ptr.name, b))

    elif coal_branch == node1 and coal_time < a:
        # move to coal branch
        ptr = tree[coal_branch]
        last = ptr
        print last.name
        while ptr.age <= times[a] and ptr.parents:
            if len(ptr.children) != 1:
                last = ptr
            ptr = ptr.parents[0]
        print last.name, a
        # may have to go up coal point is not available
        #while (last.name, a) not in states2:
        #   a += 1
        return states2.index((last.name, a))

    elif coal_branch == node1 and coal_time == a:
        raise Exception("transition is not deterministic")

    else:
        # either unaffected or moved to next node with 0 or 2 children
        if node1 in tree:
            ptr = tree[node1]
        else:
            ptr = tree.root
        while len(ptr.children) == 1:
            ptr = ptr.children[0]
        return states2.index((ptr.name, a))

'''

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
local blocks   (0, 3) (3, 5)

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
                self.state_spaces.append(j)
            else:
                self.states.append(self.states[-1])
                self.state_spaces.append(j)

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
        return (self.recomb_pos[space-1]+1, self.recomb_pos[space]+1)
        


    def check_local_tree(self, pos, force=False):

        # update local block
        if force or not (self.local_block[0] <=  pos < self.local_block[1]):

            # get new local information
            self.local_tree = self.arg.get_marginal_tree(pos-.5)
            self.local_block = self.get_local_block(self.state_spaces[pos])
            self.nlineages = get_nlineages(self.local_tree, self.times)
            
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
                log((1 - exp(- self.time_steps[b-1] * self.nlineages[b-1] /
                         (2.0 * self.popsizes[b-1]))) / self.nlineages[b-1] *
                     exp(-sum(self.time_steps[m] * self.nlineages[m] /
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

        if node.parents:
            parent = node.parents[0]
            parent_age = parent.age

            if not parent.parents:
                # unwrap top branch
                c = parent.children
                sib = c[1] if node == c[0] else c[1]
                
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
            parent = None
            parent_age = None

            # adjust time by unwrapping branch
            time = 2 * time - node.age

            v = self.seqs[self.new_name][pos]
            x = self.local_site[node.name]
            p = x

        time = max(time, mintime)

        #print pos, v, x, p

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
                t2a = max(self.times[-1].age - time, mintime)
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



def forward_algorithm(model, n, verbose=False):

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

        util.tic("switch")
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
        util.toc()

        # do rest of block quickly
        space = model.get_state_space(i)
        block = model.get_local_block(space)
        blocklen = block[1] - i

        if i > block[0] and blocklen > 4:
            #util.tic("prep")
            nstates = model.get_num_states(i)

            #util.tic("emit")
            # setup tree and states
            tree = model.arg.get_marginal_tree(i-.5)
            tree2 = tree.get_tree()
            ptree, nodes, nodelookup = make_ptree(tree2)
            int_states = [[nodelookup[tree2[node]], timei]
                          for node, timei in model.states[i]]
            ages = [tree[node.name].age for node in nodes]
            seqs = [model.seqs[node.name][i-1:block[1]]
                    for node in nodes if node.is_leaf()]
            seqs.append(model.seqs[model.new_name][i-1:block[1]])
            seqlen = blocklen + 1
            
            #emit = c_matrix(
            #    c_double,
            #    [[model.prob_emission(pos, k) for k in xrange(nstates)]
            #     for pos in xrange(i-1, block[1])])
            
            emit = new_emissions(
                ((c_int * 2) * nstates)
                (* ((c_int * 2)(n, t) for n, t in int_states)), nstates, 
                ptree, len(ptree), ages,
                (c_char_p * len(seqs))(*seqs), len(seqs), seqlen,
                model.times, len(model.times), model.mu)
            #util.toc()

            #util.tic("trans")
            trans = c_matrix(
                c_double,
                [[model.prob_transition(i-1, j, i, k)
                  for k in xrange(nstates)] for j in xrange(nstates)])
            #util.toc()
            
            fw = [probs[-1]]
            for pos in xrange(i, block[1]):
                fw.append([0.0 for k in xrange(nstates)])
            #util.toc()
                
            #util.tic("fw %d" % blocklen)
            forward_alg(blocklen+1, nstates, nstates, fw, trans, emit)
            #util.toc()

            delete_emissions(emit, blocklen)

            for col in fw[1:]:
                probs.append(col[:nstates])
            nstates1 = nstates
            i = block[1]
            
    return probs




def iter_forward_algorithm(model, n, verbose=False):
    
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

        #space = model.get_space(i)
        #block = model.get_local_block(space)

        #blocklen = block[1] - i

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

'''
def forward_step(i, col1, col2, nstates1, nstates2, trans, emit):

    # find total transition and emission
    for k in xrange(nstates2):
        tot = -util.INF
        for j in xrange(nstates1):
            p = col1[j] + trans[j][k] + emit[k]
            tot = logadd(tot, p)
        col2[k] = tot
'''


def backward_algorithm(model, n, verbose=False):

    probs = []

    # calc last position
    nstates = model.get_num_states(n-1)
    for i in xrange(n):
        probs.append(None)
    probs[n-1] = [model.prob_prior(n-1, j) + model.prob_emission(n-1, j)
                  for j in xrange(nstates)]
    
    if n > 20:
        step = (n // 20)
    else:
        step = 1
    
    # loop through positions
    nstates1 = nstates
    i = n-2
    next_print = n-step
    while i > -1:
        if verbose and i < next_print:
            next_print -= step
            print " backward iter=%d/%d" % (i+1, n)

        # do first position manually
        nstates1 = model.get_num_states(i)
        nstates2 = model.get_num_states(i+1)
        col2 = probs[i+1]

        model.check_local_tree(i+1)
        if i+1 == model.local_block[0] and model.transmat_switch:
            trans = model.transmat_switch
        else:
            trans = model.transmat


        # find total transition and emission
        col1 = []
        emit = [model.prob_emission(i+1, k) for k in xrange(nstates2)]
        for j in xrange(nstates1):
            tot = -util.INF
            for k in xrange(nstates2):
                p = col2[k] + emit[k] + trans[j][k]
                tot = logadd(tot, p)
            col1.append(tot)
        probs[i] = col1
        i -= 1
        if i <= -1:
            break

        # do rest of block quickly
        space = model.get_state_space(i)
        block = model.get_local_block(space)
        blocklen = i+1 - block[0]
        if i < block[1]-1 and blocklen > 4:
            #print i, block, blocklen

            nstates = model.get_num_states(i)

            #util.tic("emit")
            # setup tree and states
            tree = model.arg.get_marginal_tree(i-.5)
            tree2 = tree.get_tree()
            ptree, nodes, nodelookup = make_ptree(tree2)
            int_states = [[nodelookup[tree2[node]], timei]
                          for node, timei in model.states[i]]
            ages = [tree[node.name].age for node in nodes]
            seqs = [model.seqs[node.name][block[0]:i+2]
                    for node in nodes if node.is_leaf()]
            seqs.append(model.seqs[model.new_name][block[0]:i+2])
            seqlen = blocklen + 1
            
            
            #emit = c_matrix(c_double,
            #    [[model.prob_emission(pos, k) for k in xrange(nstates)]
            #     for pos in xrange(block[0], i+2)])
            
            emit = new_emissions(
                ((c_int * 2) * nstates)
                (* ((c_int * 2)(n, t) for n, t in int_states)), nstates, 
                ptree, len(ptree), ages,
                (c_char_p * len(seqs))(*seqs), len(seqs), seqlen,
                model.times, len(model.times), model.mu)

            #util.toc()
            
            trans = c_matrix(c_double,
                             [[model.prob_transition(i, j, i+1, k)
                               for k in xrange(nstates)]
                              for j in xrange(nstates)])
            bw = [[0.0 for k in xrange(nstates)]
                  for pos in xrange(block[0], i+1)]
            bw.append(probs[i+1])
            backward_alg(blocklen+1, nstates, nstates, bw, trans, emit)
            for j in xrange(blocklen):
                probs[block[0]+j] = bw[j][:nstates]
            i = block[0] - 1

    return probs



def get_posterior_probs(model, n, verbose=False,
                        probs_forward=None, probs_backward=None):

    if probs_forward is None:
        probs_forward = forward_algorithm(model, n, verbose=verbose)
    if probs_backward is None:
        probs_backward = backward_algorithm(model, n, verbose=verbose)

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


def sample_posterior(model, n, probs_forward=None, verbose=False):

    path = range(n)

    # get forward probabilities
    if probs_forward is None:
        probs_forward = forward_algorithm(model, n, verbose=verbose)

    # base case i=n-1
    B = 0.0
    i = n-1
    A = [probs_forward[i][j] for j in range(model.get_num_states(i))]
    path[i] = j = stats.sample(map(exp, A))
  
    # recurse
    for i in xrange(n-2, -1, -1):
        C = []
        A = []
        for j in range(model.get_num_states(i)):
            # !$A_{j,i} = F_{i,j} C_{i,j} B_{i+1,l}$!
            C.append(
                model.prob_transition(i, j, i+1, path[i+1]) +
                model.prob_emission(i+1, path[i+1]))
            A.append(probs_forward[i][j] + C[j] + B)
        path[i] = j = stats.sample(map(exp, A))
        # !$B_{i,j} = C_{i,j} B_{i+1,l}$!
        B *= C[j]
    
    return path



#=============================================================================
# C interface functions

def make_ptree(tree, skip_single=True):
    """Make parent tree array from tree"""

    nodes = []
    nodelookup = {}
    ptree = []

    if skip_single:
        nodes = list(x for x in tree.postorder() if len(x.children) != 1)
    else:
        nodes = list(tree.postorder())
    
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
    
    assert nodes[-1] == tree.root
    
    return ptree, nodes, nodelookup


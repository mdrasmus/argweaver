#
# Ancestral Recombination Graph Hidden Markov Model (ArgHmm)
#

from math import exp, log
import random

from rasmus import hmm, util
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



#=============================================================================
# discretization


def get_time_point(i, ntimes, maxtime, delta=10):
    return (exp(i/float(ntimes) * log(1 + delta * maxtime)) - 1) / delta


def get_time_points(ntimes=30, maxtime=45000, delta=.01):
    return [get_time_point(i, ntimes, maxtime, delta)
            for i in range(ntimes+1)]


def iter_coal_states_tree(tree, times):
    ages = {}

    for node in tree.postorder():
        if node.is_leaf():
            ages[node] = 0.0
        
        i, j = util.binsearch(times, ages[node])
        ages[node.parent] = ages[node] + node.dist

        # do not offer coalescing at bottom of branch if branch is non-zero len
        if ages[node.parent] > ages[node]:
            i += 1
        
        if node.parent:
            while times[i] <= ages[node.parent]:
                yield (node, times[i])
                i += 1
        else:
            while i < len(times):
                yield (node, times[i])
                i += 1


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
            if parent.age > node.age:
                i += 1
            while i < len(times) and times[i] <= parent.age:
                yield (node.name, i)
                i += 1
        else:
            # do not coalesce at bottom of root branch
            i += 1
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
    return nlineages



def discretize_arg(arg, times):
    for node in arg:
        i, j = util.binsearch(times, node.age)
        if j is None: j = len(times) - 1
        node.age = times[j]

        if node.event == "recomb":
            node.pos = int(node.pos)


def get_default_state(last_tree, tree, state, times, states):

    node_name, timei = state
    time = times[timei]
    node2times = dict(states)
    states = set(states)

    if node_name in tree:
        node = tree[node_name]

        # walk up to see if someone coalesced beneath thread
        # TODO: this is a bug...
        node2 = node
        node3 = node
        while node2.age <= time and node2.parents:
            node2 = node2.parents[0]
            if (node2.name, timei) in states:
                node3 = node2
                break

        if node3 != node:
            #print "case 1"
            return (node3.name, timei)

        # walk down to find next coal node
        node2 = node
        while len(node2.children) == 1:
            node2 = node2.children[0]
        #print "case 2"
        while (node2.name, timei) not in states:
            timei += 1
        return (node2.name, timei)

    else:
        node = last_tree[node_name]

        # go up until we find a node that is still in tree
        while node.parents and (node.name not in tree or
                                len(tree[node.name].children) < 2):
            node = node.parents[0]
        #print "case 3"
        if not node.parents:
            node = tree.root
        timei = times.index(node.age)
        while (node.name, timei) not in states:
            timei += 1
            if timei > len(times):
                timei = times.index(last_tree.root.age)
                node.name = tree.root.name
        return (node.name, timei)


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

def calc_transition_terms(nlineages, times, time_steps, popsizes, rho):

    ntimes = len(time_steps)

    # A_{k,j} =& s'_{j-2} k_{j-2} / (2N) + \sum_{m=k}^{j-3} s'_m k_m / (2N) \\
    #         =& s'_{j-2} k_{j-2} / (2N) + A_{k,j-1}.
    
    A = [[0] * ntimes] * ntimes
    for k in xrange(ntimes):
        # A[k][k] = A[k][k+1] = 0
        for j in xrange(k+2, ntimes):
            l = j - 2
            A[k][j] = A[k][j-1] + time_steps[l] * nlineages[l] / (2.0 * popsizes[l])

    # B_{c,a} =& \sum_{k=0}^{c} \exp(- A_{k,a}) \\
    #         =& B_{c-1,a} + \exp(- A_{c,a}).

    B = [[0] * ntimes] * ntimes
    for a in xrange(ntimes):
        B[0][a] = time_steps[0] / times[max(a,1)] * exp(-A[0][a])
        for c in xrange(1, a):
            B[c][a] = B[c-1][a] + time_steps[c] / times[max(a,1)] * exp(-A[c][a])

    # S_{a,b} &= (1 - \exp(- \rho s_a)) [1 / k_{b-1}]
    #      [1/a] [ 1 - \exp(- s'_{b-1} k_{b-1} / (2N)) ] B_{min(a-1,b-1),a} \\
    
    S = [[0] * ntimes] * ntimes
    for a in xrange(1, ntimes):
        for b in xrange(1, ntimes):
            S[a][b] = log(((1 - exp(- rho * times[a])) / nlineages[b-1]
                       * (1 - exp(- time_steps[b-1] * nlineages[b-1]
                                    / (2.0 * popsizes[b-1])))
                       * B[min(a-1, b-1)][a]))

    return A, S



def calc_transition_probs(tree, states, treelen, nlineages, times,
                          time_steps, popsizes, rho):

    ntimes = len(time_steps)
    
    # A_{k,j} =& s'_{j-2} k_{j-2} / (2N) + \sum_{m=k}^{j-3} s'_m k_m / (2N) \\
    #         =& s'_{j-2} k_{j-2} / (2N) + A_{k,j-1}.
    
    A = [[0] * ntimes] * ntimes
    for k in xrange(ntimes):
        # A[k][k] = A[k][k+1] = 0
        for j in xrange(k+2, ntimes):
            l = j - 2
            A[k][j] = A[k][j-1] + time_steps[l] * nlineages[l] / (2.0 * popsizes[l])

    # B_{c,a} =& \sum_{k=0}^{c} \exp(- A_{k,a}) \\
    #         =& B_{c-1,a} + \exp(- A_{c,a}).

    B = [[0] * ntimes] * ntimes
    for a in xrange(ntimes):
        B[0][a] = time_steps[0] / times[max(a,1)] * exp(-A[0][a])
        for c in xrange(1, a):
            B[c][a] = B[c-1][a] + time_steps[c] / times[max(a,1)] * exp(-A[c][a])

    # S_{a,b} &= B_{min(a-1,b-1),a}
    
    S = [[0] * ntimes] * ntimes
    for a in xrange(1, ntimes):
        for b in xrange(1, ntimes):
            S[a][b] = B[min(a-1, b-1)][a]

    # f =\frac{[1 - \exp(- \rho (|T^{n-1}_{i-1}| + s_a))] 
    #       [1 - \exp(- s'_{b-1} k_{b-1} / (2N))]}
    #      {\exp(-\rho |T^{n-1}_{i-1}|) (|T^{n-1}_{i-1}| + s_a) k_{b-1}}
    # |T^{n-1}_{i-1} = treelen

    transprob = [[0.0] * len(states)] * len(states)
    for i, (node1, a) in enumerate(states):
        c = times.index(tree[node1].age)
        for j, (node2, b) in enumerate(states):
            treelen2 = treelen + times[a]
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
            transprob[i][j] = log(transprob[i][j])
        
    return transprob

    


def calc_transition_probs_switch(tree, last_tree, recomb_name,
                                 states1, states2,
                                 nlineages, times,
                                 time_steps, popsizes, rho):

    treelen = sum(x.dists[0] for x in last_tree if x.dists)

    # find recomb node
    recomb_node = tree[recomb_name]
    k = times.index(recomb_node.age)

    # find re-coal point
    coal = recomb_node.parents[0]
    while coal.name not in last_tree:
        coal = ptr.parents[0]
    coal_time = times.index(coal.age)

    # find coal branch in last_tree
    ptr = last_tree[coal.name]
    while len(ptr.children) == 1:
        ptr = ptr.children[0]
    coal_branch = ptr.name

    # find recomb branch in tree
    recomb = tree[recomb_name]
    while len(recomb.children) == 1:
        recomb = recomb.children[0]
    
    # compute transition probability matrix
    transprob = [[-util.INF] * len(states2)] * len(states1)
    for i, (node1, a) in enumerate(states1):        
        if (node1, a) == (coal_branch, coal_time):
            # probabilistic transition case
            
            # \frac{k^{(n-1)}_{j-1}}{k^{(n-1)}_{j-1} + 1}
            # \frac{[1 - \exp(- s'_{j-1} k^{(n)}_{j-1} / (2N))]}
            #      {[1 - \exp(- s'_{j-1} k^{(n-1)}_{j-1} / (2N))]}
            # \frac{|T^{n-1}_{i-1}|
            #       [1 - \exp(- \rho (|T^{n-1}_{i-1}| + t_{i-1}))]}
            #      {[|T^{n-1}_{i-1}| + t_{i-1}]
            #       [1 - \exp(- \rho |T^{n-1}_{i-1}|)]} 
            # \exp(- \sum_{m=k}^{j-2} s'_k / (2N))                 
            
            if (node1, a) in states2:
                j = states2.index((node1, a))

                # TODO: fix
                b = a
                kbn1 = nlineages[b-1] # - 1 if recomb br in time segment b-1
                kbn  = kbn1 + 1

                transprob[i][j] = (
                    (kb1/(kb1+1.0)) *
                    ((1.0 - exp(-time_steps[b-1] * kbn /
                                (2.0*popsizes[b-1])))/
                     (1.0 - exp(-time_steps[b-1] * kbn1 /
                                (2.0*popsizes[b-1]))))*
                    (treelen / (treelen + times[a])) *
                    ((1.0 - exp(-rho * (treelen + times[a]))) /
                     (1.0 - exp(-rho * treelen))) *
                    exp(- sum(time_steps[m] / (2.0 * popsizes[m])
                              for m in xrange(k, b-1))))

            for j, (node2, b) in enuemrate(states2):
                if node2 != recomb.name:
                    continue

                # TODO: fix
                b = a
                kbn1 = nlineages[b-1] # - 1 if recomb br in time segment b-1
                kbn  = kbn1 + 1

                transprob[i][j] = (
                    (kb1/(kb1+1.0)) *
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
            state1 = (node1, a)
            j = get_deterministic_transition(state1, states2, tree, last_tree,
                                             recomb_branch, recomb_time,
                                             coal_branch, coal_time)
            transprob[i][j] = 1.0

    for j in xrange(len(states1)):
        for j in xrange(len(states2)):
            transprob[i][j] = log(transprob[i][j])
        
    return transprob

    
def get_deterministic_transition(state1, states2, tree, last_tree,
                                 recomb_branch, recomb_time,
                                 coal_branch, coal_time):
    node1, a = state1

    if recomb_branch == node1:
        if recomb_time >= a:
            # state unaffected
            return states2.index(state1)
        else:
            # search up for parent
            ptr = tree[coal_branch]
            last = ptr
            while len(ptr.children) == 1:
                last = ptr
                ptr = ptr.parents[0]
            b = times.index[ptr.age]

            # go down to sibling
            c = ptr.children
            sib = c[0] if c[1] == last else c[1]

            return states2.index((sib.name, b))

    elif coal_branch == node1 and coal_time < a:
        # move to coal branch
        return states2.index((coal_branch, a))

    elif coal_branch == node1 and coal_time == a:
        raise Exception("transition is not deterministic")

    else:
        # state unaffected
        return states2.index(state1)
        
    


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
# simulate additional chromosome

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
        

        # determine nstates
        self.recomb_pos = [-1] + arglib.get_recomb_pos(arg)
        self.recomb_pos.append(arg.end)
        self.states = []
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
            else:
                self.states.append(self.states[-1])

        # current local tree
        self.local_block = [-1, self.recomb_pos[1]]
        self.local_blocki = 0
        self.local_tree = None
        self.default_pos = None
        self.default_states = None
        self.local_site = None
        self.last_tree = None
        self.last_pos = None
        self.A = None
        self.S = None
        self.check_local_tree(0, True)
        


    def check_local_tree(self, pos, force=False):

        # update local block
        if force or not (self.local_block[0] <  pos <= self.local_block[1]):
            self.last_tree = self.local_tree
            self.local_tree = self.arg.get_marginal_tree(pos-.5)

            if pos == self.local_block[1] + 1:
                self.local_blocki += 1
            elif pos == self.local_block[0]:
                self.local_blocki -= 1
            else:
                i, self.local_blocki = util.binsearch(self.recomb_pos, pos)
            
            self.local_block = [self.recomb_pos[self.local_blocki-1],
                                self.recomb_pos[self.local_blocki]]

            # begining of new block, setup default states
            #if (pos-1 == self.local_block[0]) and self.last_tree:
            #    self.default_states = [
            #        get_default_state(
            #            self.last_tree, self.local_tree, state1,
            #            self.times, self.states[pos])
            #        for j, state1 in enumerate(self.states[pos-1])]
            
            arglib.remove_single_lineages(self.local_tree)
            self.nlineages = get_nlineages(self.local_tree, self.times)

            self.A, self.S = calc_transition_terms(
                self.nlineages, self.times, self.time_steps,
                self.popsizes, self.rho)
            

        # update local site
        if force or pos != self.last_pos:
            self.local_site = parsimony_ancestral_seq(
                self.local_tree, self.seqs, pos)
            

        self.last_pos = pos

        
    def get_default_states(self, pos):
        
        if self.default_pos == pos:
            return self.default_states
        else:
            self.default_pos = pos
            prev_tree = self.arg.get_marginal_tree(pos-1-.5)
            local_tree = self.arg.get_marginal_tree(pos-.5)
            
            self.default_states = [
                get_default_state(
                prev_tree, local_tree, state1,
                self.times, self.states[pos])
                for j, state1 in enumerate(self.states[pos-1])]

            return self.default_states
        

    def get_num_states(self, pos):
        return len(self.states[pos])


    def prob_prior(self, pos, state):

        self.check_local_tree(pos)
        node, b = self.states[pos][state]

        # P(t_1=s_b, x_1 | T_1) =
        # [1 / k_{b-1}] [ 1 - \exp(- s'_{b-1} k_{b-1} / (2N)) ] \exp(- A_{0,b}).
        return log(((1 - exp(- self.time_steps[b-1] * self.nlineages[b-1] /
                         (2.0 * self.popsizes[b-1]))) / self.nlineages[b-1] *
                     exp(-self.A[0][b])))
        
    def prob_transition(self, pos1, state1, pos2, state2):

        self.check_local_tree(pos2)
        node1, a = self.states[pos1][state1]
        node2, b = self.states[pos2][state2]
        
        # calc default state, compare state2 to default
        if pos1 == self.local_block[0]:
            #default_state = self.default_states[state1]
            default_state = self.get_default_states(pos2)[state1]
            a = default_state[1]
            assert default_state in self.states[pos2], \
                   (default_state, self.states[pos2])
        else:
            default_state = (node1, a)
        
        if default_state == (node2, b):
            return log(exp(-self.rho * self.times[a]) + exp(self.S[a][a]))
        else:
            return self.S[a][b]


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
    for i in xrange(1, n):
        if verbose and i % step == 0:
            print " forward iter=%d/%d, lnl=%f" % (i+1, n, max(probs[i-1]))

        nstates2 = model.get_num_states(i)
        col1 = probs[i-1]

        # find total transition and emission
        col2 = []
        for k in xrange(nstates2):
            tot = -util.INF
            emit = model.prob_emission(i, k)
            for j in xrange(nstates1):
                p = col1[j] + model.prob_transition(i-1, j, i, k) + emit
                tot = logadd(tot, p)
            col2.append(tot)
                
        probs.append(col2)
        nstates1 = nstates2

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
    for i in xrange(1, n):
        if verbose and i % step == 0:
            print " forward iter=%d/%d, lnl=%f" % (i+1, n, max(col1))

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

'''
def forward_step(i, col1, col2, nstates1, nstates2, trans, emit):

    # find total transition and emission
    for k in xrange(nstates2):
        tot = -util.INF
        for j in xrange(nstates1):
            p = col1[j] + trans[k][j] + emit[k]
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
    for i in xrange(n-2, -1, -1):
        if verbose and i % step == 0:
            print " backward iter=%d/%d, lnl=%f" % (i+1, n, max(probs[i+1]))

        nstates1 = model.get_num_states(i)
        nstates2 = model.get_num_states(i+1)
        col2 = probs[i+1]

        # find total transition and emission
        col1 = []
        emit = [model.prob_emission(i+1, k) for k in xrange(nstates2)]
        for j in xrange(nstates1):
            tot = -util.INF
            for k in xrange(nstates2):
                p = col2[k] + emit[k] + model.prob_transition(i, j, i+1, k)
                tot = logadd(tot, p)
            col1.append(tot)
                
        probs[i] = col1

    return probs


def get_posterior_probs(model, n, verbose=False):

    probs_forward = forward_algorithm(model, n, verbose=verbose)
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



def sample_posterior(model, n, forward_probs=None, verbose=False):

    path = range(n)

    # get forward probabilities
    if forward_probs is None:
        forward_probs = forward_algorithm(model, n, verbose=verbose)

    # base case i=n-1
    B = 0.0
    i = n-1
    A = [forward_probs[i][j] for j in range(model.get_num_states(i))]
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
            A.append(forward_probs[i][j] + C[j] + B)
        path[i] = j = stats.sample(map(exp, A))
        # !$B_{i,j} = C_{i,j} B_{i+1,l}$!
        B *= C[j]
    
    return path



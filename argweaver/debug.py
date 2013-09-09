
from math import exp

import arghmm
from arghmm import emit
import arghmm.ctypes_export as C
import arghmmc

from compbio import arglib
from rasmus import hmm
from rasmus import stats
from rasmus import util


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
            self.times = arghmm.get_time_points(ntimes, maxtime, delta)
        else:
            self.times = times
            ntimes = len(self.times) - 1
        self.time_steps = [self.times[i] - self.times[i-1]
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
            x.pos for x in arghmm.iter_visible_recombs(arg))
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
                    list(arghmm.iter_coal_states(arg.get_marginal_tree(i-.5),
                                                 self.times)))
                self.state_spaces.append(j-1)
            else:
                self.states.append(self.states[-1])
                self.state_spaces.append(j-1)

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
        if force or not (self.local_block[0] <= pos < self.local_block[1]):

            # get new local information
            self.local_tree = self.arg.get_marginal_tree(pos-.5)
            self.local_block = self.get_local_block(self.state_spaces[pos])
            self.nlineages = arghmm.get_nlineages_recomb_coal(
                self.local_tree, self.times)

            # get new transition matrices
            self.transmat = arghmm.calc_transition_probs(
                self.local_tree, self.states[pos], self.nlineages,
                self.times, self.time_steps, self.popsizes, self.rho)

            assert len(self.transmat) == len(self.states[pos])
            assert len(self.transmat[0]) == len(self.states[pos])

            # get switch matrix for beginning of block
            start = self.local_block[0]
            recomb = arghmm.find_tree_next_recomb(self.arg, start - 1)
            if start > 0 and recomb is not None:
                last_tree = self.arg.get_marginal_tree(start-1-.5)
                self.transmat_switch = arghmm.calc_transition_probs_switch(
                    self.local_tree, last_tree, recomb.name,
                    self.states[start-1], self.states[start],
                    self.nlineages, self.times,
                    self.time_steps, self.popsizes, self.rho)

                assert len(self.transmat_switch) == len(self.states[start-1])
                assert len(self.transmat_switch[0]) == len(self.states[start])
            else:
                self.transmat_switch = None

            # get prior matrix if needed
            self.priormat = arghmm.calc_state_priors(
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


def iter_posterior_times(model, probs, perc=.5):

    times = model.times

    for pos, probcol in enumerate(probs):
        col = [0.0] * len(times)

        for j, p in enumerate(probcol):
            node, timei = model.states[pos][j]
            col[timei] += exp(p)

        # normalize column
        s = sum(col)
        col = [x / s for x in col]

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
        nlineages = arghmm.get_nlineages_recomb_coal(tree, model.times)
        nbranches, nrecombs, ncoals = nlineages
        times_lookup = dict((t, i) for i, t in enumerate(model.times))
        tree2 = tree.get_tree()
        ptree, nodes, nodelookup = arghmmc.make_ptree(tree2)
        int_states = [[nodelookup[tree2[node]], timei]
                      for node, timei in model.states[pos]]
        nstates = len(int_states)
        #ages = [tree[node.name].age for node in nodes]
        ages_index = [times_lookup[tree[node.name].age]
                      for node in nodes]
        #treelen = sum(x.dist for x in tree2)
        treelen = arghmm.get_treelen(tree, model.times)

        # get new transition matrices
        transmat = arghmmc.new_transition_probs(
            len(nodes), ptree, ages_index, treelen,
            ((C.c_int * 2) * nstates)
            (* ((C.c_int * 2)(n, t) for n, t in int_states)), nstates,
            len(model.time_steps), model.times, model.time_steps,
            nbranches, nrecombs, ncoals,
            model.popsizes, model.rho)

        # get switch matrix for beginning of block
        start = block[0]
        recomb = arghmm.find_tree_next_recomb(model.arg, start - 1)
        if start > 0 and recomb is not None:
            assert last_tree is not None

            # debug:
            states1 = list(arghmm.iter_coal_states(last_tree, model.times))
            if states1 != model.states[start-1]:
                print "start", start
                print "recomb_pos", model.recomb_pos
                print "states1", states1
                print "model.states-1", model.states[start-1]
                print "model.states", model.states[start]

                assert states1 == model.states[start-1]

            #transmat_switch = arghmm.calc_transition_probs_switch(
            #    tree, last_tree, recomb.name,
            #    model.states[start-1], model.states[start],
            #    last_nlineages, model.times,
            #    model.time_steps, model.popsizes, model.rho)

            transmat_switch = arghmmc.calc_transition_probs_switch_c(
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

        emit = arghmmc.new_emissions(
            ((C.c_int * 2) * nstates)
            (* ((C.c_int * 2)(n, t) for n, t in int_states)), nstates,
            ptree, len(ptree), ages_index,
            (C.c_char_p * len(seqs))(*seqs), len(seqs), len(seqs[0]),
            model.times, len(model.times), model.mu)

        last_tree = tree
        last_nlineages = nlineages

        yield block, nstates, transmat, transmat_switch, emit


def forward_algorithm(model, n, verbose=False, matrices=None,
                      prior=[], internal=False):

    probs = []

    if verbose:
        util.tic("forward")

    trees, names = arghmmc.arg2ctrees(model.arg, model.times)

    (ptrees, ages, sprs, blocks), all_nodes = arghmmc.get_treeset(
        model.arg, model.times)
    seqs = [model.seqs[node] for node in all_nodes[0]
            if model.arg[node].is_leaf()]
    seqs.append(model.seqs[model.new_name])
    seqlen = len(seqs[0])

    fw = arghmmc.arghmm_forward_alg(trees, model.times, len(model.times),
                                    model.popsizes, model.rho, model.mu,
                                    (C.c_char_p * len(seqs))(*seqs), len(seqs),
                                    seqlen, len(prior) > 0, prior, internal)

    # map states to python state space
    all_states = arghmmc.get_state_spaces(trees, len(model.times), internal)

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

    arghmmc.delete_state_spaces(all_states, len(ptrees))
    arghmmc.delete_forward_matrix(fw, seqlen)

    if verbose:
        util.toc()

    return probs


def backward_algorithm(model, n, verbose=False, matrices=None):

    probs = []

    if verbose:
        util.tic("backward")

    # get prior matrix
    local_tree = model.arg.get_marginal_tree(n-.5)
    nlineages = arghmm.get_nlineages_recomb_coal(local_tree, model.times)
    priors = arghmm.calc_state_priors(
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

        arghmmc.backward_alg(blocklen, nstates, transmat, emit, bw)

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
            arghmmc.delete_emissions(emit, blocklen)
            arghmmc.delete_transition_probs(transmat, nstates)

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
        arghmmc.delete_trans_emit_matrices(matrices)

    total_prob = -util.INF
    for j in xrange(model.get_num_states(0)):
        total_prob = stats.logadd(total_prob,
                                  model.prob_prior(0, j) +
                                  model.prob_emission(0, j) +
                                  probs_backward[0][j])

    probs_post = [
        [probs_forward[i][j] + probs_backward[i][j] - total_prob
         for j in xrange(model.get_num_states(i))]
        for i in xrange(n)]

    return probs_post

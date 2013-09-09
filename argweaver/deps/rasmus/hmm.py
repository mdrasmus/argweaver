"""

  Light weight generic HMM algorithms


  Methods need the following user-defined functions:

    get_num_states(pos)
    prob_prior(pos, state)
    prob_emission(pos, state)
    prob_transition(pos1, state1, pos2, state2)


    probs = [[j for j in get_num_states(i)] 
             for i in xrange(npositions)]


"""

import random
from math import log, exp

from rasmus import util, stats
from stats import logadd


class HMM (object):
    """
    Base class for defining an Hidden Markov Model (HMM)
    """
    
    def __init__(self):
        pass

    def set_callbacks(self, get_num_states=None,
                      prob_prior=None,
                      prob_emission=None,
                      prob_transition=None,
                      emit=None):
        if get_num_states:
            self.get_num_states = get_num_states
        if prob_prior:
            self.prob_prior = prob_prior
        if prob_emission:
            self.prob_emission = prob_emission
        if prob_transition:
            self.prob_transition = prob_transition
        if emit:
            self.emit = emit


    def get_num_states(self, pos):
        """Returns the number of states at position 'pos'"""
        return 0

    def prob_prior(self, pos, state):
        """Returns the prior probability of a state"""
        return 0.0

    def prob_emission(self, pos, state):
        """
        Returns the emission probability at position 'pos' with state 'state'
        """
        return 0.0

    def prob_transition(self, pos1, state1, pos2, state2):
        """
        Returns transition probability for transitioning between states
        'state1' and 'state2' between position 'pos1' and position 'pos2'
        """
        return 0.0

    def emit(self, pos, state):
        """
        Returns emission data given state 'state'
        """
        return None



def sample_hmm_first_state(model):
    state = 0
    nstates = model.get_num_states(0)
    p = model.prob_prior(0, state)
    pick = log(random.random())
    while pick > p and state < nstates:
        state += 1
        p = logadd(p, model.prob_prior(0, state))
    return state


def sample_hmm_next_state(model, pos, state):
    nstates = model.get_num_states(pos)
    state2 = 0
    p = model.prob_transition(pos-1, state, pos, state2)
    pick = log(random.random())
    while pick > p and state2 < nstates:
        state2 += 1
        p = logadd(p, model.prob_transition(pos-1, state, pos, state2))
    return state2
        


def sample_hmm_states(model):
    
    # sample first state
    pos = 0
    state = sample_hmm_first_state(model)
    yield state

    # sample next states
    pos = 1
    while True:
        state = sample_hmm_next_state(model, pos, state)
        yield state
        pos += 1


def sample_hmm_data(model, states=None):
    
    if states is None:
        states = sample_hmm_states(model)

    for i, state in enumerate(states):
        yield model.emit(i, state)



def viterbi(model, n, verbose=False):
    """
    Compute argmax_path P(path|data)
    """

    probs = []
    ptrs = []

    # calc first position
    nstates = model.get_num_states(0)
    probs.append([model.prob_prior(0, j) + model.prob_emission(0, j)
                  for j in xrange(nstates)])
    ptrs.append([-1] * nstates)
    
    if n > 20:
        step = (n // 20)
    else:
        step = 1
    
    # loop through positions
    for i in xrange(1, n):
        if verbose and i % step == 0:
            print " viterbi iter=%d/%d, lnl=%f" % (i+1, n, max(probs[-1]))

        nstates1 = model.get_num_states(i-1)
        nstates2 = model.get_num_states(i)
        col1 = probs[i-1]

        # find max transition and emission
        col2 = []
        col2_ptr = []
        for k in xrange(nstates2):
            top = -util.INF
            ptr = -1
            emit = model.prob_emission(i, k)
            for j in xrange(nstates1):
                p = col1[j] + model.prob_transition(i-1, j, i, k) + emit
                if p > top:
                    top = p
                    ptr = j
            col2.append(top)
            col2_ptr.append(ptr)
                
        probs.append(col2)
        ptrs.append(col2_ptr)

    # find max traceback
    j = util.argmax(probs[-1])
    traceback = [0] * n
    traceback[n-1] = j
    for i in xrange(n-1, 0, -1):
        j = ptrs[i][j]
        traceback[i-1] = j

    return traceback



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

        # find total transition and emission
        col2 = []
        for k in xrange(nstates2):
            tot = -util.INF
            emit = model.prob_emission(i, k)
            for j in xrange(nstates1):
                p = col1[j] + model.prob_transition(i-1, j, i, k) + emit
                tot = logadd(tot, p)
            col2.append(tot)

        yield col2
        col1 = col2
        nstates1 = nstates2



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
        B += C[j]
    
    return path

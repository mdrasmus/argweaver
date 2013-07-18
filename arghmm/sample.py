#
# Sampling methods for ArgHmm
#

# python libs
from ctypes import c_double, c_char_p, c_int
from math import exp
import random

# rasmus combio libs
from compbio import arglib
from rasmus import stats
from rasmus import util

# arghmm libs
import arghmm


# TODO: update this method
def sample_recombinations_thread(model, thread, use_times=True):
    """Samples new recombination for a thread"""

    r = 0

    # assumes that recomb_pos starts with -1 and ends with arg.end
    arg_recomb = model.recomb_pos
    time_lookup = util.list2lookup(model.times)
    #minlen = model.time_steps[0]

    tree = model.arg.get_marginal_tree(-.5)
    #treelen = None #arghmm.get_treelen(tree, model.times)
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
            #treelen = arghmm.get_treelen(tree, model.times, use_basal=False)
            #treelen_b = treelen + arghmm.get_basal_length(tree, model.times)
            nlineages = arghmm.get_nlineages_recomb_coal(tree, model.times)
            nbranches, nrecombs, ncoals = nlineages

            if transmat is not None:
                arghmm.delete_transition_probs(transmat, nstates)
            transmat = arghmm.calc_transition_probs_c(
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

        if state == last_state:
            if pos > next_recomb:
                # sample the next recomb pos
                last_treelen = arghmm.get_treelen(last_tree, model.times,
                    use_basal=False)
                last_treelen2 = arghmm.get_treelen_branch(
                    last_tree, model.times, last_node, last_time,
                    use_basal=False)

                rate = max(1.0 - exp(-model.rho * (last_treelen2 - last_treelen)
                                     - selftrans), model.rho)

                next_recomb = pos + int(random.expovariate(rate))
                #next_recomb = pos + int(random.expovariate(rate2))

            if pos < next_recomb:
                continue


        next_recomb = -1
        last_treelen2 = arghmm.get_treelen_branch(
            last_tree, model.times, last_node, last_time,
            use_basal=False)
        last_treelen2_b = last_treelen2 + arghmm.get_basal_length(
            last_tree, model.times, last_node, last_time)
        #last_treelen2 = arghmm.get_treelen_branch(
        #    last_tree, model.times, last_node, last_time)
        statei = model.states[pos].index((node, timei))
        selftrans = transmat[statei][statei]

        # there must be a recombination
        # either because state changed or we choose to recombine
        recombs = []
        if node == last_node:
            # y = node, k in Sr(node)
            node_timei = time_lookup[tree[node].age]
            for k in range(node_timei, min(timei, last_timei)+1):
                recombs.append((node, k))

        # y = v, k in [0, min(timei, last_timei)]
        for k in range(0, min(timei, last_timei)+1):
            recombs.append((new_node, k))


        C = arghmm.calc_C(model.time_steps, nbranches, model.popsizes)
        j = timei
        probs = []
        for recomb in recombs:
            k = recomb[1]
            probs.append((nbranches[k] + 1) * model.time_steps[k] /
                         (ncoals[j] * (nrecombs[k] + 1.0) * last_treelen2_b) *
                         (1.0 - exp(-model.time_steps[j] * nbranches[j] /
                                    (2.0 * model.popsizes[j]))) *
                         (1.0 - exp(-model.rho * last_treelen2)) *
                         exp(-C[j] + C[k]))
        recomb_node, recomb_time = recombs[stats.sample(probs)]

        if use_times:
            recomb_time = model.times[recomb_time]
        yield (pos, recomb_node, recomb_time)






def py_sample_arg(seqs, ntimes=20, rho=1.5e-8, mu=2.5e-8, popsize=1e4,
                  verbose=False, force=False):
    """
    Sample ARG for sequences
    """

    def add_chrom(arg, new_name):
        if verbose:
            util.tic("adding %s..." % new_name)

        model = arghmm.ArgHmm(arg, seqs, new_name=new_name,
                              times=times, rho=rho, mu=mu, popsize=popsize)

        if verbose:
            util.tic("sample thread")
        path = arghmm.sample_posterior(model, arg.end)
        if verbose:
            util.toc()

        if verbose:
            util.tic("sample recombs")
        thread = list(arghmm.iter_thread_from_path(model, path))
        recombs = list(sample_recombinations_thread(
            model, thread))
        if verbose:
            util.toc()

        if verbose:
            util.tic("add thread")
        arg = arghmm.add_arg_thread(arg, new_name, thread, recombs)
        if verbose:
            util.toc()

        if verbose:
            util.toc()
        return arg, thread

    names = seqs.keys()
    length = len(seqs[names[0]])
    arg = arghmm.make_trunk_arg(0, length, name=names[0])
    times = arghmm.get_time_points(ntimes=ntimes, maxtime=80000, delta=.01)

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


def py_resample_arg(arg, seqs, new_name,
                 ntimes=20, rho=1.5e-8, mu=2.5e-8, popsize=1e4,
                 times=None,
                 verbose=False):

    keep = [x for x in seqs.keys() if x != new_name]
    arglib.subarg_by_leaf_names(arg, keep)
    arg = arglib.smcify_arg(arg)
    if times is None:
        times = arghmm.get_time_points(ntimes=ntimes, maxtime=80000, delta=.01)

    if verbose:
        util.tic("adding %s..." % new_name)
    model = arghmm.ArgHmm(arg, seqs, new_name=new_name,
                              times=times, rho=rho, mu=mu, popsize=popsize)

    if verbose:
        util.tic("sample thread")
    path = arghmm.sample_posterior(model, arg.end)
    if verbose:
        util.toc()

    if verbose:
        util.tic("sample recombs")
    thread = list(arghmm.iter_thread_from_path(model, path))
    recombs = list(arghmm.sample_recombinations_thread(
        model, thread))
    if verbose:
        util.toc()

    if verbose:
        util.tic("add thread")
    arg = arghmm.add_arg_thread(arg, new_name, thread, recombs)
    if verbose:
        util.toc()

    if verbose:
        util.toc()

    return arg


def py_resample_arg_all(arg, seqs, ntimes=20,
                        rho=1.5e-8, mu=2.5e-8, popsize=1e4,
                        times=None, verbose=False):
    if verbose:
        util.tic("resample all chromosomes")
    for new_name in sorted(seqs.keys()):
        arg = arghmm.resample_arg(arg, seqs, new_name,
                                  ntimes=ntimes, rho=rho, mu=mu,
                                  popsize=popsize,
                                  times=times, verbose=verbose)
    if verbose:
        util.toc()
    return arg


'''
def py_resample_arg_max(arg, seqs, new_name,
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

'''




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
                tot = stats.logadd(tot, p)
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
            ptree, nodes, nodelookup = arghmm.make_ptree(tree2)
            int_states = [[nodelookup[tree2[node]], timei]
                          for node, timei in model.states[i]]
            ages = [model.times.index(tree[node.name].age) for node in nodes]
            seqs = [model.seqs[node.name][i-1:block[1]]
                    for node in nodes if node.is_leaf()]
            seqs.append(model.seqs[model.new_name][i-1:block[1]])
            seqlen = blocklen + 1

            emit = arghmm.new_emissions(
                ((c_int * 2) * nstates)
                (* ((c_int * 2)(n, t) for n, t in int_states)), nstates,
                ptree, len(ptree), ages,
                (c_char_p * len(seqs))(*seqs), len(seqs), seqlen,
                model.times, len(model.times), model.mu)

            trans = arghmm.c_matrix(
                c_double,
                [[model.prob_transition(i-1, j, i, k)
                  for k in xrange(nstates)] for j in xrange(nstates)])

            fw = [probs[-1]]
            for pos in xrange(i, block[1]):
                fw.append([0.0 for k in xrange(nstates)])

            arghmm.forward_alg(blocklen+1, nstates, trans, emit, fw)

            arghmm.delete_emissions(emit, blocklen)

            for col in fw[1:]:
                probs.append(col[:nstates])
            nstates1 = nstates
            i = block[1]

    return probs


def py_forward_algorithm2(model, n, verbose=False, matrices=None):
    """
    Python-based forward algorithm

    Uses matrix iteration strategy
    """
    probs = []

    if verbose:
        util.tic("forward")

    # get prior matrix
    local_tree = model.arg.get_marginal_tree(-.5)
    nlineages = arghmm.get_nlineages_recomb_coal(local_tree, model.times)
    priors = arghmm.calc_state_priors(
        local_tree, model.states[0], nlineages,
        model.times, model.time_steps, model.popsizes, model.rho)
    probs.append(priors)

    matrices_given = matrices is not None
    if matrices is None:
        matrices = arghmm.iter_trans_emit_matrices(model, n)

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

        arghmm.forward_alg(blocklen, nstates, transmat, emit, fw)

        if not matrices_given:
            arghmm.delete_emissions(emit, blocklen)
            arghmm.delete_transition_probs(transmat, nstates)

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


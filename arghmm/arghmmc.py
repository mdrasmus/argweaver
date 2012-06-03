#
# C interface for ArgHmm
#

# python imports
import time

# rasmus compbio libs
from rasmus import treelib, util
from compbio import arglib

# import arghmm C lib
import arghmm
from arghmm.ctypes_export import *
arghmmclib = load_library(["..", "lib"], "libarghmm.so")


#=============================================================================
# export c functions

ex = Exporter(globals())
export = ex.export


if arghmmclib:
    # replace python function with c
    
    export(arghmmclib, "forward_alg", c_int,
           [c_int, "n", c_int, "nstates",
            c_double_p_p, "trans", c_double_p_p, "emit",
            c_out(c_double_matrix), "fw"])

    export(arghmmclib, "backward_alg", c_int,
           [c_int, "n", c_int, "nstates",
            c_double_p_p, "trans", c_double_p_p, "emit",
            c_out(c_double_matrix), "bw"])

    export(arghmmclib, "sample_hmm_posterior", c_int,
           [c_int, "n", c_int, "nstates",
            c_double_p_p, "trans", c_double_p_p, "emit", 
            c_out(c_double_matrix), "fw", c_out(c_int_list), "path"])


    export(arghmmclib, "new_transition_probs", c_double_p_p,
           [c_int, "nnodes", c_int_list, "ptree",
            c_int_list, "ages_index", c_double, "treelen",
            POINTER(c_int * 2), "states", c_int, "nstates",
            c_int, "ntimes", c_double_list, "times",
            c_double_list, "time_steps",
            c_int_list, "nbranches", c_int_list, "nrecombs",
            c_int_list, "ncoals", 
            c_double_list, "popsizes", c_double, "rho"])

    export(arghmmclib, "new_transition_probs_switch", c_double_p_p,
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

    export(arghmmclib, "delete_transition_probs", c_int,
           [c_double_p_p, "transition_probs", c_int, "nstates"])

    export(arghmmclib, "new_emissions", c_double_p_p,
           [POINTER(c_int * 2), "states",
            c_int, "nstates", 
            c_int_list, "ptree", c_int, "nnodes", c_int_list, "ages",
            c_char_p_p, "seqs", c_int, "nseqs", c_int, "seqlen",
            c_double_list, "times", c_int, "ntimes",
            c_double, "mu"])

    export(arghmmclib, "delete_emissions", c_int,
           [c_double_p_p, "emit", c_int, "seqlen"])


    export(arghmmclib, "arghmm_forward_alg", c_double_p_p,
           [c_int_matrix, "ptrees", c_int_matrix, "ages",
            c_int_matrix, "sprs", c_int_list, "blocklens",
            c_int, "ntrees", c_int, "nnodes", 
            c_double_list, "times", c_int, "ntimes",
            c_double_list, "popsizes", c_double, "rho", c_double, "mu",
            c_char_p_p, "seqs", c_int, "nseqs", c_int, "seqlen"])

    export(arghmmclib, "delete_double_matrix", c_int,
           [c_double_p_p, "mat", c_int, "nrows"])

    export(arghmmclib, "arghmm_sample_posterior", POINTER(c_int *2),
           [c_int_matrix, "ptrees", c_int_matrix, "ages",
            c_int_matrix, "sprs", c_int_list, "blocklens",
            c_int, "ntrees", c_int, "nnodes", 
            c_double_list, "times", c_int, "ntimes",
            c_double_list, "popsizes", c_double, "rho", c_double, "mu",
            c_char_p_p, "seqs", c_int, "nseqs", c_int, "seqlen",
            POINTER(POINTER(c_int *2)), "path"])

    export(arghmmclib, "arghmm_sample_thread", c_void_p,
           [c_int_matrix, "ptrees", c_int_matrix, "ages",
            c_int_matrix, "sprs", c_int_list, "blocklens",
            c_int, "ntrees", c_int, "nnodes", 
            c_double_list, "times", c_int, "ntimes",
            c_double_list, "popsizes", c_double, "rho", c_double, "mu",
            c_char_p_p, "seqs", c_int, "nseqs", c_int, "seqlen"])

    export(arghmmclib, "arghmm_sample_arg_seq", c_void_p,
           [c_double_list, "times", c_int, "ntimes",
            c_double_list, "popsizes", c_double, "rho", c_double, "mu",
            c_char_p_p, "seqs", c_int, "nseqs", c_int, "seqlen"])

    export(arghmmclib, "arghmm_sample_arg_refine", c_void_p,
           [c_double_list, "times", c_int, "ntimes",
            c_double_list, "popsizes", c_double, "rho", c_double, "mu",
            c_char_p_p, "seqs", c_int, "nseqs", c_int, "seqlen",
            c_int, "niters"])

    export(arghmmclib, "arghmm_resample_arg", c_void_p,
           [c_int_matrix, "ptrees", c_int_matrix, "ages",
            c_int_matrix, "sprs", c_int_list, "blocklens",
            c_int, "ntrees", c_int, "nnodes", 
            c_double_list, "times", c_int, "ntimes",
            c_double_list, "popsizes", c_double, "rho", c_double, "mu",
            c_char_p_p, "seqs", c_int, "nseqs", c_int, "seqlen",
            c_int, "niters", c_int, "nremove"])

    export(arghmmclib, "arghmm_remax_arg", c_void_p,
           [c_int_matrix, "ptrees", c_int_matrix, "ages",
            c_int_matrix, "sprs", c_int_list, "blocklens",
            c_int, "ntrees", c_int, "nnodes", 
            c_double_list, "times", c_int, "ntimes",
            c_double_list, "popsizes", c_double, "rho", c_double, "mu",
            c_char_p_p, "seqs", c_int, "nseqs", c_int, "seqlen",
            c_int, "niters", c_int, "nremove"])

    export(arghmmclib, "arghmm_likelihood", c_double,
           [c_int_matrix, "ptrees", c_int_matrix, "ages",
            c_int_matrix, "sprs", c_int_list, "blocklens",
            c_int, "ntrees", c_int, "nnodes", 
            c_double_list, "times", c_int, "ntimes",
            c_double, "mu",
            c_char_p_p, "seqs", c_int, "nseqs", c_int, "seqlen"])

    export(arghmmclib, "arghmm_joint_prob", c_double,
           [c_int_matrix, "ptrees", c_int_matrix, "ages",
            c_int_matrix, "sprs", c_int_list, "blocklens",
            c_int, "ntrees", c_int, "nnodes", 
            c_double_list, "times", c_int, "ntimes", c_double_list, "popsizes",
            c_double, "mu", c_double, "rho",
            c_char_p_p, "seqs", c_int, "nseqs", c_int, "seqlen"])

    export(arghmmclib, "arghmm_prior_prob", c_double,
           [c_int_matrix, "ptrees", c_int_matrix, "ages",
            c_int_matrix, "sprs", c_int_list, "blocklens",
            c_int, "ntrees", c_int, "nnodes", 
            c_double_list, "times", c_int, "ntimes", c_double_list, "popsizes",
            c_double, "rho"])

    export(arghmmclib, "get_local_trees_ntrees", c_int,
           [c_void_p, "trees"])
    export(arghmmclib, "get_local_trees_nnodes", c_int,
           [c_void_p, "trees"])
    export(arghmmclib, "get_local_trees_ptrees", c_int,
           [c_void_p, "trees", c_out(c_int_matrix), "ptrees",
            c_out(c_int_matrix), "ages",
            c_out(c_int_matrix), "sprs", c_out(c_int_list), "blocklens"])
    export(arghmmclib, "delete_local_trees", c_int,
           [c_void_p, "trees"])

    export(arghmmclib, "delete_path", c_int,
           [POINTER(c_int * 2), "path"])


    export(arghmmclib, "get_state_spaces", POINTER(POINTER(c_int * 2)),
           [c_int_matrix, "ptrees", c_int_matrix, "ages",
            c_int_matrix, "sprs", c_int_list, "blocklens",
            c_int, "ntrees", c_int, "nnodes", c_int, "ntimes"])

    export(arghmmclib, "delete_state_spaces", c_int,
           [POINTER(POINTER(c_int * 2)), "all_states", c_int, "ntrees"])


# by default use a random seed
arghmmclib.srand(int(time.time()*1000))

def set_random_seed(num):
    """Set the C random number generator seed"""
    arghmmclib.srand(num)
    

#=============================================================================
# HMM transition and emission matrices


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
    treelen = arghmm.get_treelen(tree, times)
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
        arghmm.find_recomb_coal(tree, last_tree, recomb_name=recomb_name)
    
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



def delete_trans_emit_matrices(matrices):
    """Delete matrices"""
    for block, nstates, transmat, transmat_switch, emit in matrices:
        delete_emissions(emit, block[1] - block[0])
        delete_transition_probs(transmat, nstates)
        


#=============================================================================
# sampling functions


def sample_thread(model, n, probs_forward=None, verbose=False, matrices=None):

    if verbose:
        util.tic("sample thread")

    (ptrees, ages, sprs, blocks), all_nodes = get_treeset(
        model.arg, model.times)
    blocklens = [x[1] - x[0] for x in blocks]
    seqlen = sum(blocklens)

    seqs = [model.seqs[node] for node in all_nodes[0]
            if model.arg[node].is_leaf()]
    seqs.append(model.seqs[model.new_name])
    
    trees = arghmm_sample_thread(
        ptrees, ages, sprs, blocklens,
        len(ptrees), len(ptrees[0]), 
        model.times, len(model.times),
        model.popsizes, model.rho, model.mu,
        (c_char_p * len(seqs))(*seqs), len(seqs),
        seqlen, None)

    nnodes = get_local_trees_nnodes(trees)
    ntrees = get_local_trees_ntrees(trees)
    
    ptrees = []
    ages = []
    sprs = []
    blocklens = [0] * ntrees
    for i in range(ntrees):
        ptrees.append([0] * nnodes)
        ages.append([0] * nnodes)
        sprs.append([0, 0, 0, 0])

    get_local_trees_ptrees(trees, ptrees, ages, sprs, blocklens)
    delete_local_trees(trees)

    for i in range(ntrees):
        ptrees[i] = ptrees[i][:nnodes]
        ages[i] = ages[i][:nnodes]
        sprs[i] = sprs[i][:4]
    
    
    blocks = []
    start = 0
    for blocklen in blocklens:
        end = start + blocklen
        blocks.append((start, end))
        start = end

    names = all_nodes[0][:len(seqs)-1]
    names.append(model.new_name)
    assert len(names) == ((nnodes + 1) / 2)

    util.tic("convert arg")
    arg = treeset2arg(ptrees, ages, sprs, blocks, names, model.times)
    util.toc()
    
    #delete_local_trees(trees)

    if verbose:
        util.toc()

    return arg



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



def sample_arg(seqs, ntimes=20, rho=1.5e-8, mu=2.5e-8, popsize=1e4,
               refine=0, times=None, verbose=False):
    """
    Sample ARG for sequences
    """
    if times is None:
        times = arghmm.get_time_points(ntimes=ntimes, maxtime=80000, delta=.01)
    popsizes = [popsize] * len(times)

    if verbose:
        util.tic("sample arg")

    seqs2 = seqs.values()
    
    trees = arghmm_sample_arg_refine(
        times, len(times),
        popsizes, rho, mu,
        (c_char_p * len(seqs))(*seqs2), len(seqs), len(seqs2[0]), refine)


    if verbose:
        util.tic("convert arg")

    nnodes = get_local_trees_nnodes(trees)
    ntrees = get_local_trees_ntrees(trees)
    
    ptrees = []
    ages = []
    sprs = []
    blocklens = [0] * ntrees
    for i in range(ntrees):
        ptrees.append([0] * nnodes)
        ages.append([0] * nnodes)
        sprs.append([0, 0, 0, 0])

    get_local_trees_ptrees(trees, ptrees, ages, sprs, blocklens)
    delete_local_trees(trees)

    for i in range(ntrees):
        ptrees[i] = ptrees[i][:nnodes]
        ages[i] = ages[i][:nnodes]
        sprs[i] = sprs[i][:4]
    
    
    blocks = []
    start = 0
    for blocklen in blocklens:
        end = start + blocklen
        blocks.append((start, end))
        start = end

    names = seqs.keys()
    assert len(names) == ((nnodes + 1) / 2)

    arg = treeset2arg(ptrees, ages, sprs, blocks, names, times)
    if verbose:
        util.toc()
        util.toc()
    
    return arg


def resample_arg(arg, seqs, ntimes=20, rho=1.5e-8, mu=2.5e-8, popsize=1e4,
                 refine=1, nremove=1, times=None, verbose=False):
    """
    Sample ARG for sequences
    """
    if times is None:
        times = arghmm.get_time_points(ntimes=ntimes, maxtime=80000, delta=.01)
    popsizes = [popsize] * len(times)

    if verbose:
        util.tic("resample arg")
        util.tic("convert arg")
    seqs2 = seqs.values()

    (ptrees, ages, sprs, blocks), all_nodes = get_treeset(
        arg, times)
    blocklens = [x[1] - x[0] for x in blocks]
    seqlen = sum(blocklens)
    
    names = []
    seqs2 = []
    for node in all_nodes[0]:
        if arg[node].is_leaf():
            names.append(node)
            seqs2.append(seqs[node])

    # add all other sequences not in arg yet
    leaves = set(arg.leaf_names())
    for name, seq in seqs.items():
        if name not in leaves:
            names.append(name)
            seqs2.append(seq)
            
    
    if verbose:
        util.toc()
    
    trees = arghmm_resample_arg(
        ptrees, ages, sprs, blocklens,
        len(ptrees), len(ptrees[0]), 
        times, len(times),
        popsizes, rho, mu,
        (c_char_p * len(seqs2))(*seqs2), len(seqs2),
        seqlen, refine, nremove)

    if verbose:
        util.tic("convert arg 2")
    
    nnodes = get_local_trees_nnodes(trees)
    ntrees = get_local_trees_ntrees(trees)
    
    ptrees = []
    ages = []
    sprs = []
    blocklens = [0] * ntrees
    for i in range(ntrees):
        ptrees.append([0] * nnodes)
        ages.append([0] * nnodes)
        sprs.append([0, 0, 0, 0])

    get_local_trees_ptrees(trees, ptrees, ages, sprs, blocklens)
    delete_local_trees(trees)

    for i in range(ntrees):
        ptrees[i] = ptrees[i][:nnodes]
        ages[i] = ages[i][:nnodes]
        sprs[i] = sprs[i][:4]
    
    
    blocks = []
    start = 0
    for blocklen in blocklens:
        end = start + blocklen
        blocks.append((start, end))
        start = end

    assert len(names) == ((nnodes + 1) / 2)

    arg = treeset2arg(ptrees, ages, sprs, blocks, names, times)
    if verbose:
        util.toc()
        util.toc()
    
    return arg


def remax_arg(arg, seqs, ntimes=20, rho=1.5e-8, mu=2.5e-8, popsize=1e4,
                 refine=1, nremove=1, times=None, verbose=False):
    """
    Sample ARG for sequences
    """
    if times is None:
        times = arghmm.get_time_points(ntimes=ntimes, maxtime=80000, delta=.01)
    popsizes = [popsize] * len(times)

    if verbose:
        util.tic("remax arg")
        util.tic("convert arg")
    seqs2 = seqs.values()

    (ptrees, ages, sprs, blocks), all_nodes = get_treeset(
        arg, times)
    blocklens = [x[1] - x[0] for x in blocks]
    seqlen = sum(blocklens)
    
    names = []
    seqs2 = []
    for node in all_nodes[0]:
        if arg[node].is_leaf():
            names.append(node)
            seqs2.append(seqs[node])

    # add all other sequences not in arg yet
    leaves = set(arg.leaf_names())
    for name, seq in seqs.items():
        if name not in leaves:
            names.append(name)
            seqs2.append(seq)
            
    
    if verbose:
        util.toc()
    
    trees = arghmm_remax_arg(
        ptrees, ages, sprs, blocklens,
        len(ptrees), len(ptrees[0]), 
        times, len(times),
        popsizes, rho, mu,
        (c_char_p * len(seqs2))(*seqs2), len(seqs2),
        seqlen, refine, nremove)

    if verbose:
        util.tic("convert arg 2")
    
    nnodes = get_local_trees_nnodes(trees)
    ntrees = get_local_trees_ntrees(trees)
    
    ptrees = []
    ages = []
    sprs = []
    blocklens = [0] * ntrees
    for i in range(ntrees):
        ptrees.append([0] * nnodes)
        ages.append([0] * nnodes)
        sprs.append([0, 0, 0, 0])

    get_local_trees_ptrees(trees, ptrees, ages, sprs, blocklens)
    delete_local_trees(trees)

    for i in range(ntrees):
        ptrees[i] = ptrees[i][:nnodes]
        ages[i] = ages[i][:nnodes]
        sprs[i] = sprs[i][:4]
    
    
    blocks = []
    start = 0
    for blocklen in blocklens:
        end = start + blocklen
        blocks.append((start, end))
        start = end

    assert len(names) == ((nnodes + 1) / 2)

    arg = treeset2arg(ptrees, ages, sprs, blocks, names, times)
    if verbose:
        util.toc()
        util.toc()
    
    return arg




def calc_likelihood(arg, seqs, ntimes=20, mu=2.5e-8, 
                    times=None, verbose=False):
    """
    Calculate arg_likelihood
    """
    if times is None:
        times = arghmm.get_time_points(ntimes=ntimes, maxtime=80000, delta=.01)

    if verbose:
        util.tic("calc likelihood")
        util.tic("convert arg")
    seqs2 = seqs.values()

    (ptrees, ages, sprs, blocks), all_nodes = get_treeset(
        arg, times)
    blocklens = [x[1] - x[0] for x in blocks]
    seqlen = sum(blocklens)
    
    #names = []
    seqs2 = []
    for node in all_nodes[0]:
        if arg[node].is_leaf():
            #names.append(node)
            seqs2.append(seqs[node])

    # add all other sequences not in arg yet
    #leaves = set(arg.leaf_names())
    #for name, seq in seqs.items():
    #    if name not in leaves:
    #        names.append(name)
    #        seqs2.append(seq)
            
    
    if verbose:
        util.toc()
    
    lk = arghmm_likelihood(
        ptrees, ages, sprs, blocklens, len(ptrees), len(ptrees[0]), 
        times, len(times), mu,
        (c_char_p * len(seqs2))(*seqs2), len(seqs2), seqlen)

    if verbose:
        util.toc()
    
    return lk



def calc_joint_prob(arg, seqs, ntimes=20, mu=2.5e-8, rho=1.5e-8, popsize=1e4,
                    times=None, verbose=False):
    """
    Calculate arg_joint_prob
    """
    if times is None:
        times = arghmm.get_time_points(ntimes=ntimes, maxtime=80000, delta=.01)
    popsizes = [popsize] * len(times)

    if verbose:
        util.tic("calc likelihood")
        util.tic("convert arg")

    (ptrees, ages, sprs, blocks), all_nodes = get_treeset(
        arg, times)
    blocklens = [x[1] - x[0] for x in blocks]
    seqlen = sum(blocklens)
    
    #names = []
    seqs2 = []
    for node in all_nodes[0]:
        if arg[node].is_leaf():
            #names.append(node)
            seqs2.append(seqs[node])
    
    # add all other sequences not in arg yet
    #leaves = set(arg.leaf_names())
    #for name, seq in seqs.items():
    #    if name not in leaves:
    #        #names.append(name)
    #        seqs2.append(seq)

    #print [arg[node] for node in all_nodes[0] if arg[node].is_leaf()]
    
    if verbose:
        util.toc()
    
    lk = arghmm_joint_prob(
        ptrees, ages, sprs, blocklens, len(ptrees), len(ptrees[0]), 
        times, len(times), popsizes, mu, rho,
        (c_char_p * len(seqs2))(*seqs2), len(seqs2), seqlen)

    if verbose:
        util.toc()
    
    return lk


def calc_prior_prob(arg, ntimes=20, rho=1.5e-8, popsize=1e4,
                    times=None, verbose=False):
    """
    Calculate arg_joint_prob
    """
    if times is None:
        times = arghmm.get_time_points(ntimes=ntimes, maxtime=80000, delta=.01)
    popsizes = [popsize] * len(times)

    if verbose:
        util.tic("calc likelihood")
        util.tic("convert arg")

    (ptrees, ages, sprs, blocks), all_nodes = get_treeset(
        arg, times)
    blocklens = [x[1] - x[0] for x in blocks]
    
    if verbose:
        util.toc()
    
    lk = arghmm_prior_prob(
        ptrees, ages, sprs, blocklens, len(ptrees), len(ptrees[0]), 
        times, len(times), popsizes, rho)

    if verbose:
        util.toc()
    
    return lk



#=============================================================================
# tree functions


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
# passing ARG through C interface


def iter_arg_sprs(arg, start=None, end=None):
    """
    Iterates through the SPRs of an ARG
    """

    if start is None:
        start = arg.start
    if end is None:
        end = arg.end

    last_tree_full = None
    last_tree = None
    for block, tree_full in arglib.iter_tree_tracks(arg, start, end):
        if last_tree_full:
            recomb = (x for x in tree_full if x.pos == block[0]).next()
            spr = arghmm.find_recomb_coal(tree_full, last_tree_full,
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

        tree2 = tree.get_tree()
        #phylo.hash_order_tree(tree2)
        #print block[0], tree2.get_one_line_newick()
        
        if last_tree is None:
            # get frist ptree
            ptree, nodes, nodelookup = make_ptree(tree2)
            ispr = [-1, -1, -1, -1]
            age = [times_lookup[tree[x.name].age] for x in nodes]

        else:
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

            #print last_tree2[rname].leaf_names()

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


def treeset2arg(ptrees, ages, sprs, blocks, names, times):

    seqlen = blocks[-1][1]
    arg = arglib.ARG(0, seqlen)

    
    # build first tree
    lookup = {}
    for i, p in enumerate(ptrees[0]):
        if i < len(names):
            # make leaf
            lookup[i] = arg.new_node(names[i], age=times[ages[0][i]],
                                     event="gene")
        else:
            lookup[i] = arg.new_node(age=times[ages[0][i]], event="coal")
    
    # set parents of new tree
    for i, p in enumerate(ptrees[0]):
        node = lookup[i]
        if p != -1:
            node.parents.append(lookup[p])
            node.parents[0].children.append(node)
        else:
            arg.root = node

    
    # convert sprs
    sprs2 = []
    for i, (rinode, ritime, cinode, citime) in enumerate(sprs):
        pos = blocks[i][0] - 1

        # check for null spr
        if rinode == -1:
            continue

        # make local tree
        ptree = ptrees[i-1]
        tree = treelib.Tree()
        lookup = []
        for j in range(len(ptree)):
            if j < len(names):
                lookup.append(tree.new_node(names[j]))
            else:
                lookup.append(tree.new_node())
        for j in range(len(ptree)):
            if ptree[j] != -1:
                parent = lookup[ptree[j]]
                tree.add_child(parent, lookup[j])
            else:
                tree.root = lookup[j]
        
        #phylo.hash_order_tree(tree)
        #print pos+1, tree.get_one_line_newick()
        
        # get leaf sets
        rleaves = lookup[rinode].leaf_names()
        cleaves = lookup[cinode].leaf_names()
        assert ritime >= ages[i-1][rinode], (pos, ritime, ages[i-1][rinode])
        assert citime >= ages[i-1][cinode], (pos, citime, ages[i-1][cinode])
        
        sprs2.append((pos, (rleaves, times[ritime]), (cleaves, times[citime])))

    #assert against local ptree and leading edge of ARG.

    arglib.make_arg_from_sprs(arg, sprs2)
    arglib.assert_arg(arg)
    return arg
    



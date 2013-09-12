#
# C interface for ARGweaver
#

# python imports
import time

# rasmus compbio libs
from compbio import arglib
from rasmus import stats
from rasmus import treelib
from rasmus import util

# import argweaver C lib
import argweaver
import argweaver.ctypes_export as C
from argweaver.ctypes_export import POINTER

argweaverclib = C.load_library(["..", "lib"], "libargweaver.so")


#=============================================================================
# export c functions

ex = C.Exporter(globals())
export = ex.export


if argweaverclib:
    # replace python function with c

    # Sequences.
    argweaver_read_sites = export(
        argweaverclib, "arghmm_read_sites", C.c_void_p,
        [C.c_char_p, "filename", C.c_int, "subregion_start",
         C.c_int, "subregion_end"])
    argweaver_delete_sites = export(
        argweaverclib, "arghmm_delete_sites", C.c_int,
        [C.c_void_p, "sites"])

    # Basic HMM functions.
    forward_alg = export(
        argweaverclib, "forward_alg", C.c_int,
        [C.c_int, "n", C.c_int, "nstates",
         C.c_double_p_p, "trans", C.c_double_p_p, "emit",
         C.c_out(C.c_double_matrix), "fw"])
    backward_alg = export(
        argweaverclib, "backward_alg", C.c_int,
        [C.c_int, "n", C.c_int, "nstates",
         C.c_double_p_p, "trans", C.c_double_p_p, "emit",
         C.c_out(C.c_double_matrix), "bw"])
    sample_hmm_posterior = export(
        argweaverclib, "sample_hmm_posterior", C.c_int,
        [C.c_int, "n", C.c_int, "nstates",
         C.c_double_p_p, "trans", C.c_double_p_p, "emit",
         C.c_out(C.c_double_matrix), "fw", C.c_out(C.c_int_list), "path"])

    # Transition matrices calculation.
    new_transition_probs = export(
        argweaverclib, "new_transition_probs", C.c_double_p_p,
        [C.c_int, "nnodes", C.c_int_list, "ptree",
         C.c_int_list, "ages_index", C.c_double, "treelen",
         POINTER(C.c_int * 2), "states", C.c_int, "nstates",
         C.c_int, "ntimes", C.c_double_list, "times",
         C.c_double_list, "time_steps",
         C.c_int_list, "nbranches", C.c_int_list, "nrecombs",
         C.c_int_list, "ncoals",
         C.c_double_list, "popsizes", C.c_double, "rho"])
    new_transition_probs_switch = export(
        argweaverclib, "new_transition_probs_switch", C.c_double_p_p,
        [C.c_int_list, "ptree", C.c_int_list, "last_ptree",
         C.c_int, "nnodes",
         C.c_int, "recomb_name", C.c_int, "recomb_time",
         C.c_int, "coal_name", C.c_int, "coal_time",
         C.c_int_list, "ages_index", C.c_int_list, "last_ages_index",
         C.c_double, "treelen", C.c_double, "last_treelen",
         POINTER(C.c_int * 2), "states1", C.c_int, "nstates1",
         POINTER(C.c_int * 2), "states2", C.c_int, "nstates2",
         C.c_int, "ntimes", C.c_double_list, "times",
         C.c_double_list, "time_steps",
         C.c_int_list, "nbranches", C.c_int_list, "nrecombs",
         C.c_int_list, "ncoals",
         C.c_double_list, "popsizes", C.c_double, "rho"])
    delete_transition_probs = export(
        argweaverclib, "delete_transition_probs", C.c_int,
        [C.c_double_p_p, "transition_probs", C.c_int, "nstates"])
    argweaver_assert_transmat = export(
        argweaverclib, "arghmm_assert_transmat", C.c_bool,
        [C.c_int, "nnodes", C.c_int_list, "ptree", C.c_int_list, "ages",
         C.c_int, "ntimes", C.c_double_list, "times",
         C.c_double_list, "popsizes", C.c_double, "rho"])
    argweaver_assert_transmat_switch = export(
        argweaverclib, "arghmm_assert_transmat_switch", C.c_bool,
        [C.c_int, "nnodes", C.c_int_list, "ptree", C.c_int_list, "ages",
         C.c_int, "recomb_name", C.c_int, "recomb_time",
         C.c_int, "coal_name", C.c_int, "coal_time",
         C.c_int, "ntimes", C.c_double_list, "times",
         C.c_double_list, "popsizes", C.c_double, "rho"])
    argweaver_assert_transmat_internal = export(
        argweaverclib, "arghmm_assert_transmat_internal", C.c_bool,
        [C.c_int, "nnodes", C.c_int_list, "ptree", C.c_int_list, "ages",
         C.c_int, "ntimes", C.c_double_list, "times",
         C.c_double_list, "popsizes", C.c_double, "rho"])
    argweaver_assert_transmat_switch_internal = export(
        argweaverclib, "arghmm_assert_transmat_switch_internal", C.c_bool,
        [C.c_void_p, "trees", C.c_int, "ntimes", C.c_double_list, "times",
         C.c_double_list, "popsizes", C.c_double, "rho"])

    # Emission calculation.
    new_emissions = export(
        argweaverclib, "new_emissions", C.c_double_p_p,
        [POINTER(C.c_int * 2), "states",
         C.c_int, "nstates",
         C.c_int_list, "ptree", C.c_int, "nnodes", C.c_int_list, "ages",
         C.c_char_p_p, "seqs", C.c_int, "nseqs", C.c_int, "seqlen",
         C.c_double_list, "times", C.c_int, "ntimes",
         C.c_double, "mu"])
    delete_emissions = export(
        argweaverclib, "delete_emissions", C.c_int,
        [C.c_double_p_p, "emit", C.c_int, "seqlen"])
    argweaver_assert_emit = export(
        argweaverclib, "arghmm_assert_emit", C.c_bool,
        [C.c_void_p, "trees", C.c_int, "ntimes", C.c_double_list, "times",
         C.c_double, "mu",
         C.c_char_p_p, "seqs", C.c_int, "nseqs", C.c_int, "seqlen"])
    argweaver_assert_emit_internal = export(
        argweaverclib, "arghmm_assert_emit_internal", C.c_bool,
        [C.c_void_p, "trees", C.c_int, "ntimes", C.c_double_list, "times",
         C.c_double, "mu",
         C.c_char_p_p, "seqs", C.c_int, "nseqs", C.c_int, "seqlen"])

    # Argweaver Forward algorithm.
    argweaver_forward_alg = export(
        argweaverclib, "arghmm_forward_alg", C.c_double_p_p,
        [C.c_void_p, "trees", C.c_double_list, "times", C.c_int, "ntimes",
         C.c_double_list, "popsizes", C.c_double, "rho", C.c_double, "mu",
         C.c_char_p_p, "seqs", C.c_int, "nseqs", C.c_int, "seqlen",
         C.c_bool, "prior_given", C.c_double_list, "prior",
         C.c_bool, "internal", C.c_bool, "slow"])
    delete_double_matrix = export(
        argweaverclib, "delete_double_matrix", C.c_int,
        [C.c_double_p_p, "mat", C.c_int, "nrows"])
    delete_forward_matrix = export(
        argweaverclib, "delete_forward_matrix", C.c_int,
        [C.c_double_p_p, "mat", C.c_int, "nrows"])

    # ARG thread sampling.
    argweaver_sample_posterior = export(
        argweaverclib, "arghmm_sample_posterior", POINTER(C.c_int * 2),
        [C.c_int_matrix, "ptrees", C.c_int_matrix, "ages",
         C.c_int_matrix, "sprs", C.c_int_list, "blocklens",
         C.c_int, "ntrees", C.c_int, "nnodes",
         C.c_double_list, "times", C.c_int, "ntimes",
         C.c_double_list, "popsizes", C.c_double, "rho", C.c_double, "mu",
         C.c_char_p_p, "seqs", C.c_int, "nseqs", C.c_int, "seqlen",
         POINTER(POINTER(C.c_int * 2)), "path"])
    argweaver_sample_thread = export(
        argweaverclib, "arghmm_sample_thread", C.c_void_p,
        [C.c_void_p, "trees", C.c_double_list, "times", C.c_int, "ntimes",
         C.c_double_list, "popsizes", C.c_double, "rho", C.c_double, "mu",
         C.c_char_p_p, "seqs", C.c_int, "nseqs", C.c_int, "seqlen"])
    argweaver_sample_arg_thread_internal = export(
        argweaverclib, "arghmm_sample_arg_thread_internal", C.c_int,
        [C.c_void_p, "trees", C.c_double_list, "times", C.c_int, "ntimes",
         C.c_double_list, "popsizes", C.c_double, "rho", C.c_double, "mu",
         C.c_char_p_p, "seqs", C.c_int, "nseqs", C.c_int, "seqlen",
         C.c_out(C.c_int_list), "thread_path"])

    # ARG sampling.
    argweaver_sample_arg_seq = export(
        argweaverclib, "arghmm_sample_arg_seq", C.c_void_p,
        [C.c_double_list, "times", C.c_int, "ntimes",
         C.c_double_list, "popsizes", C.c_double, "rho", C.c_double, "mu",
         C.c_char_p_p, "seqs", C.c_int, "nseqs", C.c_int, "seqlen"])
    argweaver_sample_arg_refine = export(
        argweaverclib, "arghmm_sample_arg_refine", C.c_void_p,
        [C.c_double_list, "times", C.c_int, "ntimes",
         C.c_double_list, "popsizes", C.c_double, "rho", C.c_double, "mu",
         C.c_char_p_p, "seqs", C.c_int, "nseqs", C.c_int, "seqlen",
         C.c_int, "niters", C.c_int, "nremove"])
    argweaver_resample_arg = export(
        argweaverclib, "arghmm_resample_arg", C.c_void_p,
        [C.c_void_p, "trees", C.c_double_list, "times", C.c_int, "ntimes",
         C.c_double_list, "popsizes", C.c_double, "rho", C.c_double, "mu",
         C.c_char_p_p, "seqs", C.c_int, "nseqs", C.c_int, "seqlen",
         C.c_int, "niters", C.c_int, "nremove"])
    argweaver_resample_all_arg = export(
        argweaverclib, "arghmm_resample_all_arg", C.c_void_p,
        [C.c_void_p, "trees", C.c_double_list, "times", C.c_int, "ntimes",
         C.c_double_list, "popsizes", C.c_double, "rho", C.c_double, "mu",
         C.c_char_p_p, "seqs", C.c_int, "nseqs", C.c_int, "seqlen",
         C.c_int, "niters", C.c_double, "prob_path_switch"])
    argweaver_resample_mcmc_arg = export(
        argweaverclib, "arghmm_resample_mcmc_arg", C.c_void_p,
        [C.c_void_p, "trees", C.c_double_list, "times", C.c_int, "ntimes",
         C.c_double_list, "popsizes", C.c_double, "rho", C.c_double, "mu",
         C.c_char_p_p, "seqs", C.c_int, "nseqs", C.c_int, "seqlen",
         C.c_int, "niters", C.c_int, "niters2", C.c_int, "window"])
    argweaver_resample_climb_arg = export(
        argweaverclib, "arghmm_resample_climb_arg", C.c_void_p,
        [C.c_void_p, "trees", C.c_double_list, "times", C.c_int, "ntimes",
         C.c_double_list, "popsizes", C.c_double, "rho", C.c_double, "mu",
         C.c_char_p_p, "seqs", C.c_int, "nseqs", C.c_int, "seqlen",
         C.c_int, "niters", C.c_double, "recomb_preference"])
    #argweaver_resample_arg_cut = export(
    #    argweaverclib, "arghmm_resample_arg_cut", C.c_void_p,
    #    [C.c_void_p, "trees", C.c_double_list, "times", C.c_int, "ntimes",
    #     C.c_double_list, "popsizes", C.c_double, "rho", C.c_double, "mu",
    #     C.c_char_p_p, "seqs", C.c_int, "nseqs", C.c_int, "seqlen",
    #     C.c_int, "niters"])
    argweaver_resample_arg_region = export(
        argweaverclib, "arghmm_resample_arg_region", C.c_void_p,
        [C.c_void_p, "trees", C.c_double_list, "times", C.c_int, "ntimes",
         C.c_double_list, "popsizes", C.c_double, "rho", C.c_double, "mu",
         C.c_char_p_p, "seqs", C.c_int, "nseqs", C.c_int, "seqlen",
         C.c_int, "region_start", C.c_int, "region_end", C.c_int, "niters"])

    # ARG probability.
    argweaver_likelihood = export(
        argweaverclib, "arghmm_likelihood", C.c_double,
        [C.c_void_p, "trees", C.c_double_list, "times", C.c_int, "ntimes",
         C.c_double, "mu",
         C.c_char_p_p, "seqs", C.c_int, "nseqs", C.c_int, "seqlen"])
    argweaver_likelihood_parsimony = export(
        argweaverclib, "arghmm_likelihood_parsimony", C.c_double,
        [C.c_void_p, "trees", C.c_double_list, "times", C.c_int, "ntimes",
         C.c_double, "mu",
         C.c_char_p_p, "seqs", C.c_int, "nseqs", C.c_int, "seqlen"])
    argweaver_joint_prob = export(
        argweaverclib, "arghmm_joint_prob", C.c_double,
        [C.c_void_p, "trees", C.c_double_list, "times", C.c_int, "ntimes",
         C.c_double_list, "popsizes",
         C.c_double, "mu", C.c_double, "rho",
         C.c_char_p_p, "seqs", C.c_int, "nseqs", C.c_int, "seqlen"])
    argweaver_prior_prob = export(
        argweaverclib, "arghmm_prior_prob", C.c_double,
        [C.c_void_p, "trees",
         C.c_double_list, "times", C.c_int, "ntimes",
         C.c_double_list, "popsizes", C.c_double, "rho"])
    argweaver_tree_prior_prob = export(
        argweaverclib, "arghmm_tree_prior_prob", C.c_double,
        [C.c_void_p, "trees", C.c_double_list, "times", C.c_int, "ntimes",
         C.c_double_list, "popsizes"])
    #prob_coal_counts_matrix = export(
    #    argweaverclib, "prob_coal_counts_matrix", C.c_double,
    #    [C.c_int, "a", C.c_int, "b", C.c_double, "t", C.c_double, "n"])

    # Estimating population sizes.
    argweaver_est_popsizes_trees = export(
        argweaverclib, "arghmm_est_popsizes_trees", C.c_double,
        [C.c_void_p, "trees",
         C.c_double_list, "times", C.c_int, "ntimes", C.c_int, "step",
         C.c_out(C.c_double_list), "popsizes"])

    # Threading.
    argweaver_sample_arg_removal_path = export(
        argweaverclib, "arghmm_sample_arg_removal_path", C.c_int,
        [C.c_void_p, "trees", C.c_int, "node",
         C.c_out(C.c_int_list), "path"])
    argweaver_sample_arg_removal_leaf_path = export(
        argweaverclib, "arghmm_sample_arg_removal_leaf_path", C.c_int,
        [C.c_void_p, "trees", C.c_int, "node",
         C.c_out(C.c_int_list), "path"])
    argweaver_sample_arg_removal_path2 = export(
        argweaverclib, "arghmm_sample_arg_removal_path2", C.c_int,
        [C.c_void_p, "trees", C.c_int, "node", C.c_int, "pos",
         C.c_out(C.c_int_list), "path"])
    argweaver_sample_arg_removal_path_recomb = export(
        argweaverclib, "arghmm_sample_arg_removal_path_recomb", C.c_int,
        [C.c_void_p, "trees", C.c_double, "recomb_preference",
         C.c_out(C.c_int_list), "path"])
    argweaver_remove_arg_thread_path = export(
        argweaverclib, "arghmm_remove_arg_thread_path", C.c_int,
        [C.c_void_p, "trees", C.c_int_list, "path", C.c_int, "maxtime"])
    argweaver_remove_arg_thread_path2 = export(
        argweaverclib, "arghmm_remove_arg_thread_path2", C.c_int,
        [C.c_void_p, "trees", C.c_int_list, "path", C.c_int, "maxtime",
         C.c_out(C.c_int_list), "original_thread"])
    argweaver_get_thread_times = export(
        argweaverclib, "arghmm_get_thread_times", C.c_int,
        [C.c_void_p, "trees", C.c_int, "ntimes", C.c_int_list, "path",
         C.c_out(C.c_int_list), "path_times"])

    # ARG data structure API.
    argweaver_new_trees = export(
        argweaverclib, "arghmm_new_trees", C.c_void_p,
        [C.c_int_matrix, "ptrees", C.c_int_matrix, "ages",
         C.c_int_matrix, "sprs", C.c_int_list, "blocklens",
         C.c_int, "ntrees", C.c_int, "nnodes", C.c_int, "start_coord"])
    argweaver_copy_trees = export(
        argweaverclib, "arghmm_copy_trees", C.c_void_p,
        [C.c_void_p, "trees"])
    get_local_trees_ntrees = export(
        argweaverclib, "get_local_trees_ntrees", C.c_int,
        [C.c_void_p, "trees"])
    get_local_trees_nnodes = export(
        argweaverclib, "get_local_trees_nnodes", C.c_int,
        [C.c_void_p, "trees"])
    get_local_trees_ptrees = export(
        argweaverclib, "get_local_trees_ptrees", C.c_int,
        [C.c_void_p, "trees", C.c_out(C.c_int_matrix), "ptrees",
         C.c_out(C.c_int_matrix), "ages",
         C.c_out(C.c_int_matrix), "sprs",
         C.c_out(C.c_int_list), "blocklens"])
    get_treelens = export(
        argweaverclib, "get_treelens", C.c_int,
        [C.c_void_p, "trees", C.c_double_list, "times", C.c_int, "ntimes",
         C.c_out(C.c_double_list), "treelens"])
    get_local_trees_blocks = export(
        argweaverclib, "get_local_trees_blocks", C.c_int,
        [C.c_void_p, "trees", C.c_out(C.c_int_list), "starts",
         C.c_out(C.c_int_list), "ends"])
    delete_local_trees = export(
        argweaverclib, "delete_local_trees", C.c_int,
        [C.c_void_p, "trees"])
    write_local_trees = export(
        argweaverclib, "write_local_trees", C.c_int,
        [C.c_char_p, "filename", C.c_void_p, "trees",
         C.c_char_p_list, "names",
         C.c_double_list, "times", C.c_int, "ntimes"])
    read_local_trees = export(
        argweaverclib, "read_local_trees", C.c_void_p,
        [C.c_char_p, "filename",
         C.c_double_list, "times", C.c_int, "ntimes"])

    # Thread data structures.
    delete_path = export(
        argweaverclib, "delete_path", C.c_int,
        [POINTER(C.c_int * 2), "path"])
    argweaver_get_nstates = export(
        argweaverclib, "arghmm_get_nstates", int,
        [C.c_void_p, "trees",  C.c_int, "ntimes", C.c_bool, "internal",
         C.c_out(C.c_int_list), "nstates"])
    get_state_spaces = export(
        argweaverclib, "get_state_spaces", POINTER(POINTER(C.c_int * 2)),
        [C.c_void_p, "trees",  C.c_int, "ntimes", C.c_bool, "internal"])
    delete_state_spaces = export(
        argweaverclib, "delete_state_spaces", C.c_int,
        [POINTER(POINTER(C.c_int * 2)), "all_states", C.c_int, "ntrees"])

    setLogLevel = export(
        argweaverclib, "setLogLevel", C.c_int,
        [C.c_int, "level"])


# By default use a random seed.
if argweaverclib:
    argweaverclib.srand(int((time.time() * 1000) % 1e9))
    argweaverclib.setLogLevel(1)


def set_random_seed(num):
    """Set the C random number generator seed"""
    argweaverclib.srand(num)


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
    treelen = argweaver.get_treelen(tree, times)
    transmat = new_transition_probs(
        len(nodes), ptree, ages_index, treelen,
        ((C.c_int * 2) * nstates)
        (* ((C.c_int * 2)(n, t) for n, t in int_states)), nstates,
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
        argweaver.find_recomb_coal(tree, last_tree, recomb_name=recomb_name)

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
    assert nodes[i] is None
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
        ((C.c_int * 2) * nstates1)
        (* ((C.c_int * 2)(n, t) for n, t in int_states1)), nstates1,
        ((C.c_int * 2) * nstates2)
        (* ((C.c_int * 2)(n, t) for n, t in int_states2)), nstates2,

        len(time_steps), times, time_steps,
        nbranches, nrecombs, ncoals,
        popsizes, rho)

    if raw:
        return transmat
    else:
        transmat2 = [transmat[j][:nstates2]
                     for j in range(nstates1)]
        delete_transition_probs(transmat, nstates1)
        return transmat2


def delete_trans_emit_matrices(matrices):
    """Delete matrices"""
    for block, nstates, transmat, transmat_switch, emit in matrices:
        delete_emissions(emit, block[1] - block[0])
        delete_transition_probs(transmat, nstates)


def assert_transition_probs(tree, times, popsizes, rho):

    times_lookup = dict((t, i) for i, t in enumerate(times))
    tree2 = tree.get_tree()
    ptree, nodes, nodelookup = make_ptree(tree2)
    ages = [times_lookup[tree[node.name].age] for node in nodes]

    return argweaver_assert_transmat(len(ptree), ptree, ages,
                                     len(times), times, popsizes, rho)


def assert_transition_probs_internal(tree, times, popsizes, rho):

    times_lookup = dict((t, i) for i, t in enumerate(times))
    tree2 = tree.get_tree()
    ptree, nodes, nodelookup = make_ptree(tree2)
    ages = [times_lookup[tree[node.name].age] for node in nodes]

    return argweaver_assert_transmat_internal(len(ptree), ptree, ages,
                                              len(times), times, popsizes, rho)


def assert_transition_switch_probs(tree, spr, times, popsizes, rho):

    times_lookup = dict((t, i) for i, t in enumerate(times))
    tree2 = tree.get_tree()
    ptree, nodes, nodelookup = make_ptree(tree2)
    ages = [times_lookup[tree[node.name].age] for node in nodes]

    (r, rt), (c, ct) = spr
    recomb_name = nodelookup[tree2[r]]
    coal_name = nodelookup[tree2[c]]
    recomb_time = times_lookup[rt]
    coal_time = times_lookup[ct]

    return argweaver_assert_transmat_switch(
        len(ptree), ptree, ages,
        recomb_name, recomb_time, coal_name, coal_time,
        len(times), times, popsizes, rho)


def assert_transition_probs_switch_internal(trees, times, popsizes, rho):

    return argweaver_assert_transmat_switch_internal(
        trees, len(times), times, popsizes, rho)


#=============================================================================


def argweaver_forward_algorithm(arg, seqs, rho=1.5e-8,
                                mu=2.5e-8, popsizes=1e4, times=None,
                                ntimes=20, maxtime=180000,
                                verbose=False,
                                prior=[], internal=False, slow=False):
    if times is None:
        times = argweaver.get_time_points(
            ntimes=ntimes, maxtime=maxtime, delta=.01)
    if isinstance(popsizes, float) or isinstance(popsizes, int):
        popsizes = [popsizes] * len(times)

    probs = []

    if verbose:
        util.tic("forward")

    if is_carg(arg):
        trees, names = arg
    else:
        trees, names = arg2ctrees(arg, times)

    seqs2 = [seqs[node] for node in names]
    for name in seqs.keys():
        if name not in names:
            seqs2.append(seqs[name])
    seqlen = len(seqs2[0])

    fw = argweaver_forward_alg(trees, times, len(times),
                               popsizes, rho, mu,
                               (C.c_char_p * len(seqs2))(*seqs2), len(seqs2),
                               seqlen, len(prior) > 0, prior, internal,
                               slow)

    nstates = [0] * seqlen
    argweaver_get_nstates(trees, len(times), internal, nstates)

    probs = [row[:n] for row, n in zip(fw, nstates)]

    delete_forward_matrix(fw, seqlen)

    if verbose:
        util.toc()

    return probs


#=============================================================================
# sampling ARG threads


def sample_thread(arg, seqs, rho=1.5e-8, mu=2.5e-8, popsize=1e4,
                  times=None, ntimes=20, maxtime=200000, verbose=False):

    if times is None:
        times = argweaver.get_time_points(
            ntimes=ntimes, maxtime=maxtime, delta=.01)
    popsizes = [popsize] * len(times)

    if verbose:
        util.tic("sample thread")

    trees, names = arg2ctrees(arg, times)

    seqs2 = [seqs[name] for name in names]

    new_name = [x for x in seqs.keys() if x not in names][0]
    names.append(new_name)
    seqs2.append(seqs[new_name])
    seqlen = len(seqs2[0])

    trees = argweaver_sample_thread(
        trees, times, len(times),
        popsizes, rho, mu,
        (C.c_char_p * len(seqs2))(*seqs2), len(seqs2), seqlen, None)
    arg = ctrees2arg(trees, names, times, verbose=verbose)

    if verbose:
        util.toc()

    return arg


'''
def sample_posterior(model, n, verbose=False):

    if verbose:
        util.tic("sample thread")

    (ptrees, ages, sprs, blocks), all_nodes = get_treeset(
        model.arg, model.times)
    blocklens = [x[1] - x[0] for x in blocks]
    seqlen = sum(blocklens)

    seqs = [model.seqs[node] for node in all_nodes[0]
            if model.arg[node].is_leaf()]
    seqs.append(model.seqs[model.new_name])
    path = argweaver_sample_posterior(
        ptrees, ages, sprs, blocklens,
        len(ptrees), len(ptrees[0]),
        model.times, len(model.times),
        model.popsizes, model.rho, model.mu,
        (C.c_char_p * len(seqs))(*seqs), len(seqs),
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
'''


#=============================================================================
# ARG sampling


def sample_arg(seqs, ntimes=20, rho=1.5e-8, mu=2.5e-8, popsizes=1e4,
               refine=0, nremove=1, times=None, verbose=False,
               carg=False):
    """
    Sample ARG for sequences
    """
    if times is None:
        times = argweaver.get_time_points(
            ntimes=ntimes, maxtime=80000, delta=.01)
    if isinstance(popsizes, float) or isinstance(popsizes, int):
        popsizes = [popsizes] * len(times)

    if verbose:
        util.tic("sample arg")

    names = []
    seqs2 = []
    for name, seq in seqs.items():
        names.append(name)
        seqs2.append(seq)

    # sample arg
    trees = argweaver_sample_arg_refine(
        times, len(times),
        popsizes, rho, mu,
        (C.c_char_p * len(seqs))(*seqs2), len(seqs), len(seqs2[0]), refine,
        nremove)

    if carg:
        arg = (trees, names)
    else:
        # convert to python
        arg = ctrees2arg(trees, names, times, verbose=verbose)

    if verbose:
        util.toc()

    return arg


def resample_arg(arg, seqs, ntimes=20, rho=1.5e-8, mu=2.5e-8, popsizes=1e4,
                 refine=1, nremove=1, times=None, verbose=False, carg=False):
    """
    Sample ARG for sequences
    """
    if times is None:
        times = argweaver.get_time_points(
            ntimes=ntimes, maxtime=80000, delta=.01)
    if isinstance(popsizes, float) or isinstance(popsizes, int):
        popsizes = [popsizes] * len(times)

    if verbose:
        util.tic("resample arg")

    # convert arg to c++
    if verbose:
        util.tic("convert arg")
    trees, names = arg2ctrees(arg, times)
    if verbose:
        util.toc()

    # get sequences in same order
    # and add all other sequences not in arg yet
    seqs2 = [seqs[name] for name in names]
    leaves = set(names)
    names = list(names)
    for name, seq in seqs.items():
        if name not in leaves:
            names.append(name)
            seqs2.append(seq)

    # resample arg
    seqlen = len(seqs[names[0]])
    trees = argweaver_resample_arg(
        trees, times, len(times),
        popsizes, rho, mu,
        (C.c_char_p * len(seqs2))(*seqs2), len(seqs2),
        seqlen, refine, nremove)

    if carg:
        arg = (trees, names)
    else:
        # convert arg back to python
        arg = ctrees2arg(trees, names, times, verbose=verbose)

    if verbose:
        util.toc()

    return arg


def sample_all_arg(seqs, ntimes=20, rho=1.5e-8, mu=2.5e-8, popsizes=1e4,
                   refine=1, times=None, verbose=False, carg=False,
                   prob_path_switch=.1):
    """
    Sample ARG for sequences
    """
    if times is None:
        times = argweaver.get_time_points(
            ntimes=ntimes, maxtime=80000, delta=.01)
    if isinstance(popsizes, float) or isinstance(popsizes, int):
        popsizes = [popsizes] * len(times)

    if verbose:
        util.tic("resample arg")

    # convert arg to c++
    if verbose:
        util.tic("convert arg")

    arg = argweaver.make_trunk_arg(
        0, len(seqs.values()[0]), name=seqs.keys()[0])
    trees, names = arg2ctrees(arg, times)
    if verbose:
        util.toc()

    # get sequences in same order
    # and add all other sequences not in arg yet
    seqs2 = [seqs[name] for name in names]
    leaves = set(names)
    for name, seq in seqs.items():
        if name not in leaves:
            names.append(name)
            seqs2.append(seq)

    # resample arg
    seqlen = len(seqs[names[0]])
    trees = argweaver_resample_all_arg(
        trees, times, len(times),
        popsizes, rho, mu,
        (C.c_char_p * len(seqs2))(*seqs2), len(seqs2),
        seqlen, refine, prob_path_switch)

    if carg:
        arg = (trees, names)
    else:
        # convert arg back to python
        arg = ctrees2arg(trees, names, times, verbose=verbose)

    if verbose:
        util.toc()

    return arg


def resample_all_arg(arg, seqs, ntimes=20, rho=1.5e-8, mu=2.5e-8, popsizes=1e4,
                     refine=1, times=None, verbose=False, carg=False,
                     prob_path_switch=.1):
    """
    Sample ARG for sequences
    """
    if times is None:
        times = argweaver.get_time_points(
            ntimes=ntimes, maxtime=80000, delta=.01)
    if isinstance(popsizes, float) or isinstance(popsizes, int):
        popsizes = [popsizes] * len(times)

    if verbose:
        util.tic("resample arg")

    # convert arg to c++
    if verbose:
        util.tic("convert arg")
    trees, names = arg2ctrees(arg, times)
    if verbose:
        util.toc()

    # get sequences in same order
    # and add all other sequences not in arg yet
    leaves = set(names)
    names = list(names)
    for name in seqs:
        if name not in leaves:
            names.append(name)
    seqs2, nseqs, seqlen = seqs2cseqs(seqs, names)

    # resample arg
    trees = argweaver_resample_all_arg(
        trees, times, len(times),
        popsizes, rho, mu,
        seqs2, nseqs, seqlen, refine, prob_path_switch)

    if carg:
        arg = (trees, names)
    else:
        # convert arg back to python
        arg = ctrees2arg(trees, names, times, verbose=verbose)

    if verbose:
        util.toc()

    return arg


def resample_climb_arg(arg, seqs, ntimes=20, rho=1.5e-8, mu=2.5e-8,
                       popsizes=1e4, refine=1, recomb_pref=.7,
                       times=None, verbose=False, carg=False):
    """
    Sample ARG for sequences
    """
    if times is None:
        times = argweaver.get_time_points(
            ntimes=ntimes, maxtime=80000, delta=.01)
    if isinstance(popsizes, float) or isinstance(popsizes, int):
        popsizes = [popsizes] * len(times)

    if verbose:
        util.tic("resample arg")

    # convert arg to c++
    if verbose:
        util.tic("convert arg")
    trees, names = arg2ctrees(arg, times)
    if verbose:
        util.toc()

    # get sequences in same order
    # and add all other sequences not in arg yet
    leaves = set(names)
    for name, seq in seqs.items():
        if name not in leaves:
            names.append(name)
    seqs2, nseqs, seqlen = seqs2cseqs(seqs, names)

    # resample arg
    trees = argweaver_resample_climb_arg(
        trees, times, len(times),
        popsizes, rho, mu, seqs2, nseqs, seqlen, refine, recomb_pref)

    if carg:
        arg = (trees, names)
    else:
        # convert arg back to python
        arg = ctrees2arg(trees, names, times, verbose=verbose)

    if verbose:
        util.toc()

    return arg


def resample_mcmc_arg(arg, seqs, ntimes=20,
                      rho=1.5e-8, mu=2.5e-8, popsizes=1e4,
                      refine=1, times=None, verbose=False, carg=False,
                      window=200000, niters2=5):
    """
    Sample ARG for sequences
    """
    if times is None:
        times = argweaver.get_time_points(
            ntimes=ntimes, maxtime=80000, delta=.01)
    if isinstance(popsizes, float) or isinstance(popsizes, int):
        popsizes = [popsizes] * len(times)

    if verbose:
        util.tic("resample arg")

    # convert arg to c++
    if verbose:
        util.tic("convert arg")
    trees, names = arg2ctrees(arg, times)
    if verbose:
        util.toc()

    # get sequences in same order
    # and add all other sequences not in arg yet
    leaves = set(names)
    names = list(names)
    for name in seqs:
        if name not in leaves:
            names.append(name)
    seqs2, nseqs, seqlen = seqs2cseqs(seqs, names)

    # resample arg
    trees = argweaver_resample_mcmc_arg(
        trees, times, len(times),
        popsizes, rho, mu,
        seqs2, nseqs, seqlen, refine, niters2, window)

    if carg:
        arg = (trees, names)
    else:
        # convert arg back to python
        arg = ctrees2arg(trees, names, times, verbose=verbose)

    if verbose:
        util.toc()

    return arg


'''
def resample_arg_cut(
        arg, seqs, ntimes=20, rho=1.5e-8, mu=2.5e-8, popsizes=1e4,
        refine=1, times=None, verbose=False, carg=False):
    """
    Sample ARG for sequences
    """
    if times is None:
        times = argweaver.get_time_points(
            ntimes=ntimes, maxtime=80000, delta=.01)
    if isinstance(popsizes, float) or isinstance(popsizes, int):
        popsizes = [popsizes] * len(times)

    if verbose:
        util.tic("resample arg")

    # convert arg to c++
    if verbose:
        util.tic("convert arg")
    trees, names = arg2ctrees(arg, times)
    if verbose:
        util.toc()

    # get sequences in same order
    # and add all other sequences not in arg yet
    leaves = set(names)
    names = list(names)
    for name in seqs:
        if name not in leaves:
            names.append(name)
    seqs2, nseqs, seqlen = seqs2cseqs(seqs, names)

    # resample arg
    trees = argweaver_resample_arg_cut(trees, times, len(times),
                                       popsizes, rho, mu,
                                       seqs2, nseqs, seqlen, refine)

    if carg:
        arg = (trees, names)
    else:
        # convert arg back to python
        arg = ctrees2arg(trees, names, times, verbose=verbose)

    if verbose:
        util.toc()

    return arg
'''


def resample_arg_region(arg, seqs, region_start, region_end,
                        ntimes=20, rho=1.5e-8, mu=2.5e-8,
                        popsizes=1e4, times=None, carg=False,
                        refine=1, verbose=False):
    """
    Sample ARG for sequences
    """
    if times is None:
        times = argweaver.get_time_points(
            ntimes=ntimes, maxtime=80000, delta=.01)
    if isinstance(popsizes, float) or isinstance(popsizes, int):
        popsizes = [popsizes] * len(times)

    if verbose:
        util.tic("resample arg")

    # convert arg to c++
    if verbose:
        util.tic("convert arg")
    trees, names = arg2ctrees(arg, times)
    if verbose:
        util.toc()

    # get sequences in same order
    # and add all other sequences not in arg yet
    leaves = set(names)
    for name, seq in seqs.items():
        if name not in leaves:
            names.append(name)
    seqs2, nseqs, seqlen = seqs2cseqs(seqs, names)

    # resample arg
    seqlen = len(seqs[names[0]])

    trees = argweaver_resample_arg_region(
        trees, times, len(times),
        popsizes, rho, mu, seqs2, nseqs, seqlen,
        region_start, region_end, refine)

    #trees = argweaver_resample_arg_region(
    #    trees, times, len(times),
    #    popsizes, rho, mu, seqs2, nseqs, seqlen,
    #    region_start, region_end)

    # convert arg back to python
    if carg:
        arg = (trees, names)
    else:
        arg = ctrees2arg(trees, names, times, verbose=verbose)

    if verbose:
        util.toc()

    return arg


def resample_arg_regions(arg, seqs, niters, width=1000,
                         ntimes=20, rho=1.5e-8, mu=2.5e-8,
                         popsize=1e4, times=None, carg=False,
                         verbose=False):
    seqlen = len(seqs.values()[0])

    if is_carg(arg):
        trees, names = arg
        arg2 = ctrees2arg(trees, names, times, verbose=verbose,
                          delete_arg=False)
        recomb_pos = list(x.pos for x in arg2 if x.event == "recomb")
    else:
        recomb_pos = list(x.pos for x in arg if x.event == "recomb")

    for it in range(niters):
        maxr = 0
        for i, j, a, b in stats.iter_window_index(recomb_pos, width):
            r = j - i + 1
            if r > maxr:
                maxr = r
                region = [max(recomb_pos[i]-10, 10),
                          min(recomb_pos[j]+10, seqlen - 10)]

        if verbose:
            util.tic("sample ARG region %s" % region)
        print arg
        arg = argweaver.resample_arg_region(arg, seqs, region[0], region[1],
                                            rho=rho, mu=mu, times=times,
                                            carg=carg, verbose=True)
        if not carg:
            recomb_pos = list(x.pos for x in arg if x.event == "recomb")
            if verbose:
                util.logger("%d: # recombs %d" % (it, len(recomb_pos)))
        if verbose:
            util.toc()

    return arg


#=============================================================================
# ARG probabilities


def calc_likelihood(arg, seqs, ntimes=20, mu=2.5e-8,
                    times=None, delete_arg=True, verbose=False):
    """
    Calculate arg_likelihood
    """
    if times is None:
        times = argweaver.get_time_points(
            ntimes=ntimes, maxtime=80000, delta=.01)

    if verbose:
        util.tic("calc likelihood")

    trees, names = arg2ctrees(arg, times)
    seqs, nseqs, seqlen = seqs2cseqs(seqs, names)

    lk = argweaver_likelihood(
        trees, times, len(times), mu, seqs, nseqs, seqlen)
    if delete_arg:
        delete_local_trees(trees)

    if verbose:
        util.toc()

    return lk


def calc_likelihood_parsimony(arg, seqs, ntimes=20, mu=2.5e-8,
                              times=None, delete_arg=True, verbose=False):
    """
    Calculate arg_likelihood
    """
    if times is None:
        times = argweaver.get_time_points(
            ntimes=ntimes, maxtime=80000, delta=.01)

    if verbose:
        util.tic("calc likelihood")

    trees, names = arg2ctrees(arg, times)
    seqs, nseqs, seqlen = seqs2cseqs(seqs, names)

    lk = argweaver_likelihood_parsimony(
        trees, times, len(times), mu, seqs, nseqs, seqlen)
    if delete_arg:
        delete_local_trees(trees)

    if verbose:
        util.toc()

    return lk


def calc_prior_prob(arg, ntimes=20, rho=1.5e-8, popsizes=1e4,
                    times=None, delete_arg=True, verbose=False):
    """
    Calculate arg_joint_prob
    """
    if times is None:
        times = argweaver.get_time_points(
            ntimes=ntimes, maxtime=80000, delta=.01)
    if isinstance(popsizes, float) or isinstance(popsizes, int):
        popsizes = [popsizes] * len(times)

    if verbose:
        util.tic("calc likelihood")

    trees, names = arg2ctrees(arg, times)

    p = argweaver_prior_prob(trees, times, len(times), popsizes, rho)
    if delete_arg:
        delete_local_trees(trees)

    if verbose:
        util.toc()

    return p


def calc_joint_prob(arg, seqs, ntimes=20, mu=2.5e-8, rho=1.5e-8, popsizes=1e4,
                    times=None, verbose=False, delete_arg=True):
    """
    Calculate arg_joint_prob
    """
    if times is None:
        times = argweaver.get_time_points(
            ntimes=ntimes, maxtime=80000, delta=.01)
    if isinstance(popsizes, float) or isinstance(popsizes, int):
        popsizes = [popsizes] * len(times)

    if verbose:
        util.tic("calc likelihood")

    trees, names = arg2ctrees(arg, times)
    seqs, nseqs, seqlen = seqs2cseqs(seqs, names)

    p = argweaver_joint_prob(
        trees, times, len(times), popsizes, mu, rho, seqs, nseqs, seqlen)
    if delete_arg:
        delete_local_trees(trees)

    if verbose:
        util.toc()

    return p


#=============================================================================
def est_popsizes_trees(arg, times, step, verbose=False):

    if verbose:
        util.tic("convert arg")
    trees, names = arg2ctrees(arg, times)

    if verbose:
        util.toc()
        util.tic("estimate popsizes")

    popsizes = [0.0] * (len(times) - 1)
    argweaver_est_popsizes_trees(trees, times, len(times), step, popsizes)

    if verbose:
        util.toc()

    if not is_carg(arg):
        delete_local_trees(trees)

    return popsizes


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


def is_carg(arg):
    return isinstance(arg, tuple)


def arg2ctrees(arg, times):

    # check to see if arg is already converted
    if is_carg(arg):
        return arg

    (ptrees, ages, sprs, blocks), all_nodes = get_treeset(
        arg, times)
    blocklens = [x[1] - x[0] for x in blocks]
    names = [node for node in all_nodes[0]
             if arg[node].is_leaf()]

    trees = argweaver_new_trees(
        ptrees, ages, sprs, blocklens,
        len(ptrees), len(ptrees[0]), arg.start)

    return trees, names


def ctrees2arg(trees, names, times, verbose=False, delete_arg=True):
    """
    Convert a C data structure for the ARG into a python ARG
    """

    if verbose:
        util.tic("convert arg")

    # get local trees info
    nnodes = get_local_trees_nnodes(trees)
    ntrees = get_local_trees_ntrees(trees)

    # allocate data structures for treeset
    ptrees = []
    ages = []
    sprs = []
    blocklens = [0] * ntrees
    for i in range(ntrees):
        ptrees.append([0] * nnodes)
        ages.append([0] * nnodes)
        sprs.append([0, 0, 0, 0])

    # populate data structures
    get_local_trees_ptrees(trees, ptrees, ages, sprs, blocklens)

    # fully convert to python
    for i in range(ntrees):
        ptrees[i] = ptrees[i][:nnodes]
        ages[i] = ages[i][:nnodes]
        sprs[i] = sprs[i][:4]

    # convert treeset to arg data structure
    blocks = []
    start = 0
    for blocklen in blocklens:
        end = start + blocklen
        blocks.append((start, end))
        start = end

    assert len(names) == ((nnodes + 1) / 2)

    arg = treeset2arg(ptrees, ages, sprs, blocks, names, times)

    if delete_arg:
        delete_local_trees(trees)

    if verbose:
        util.toc()

    return arg


def iter_local_ptrees(arg, times, start=None, end=None):
    """
    Iterate through the local trees as ptrees (parent arrays)

    Node maintain the same index across trees until it is broken
    (parent of recomb node).  The recoal node takes the index of the broken
    node.
    """
    times_lookup = dict((t, i) for i, t in enumerate(times))

    last_tree2 = None
    last_nodes = None
    last_nodelookup = {}

    for block, tree, last_tree, spr in argweaver.iter_arg_sprs(
            arg, start, end):
        # get treelib.Tree from arglib.ARG
        tree2 = tree.get_tree()

        if last_tree is None:
            # get frist ptree
            ptree, nodes, nodelookup = make_ptree(tree2)
            ispr = [-1, -1, -1, -1]
            age = [times_lookup[tree[x.name].age] for x in nodes]

        else:
            (rname, rtime), (cname, ctime) = spr

            # assert ARG is SMC-style (no bubbles)
            assert rname != cname

            # find recoal node in tree2, it does not exist in last_tree2
            recomb_parent = last_tree2[rname].parent
            recoal = [x for x in tree2 if x.name not in last_tree2][0]

            # make nodes array consistent
            nodes = [tree2.nodes.get(x.name, None) for x in last_nodes]
            i = last_nodes.index(recomb_parent)
            # recomb_parent is 'broken' and should not be in nodes list
            assert nodes[i] is None

            # recoal takes broken node's spot in nodes list
            nodes[i] = recoal

            # get ptree
            ptree, nodes, nodelookup = make_ptree(tree2, nodes=nodes)
            age = [times_lookup[tree[x.name].age] for x in nodes]

            # get integer-based spr
            recomb_name = last_nodelookup[last_tree2[rname]]
            coal_name = last_nodelookup[last_tree2[cname]]
            ispr = [recomb_name, times_lookup[rtime],
                    coal_name, times_lookup[ctime]]

        yield ptree, age, ispr, block, nodes

        # setup last tree
        last_tree2 = tree2
        last_ptree, last_nodes, last_nodelookup = ptree, nodes, nodelookup


def get_treeset(arg, times, start=None, end=None):
    """
    Convert ARG into a ptree representation.

    Returns ((ptrees, ages, sprs, blocks), all_nodes).
    """
    ptrees = []
    ages = []
    sprs = []
    blocks = []
    all_nodes = []

    for ptree, age, ispr, block, nodes in iter_local_ptrees(
            arg, times, start=start, end=end):
        ptrees.append(ptree)
        ages.append(age)
        sprs.append(ispr)
        blocks.append(block)
        all_nodes.append([x.name for x in nodes])

    return (ptrees, ages, sprs, blocks), all_nodes


def treeset2arg(ptrees, ages, sprs, blocks, names, times):
    """
    Converts a ptree representation into an ARG.
    """

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
        pos = blocks[i][0]  # - 1

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

        # get leaf sets
        rleaves = lookup[rinode].leaf_names()
        cleaves = lookup[cinode].leaf_names()
        assert ritime >= ages[i-1][rinode], (pos, ritime, ages[i-1][rinode])
        assert citime >= ages[i-1][cinode], (pos, citime, ages[i-1][cinode])

        sprs2.append((pos, (rleaves, times[ritime]), (cleaves, times[citime])))

    arglib.make_arg_from_sprs(arg, sprs2)
    return arg


def seqs2cseqs(seqs, names):
    seqs2 = [seqs[name] for name in names]
    nseqs = len(seqs2)
    seqlen = len(seqs2[0])
    return (C.c_char_p * len(seqs2))(*seqs2), nseqs, seqlen

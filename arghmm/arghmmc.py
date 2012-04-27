#
# C interface for ArgHmm
#



# import arghmm C lib
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
            c_char_p_p, "seqs", c_int, "nseqs", c_int, "seqlen",
            c_double_p_p, "fw"])

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
            c_int, "niters"])


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


#=============================================================================
# helper functions for C interface




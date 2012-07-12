#ifndef ARGHMM_THREAD_H
#define ARGHMM_THREAD_H

#include "local_tree.h"
#include "states.h"

namespace arghmm {

// assert that a thread is compatiable with an ARG
bool assert_trees_thread(LocalTrees *trees, int *thread_path, int ntimes);

void add_tree_branch(LocalTree *tree, State state);

// update an SPR and mapping after adding a new branch
void add_spr_branch(LocalTree *tree, LocalTree *last_tree, 
                    State state, State last_state,
                    Spr *spr, int *mapping,
                    int newleaf, int displaced, int newcoal);

// add a thread to an ARG
void add_arg_thread(LocalTrees *trees, int ntimes, int *thread_path, int seqid,
                    vector<int> &recomb_pos, vector<NodePoint> &recombs);

// remove a thread from an ARG
void remove_arg_thread(LocalTrees *trees, int remove_seqid);


// Add a branch to a partial ARG
void add_arg_thread_path(LocalTrees *trees, int ntimes, int *thread_path, 
                         vector<int> &recomb_pos, vector<NodePoint> &recombs);
// Removes a thread path from an ARG and returns a partial ARG
void remove_arg_thread_path(LocalTrees *trees, const int *removal_path, 
                            int maxtime, int *original_thread=NULL);
void sample_arg_removal_path(LocalTrees *trees, int node, int *path);
void sample_arg_removal_leaf_path(LocalTrees *trees, int node, int *path);



} // namespace arghmm

#endif //ARGHMM_THREAD_H

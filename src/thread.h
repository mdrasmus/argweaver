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

} // namespace arghmm

#endif //ARGHMM_THREAD_H

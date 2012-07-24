#ifndef ARGHMM_THREAD_H
#define ARGHMM_THREAD_H

#include "local_tree.h"
#include "states.h"


namespace arghmm {

// assert that a thread is compatiable with an ARG
bool assert_trees_thread(LocalTrees *trees, int *thread_path, int ntimes);

// adding/removing branches from a local tree
void add_tree_branch(LocalTree *tree, int node, int time);
void remove_tree_branch(LocalTree *tree, int remove_leaf, int *displace);

// update an SPR and mapping after adding a new branch
void add_spr_branch(LocalTree *tree, LocalTree *last_tree, 
                    State state, State last_state,
                    Spr *spr, int *mapping,
                    int newleaf, int displaced, int newcoal);
void add_spr_branch(LocalTree *tree, LocalTree *last_tree, 
                    State state, State last_state,
                    Spr *spr, int *mapping,
                    int subtree_root, int last_subtree_root);


// add a thread to an ARG
void add_arg_thread(LocalTrees *trees, int ntimes, int *thread_path, int seqid,
                    vector<int> &recomb_pos, vector<NodePoint> &recombs);

// remove a thread from an ARG
void remove_arg_thread(LocalTrees *trees, int remove_seqid);



void get_next_removal_nodes(const LocalTree *tree1, const LocalTree *tree2,
                            const Spr &spr2, const int *mapping2,
                            int node, int next_nodes[2]);
void get_prev_removal_nodes(const LocalTree *tree1, const LocalTree *tree2,
                            const Spr &spr2, const int *mapping2,
                            int node, int prev_nodes[2]);


// Add a branch to a partial ARG
void add_arg_thread_path(LocalTrees *trees, int ntimes, int *thread_path, 
                         vector<int> &recomb_pos, vector<NodePoint> &recombs);
// Removes a thread path from an ARG and returns a partial ARG
void remove_arg_thread_path(LocalTrees *trees, const int *removal_path, 
                            int maxtime, int *original_thread=NULL);
void sample_arg_removal_path(LocalTrees *trees, int node, int *path);
void sample_arg_removal_path(LocalTrees *trees, int node, int pos, int *path,
                             double prob_switch=.1);
void sample_arg_removal_leaf_path(LocalTrees *trees, int node, int *path);

void sample_arg_removal_path_recomb(LocalTrees *trees, double recomb_preference,
                                    int *path);



} // namespace arghmm

#endif //ARGHMM_THREAD_H

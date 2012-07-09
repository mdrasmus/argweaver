
#include "thread.h"
#include "trans.h"

namespace arghmm {


// assert that a thread is compatiable with an ARG
bool assert_trees_thread(LocalTrees *trees, int *thread_path, int ntimes)
{
    LocalTree *last_tree = NULL;
    States states1, states2;
    States *states = &states1;
    States *last_states = &states2;

    // loop through blocks
    int start = trees->start_coord;
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it) {
        get_coal_states(it->tree, ntimes, *states);

        // check spr
        if (last_tree) {
            assert(assert_spr(last_tree, it->tree, &it->spr, it->mapping));

            int determ[last_states->size()];
            get_deterministic_transitions(
                 it->tree, last_tree, it->spr, it->mapping,
                 *last_states, *states, ntimes, determ);

            int a = thread_path[start-1];
            int b = determ[thread_path[start-1]];
            int c = thread_path[start];

            printf(">> %d (%d,%d) --> (%d,%d) == (%d,%d)\n", start, 
                   (*last_states)[a].node, (*last_states)[a].time, 
                   (*states)[b].node, (*states)[b].time,
                   (*states)[c].node, (*states)[c].time);
        }
        
        // set last tree and state pointers
        last_tree = it->tree;
        last_states = states;
        if (states == &states1)
            states = &states2;
        else
            states = &states1;

        start += it->blocklen; //.length();
    }

    return true;
}


// rename a node from src_node to dest_node while mantaining the 
// structure of the tree
void displace_node(LocalTree *tree, int src_node, int dest_node)
{
    // special case
    if (src_node == dest_node)
        return;

    LocalNode *nodes = tree->nodes;

    // copy displaced node
    nodes[dest_node].copy(nodes[src_node]);

    // notify parent of displacement
    int parent = nodes[dest_node].parent;
    if (parent != -1) {
        int *c = nodes[parent].child;
        if (c[0] == src_node)
            c[0] = dest_node;
        else
            c[1] = dest_node;
    }

    // notify children of displacement
    int *c = nodes[dest_node].child;
    if (c[0] != -1)
        nodes[c[0]].parent = dest_node;
    if (c[1] != -1)
        nodes[c[1]].parent = dest_node;
}


// Add a single new leaf to a tree.  
// Leaf connects to the point indicated by (node, time)
// 
//   - New leaf will be named 'nleaves' (newleaf)
//   - Node currently named 'nleaves' will be displaced to 'nnodes' (displaced)
//   - New internal node will be named 'nnodes+1' (newcoal)
//
void add_tree_branch(LocalTree *tree, int node, int time)
{
    // get tree info
    LocalNode *nodes = tree->nodes;
    int nleaves = tree->get_num_leaves();
    int nnodes = tree->nnodes;
    int nnodes2 = nnodes + 2;
    
    // get major node ids
    int newleaf = nleaves;
    int displaced = nnodes;
    int newcoal = nnodes + 1;
    
    // determine node displacement
    int node2 = (node != newleaf ? node : displaced);
    int parent = nodes[node].parent;
    int parent2 = (parent != newleaf ? parent : displaced);

    // displace node
    if (newleaf < displaced)
        displace_node(tree, newleaf, displaced);

    // add new leaf
    nodes[newleaf].parent = newcoal;
    nodes[newleaf].child[0] = -1;
    nodes[newleaf].child[1] = -1;
    nodes[newleaf].age = 0;
        
    // add new coal node
    nodes[newcoal].parent = parent2;
    nodes[newcoal].child[0] = newleaf;
    nodes[newcoal].child[1] = node2;
    nodes[newcoal].age = time;

    // fix pointers
    nodes[node2].parent = newcoal;
    if (parent2 != -1) {
        int *c = nodes[parent2].child;
        if (c[0] == node2)
            c[0] = newcoal;
        else
            c[1] = newcoal;
    }

    // fix up tree data
    tree->nnodes = nnodes2;
    tree->set_root();
    //assert(assert_tree(tree));
}


void remove_tree_branch(LocalTree *tree, int remove_leaf, int *displace)
{
    LocalNode *nodes = tree->nodes;
    int nnodes = tree->nnodes;
    int last_leaf = tree->get_num_leaves() - 1;

    // remove coal node
    int remove_coal = nodes[remove_leaf].parent;
    int *c = nodes[remove_coal].child;
    int coal_child = (c[0] == remove_leaf ? c[1] : c[0]);
    int coal_parent = nodes[remove_coal].parent;
    nodes[coal_child].parent = coal_parent;
    if (coal_parent != -1) {
        c = nodes[coal_parent].child;
        if (c[0] == remove_coal)
            c[0] = coal_child;
        else
            c[1] = coal_child;
    }
        
    // displace nodes
    for (int i=0; i<nnodes; i++)
        displace[i] = i;
    displace[remove_leaf] = -1;
    displace[remove_coal] = -1;
        
    // move last leaf into remove_leaf spot
    if (last_leaf != remove_leaf) {
        displace[last_leaf] = remove_leaf;
        displace_node(tree, last_leaf, remove_leaf);
    }

    // move nodes in nnodes-2 and nnodes-1 into holes
    int hole = last_leaf;
    if (remove_coal != nnodes-2) {
        displace[nnodes-2] = hole;
        displace_node(tree, nnodes-2, hole);
        hole = remove_coal;
    }
    if (remove_coal != nnodes-1) {
        displace[nnodes-1] = hole;
        displace_node(tree, nnodes-1, hole);
    }
    
    // set tree data
    tree->nnodes -= 2;
    tree->set_root();
    assert_tree(tree);
}



// update an SPR and mapping after adding a new branch
void add_spr_branch(LocalTree *tree, LocalTree *last_tree, 
                    State state, State last_state,
                    Spr *spr, int *mapping,
                    int newleaf, int displaced, int newcoal)
{
    // get tree info
    LocalNode *nodes = tree->nodes;
    LocalNode *last_nodes = last_tree->nodes;

    // determine node displacement
    int node2 = (state.node != newleaf ? state.node : displaced);


    // update mapping due to displacement
    mapping[displaced] = mapping[newleaf];
    mapping[newleaf] = newleaf;
            
    // set default new node mapping 
    mapping[newcoal] = newcoal;
        
    for (int i=newleaf+1; i<tree->nnodes; i++) {
        if (mapping[i] == newleaf)
            mapping[i] = displaced;
    }


    // update spr due to displacement
    if (spr->recomb_node == newleaf)
        spr->recomb_node = displaced;
    if (spr->coal_node == newleaf)
        spr->coal_node = displaced;

        
    // parent of recomb node should be the recoal point
    // however, if it equals newcoal, then either (1) the recomb branch is 
    // renamed, (2) there is mediation, or (3) new branch escapes
    int recoal = nodes[mapping[spr->recomb_node]].parent;
    if (recoal == newcoal) {
        if (mapping[last_state.node] == node2) {
            // (1) recomb is above coal state, we rename spr recomb node
            spr->recomb_node = newcoal;
        } else {
            // if this is a mediated coal, then state should equal recomb
            int state_node = (state.node != newleaf) ? 
                state.node : displaced;
            if (state_node == mapping[spr->recomb_node]) {
                // (3) this is a mediated coal, rename coal node and time
                spr->coal_node = newleaf;
                spr->coal_time = state.time;
            } else {
                // (2) this is the new branch escaping
                // no other updates are necessary
            }
        }
    } else {
        // the other possibility is that newcoal is under recoal point
        // if newcoal is child of recoal, then coal is renamed
        int *c = nodes[recoal].child;
        if (c[0] == newcoal || c[1] == newcoal) {
            // we either coal above the newcoal or our existing
            // node just broke and newcoal was underneath.
                
            // if newcoal was previously above spr->coal_node
            // then we rename the spr coal node
            if (last_nodes[spr->coal_node].parent == newcoal)
                spr->coal_node = newcoal;
        }
    }
            
    // determine if mapping of new node needs to be changed
    // newcoal was parent of recomb, it is broken
    if (last_nodes[spr->recomb_node].parent == newcoal) {
        mapping[newcoal] = -1;
        int p = last_nodes[newcoal].parent;
        if (p != -1)
            mapping[p] = newcoal;
    } else {
        // newcoal was not broken
        // find child without recomb or coal on it
        int x = newcoal;
        while (true) {
            int y = last_nodes[x].child[0];
            if (y == spr->coal_node || y == spr->recomb_node)
                y = last_nodes[x].child[1];
            x = y;
            if (mapping[x] != -1)
                break;
        }
        mapping[newcoal] = nodes[mapping[x]].parent;
    }
}



// add a thread to an ARG
void add_arg_thread(LocalTrees *trees, int ntimes, int *thread_path, int seqid,
                    vector<int> &recomb_pos, vector<NodePoint> &recombs)
{
    unsigned int irecomb = 0;
    int nleaves = trees->get_num_leaves();
    int nnodes = trees->nnodes;
    int nnodes2 = nnodes + 2;
    
    // node names
    int newleaf = nleaves;
    int displaced = nnodes;
    int newcoal = nnodes + 1;

    States states;
    State last_state;
    LocalTree *last_tree = NULL;


    // update trees info
    trees->seqids.push_back(seqid);
    trees->nnodes = nnodes2;
    

    // loop through blocks
    int end = trees->start_coord;
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it) {
        LocalTree *tree = it->tree;
        Spr *spr = &(it->spr);
        int start = end;
        end += it->blocklen;
        get_coal_states(tree, ntimes, states);
        
        // add new branch to local tree
        it->ensure_capacity(nnodes2);
        State state = states[thread_path[start]];
        add_tree_branch(tree, state.node, state.time);
        
        // update mapping and spr
        int *mapping = it->mapping;
        if (mapping) {
            //printf("spr %d %d, %d %d\n", 
            //       it->spr.recomb_node, it->spr.recomb_time,
            //       it->spr.coal_node, it->spr.coal_time);
            add_spr_branch(tree, last_tree, state, last_state,
                           &it->spr, mapping,
                           newleaf, displaced, newcoal);
            assert(assert_spr(last_tree, tree, spr, mapping));
        }

        // assert new branch is where it should be
        assert(tree->nodes[newcoal].age == states[thread_path[start]].time);


        // break this block for each new recomb within this block
        for (;irecomb < recombs.size() && 
              recomb_pos[irecomb] < end; irecomb++) {
            int pos = recomb_pos[irecomb];
            LocalNode *nodes = tree->nodes;
            state = states[thread_path[pos]];
            last_state = states[thread_path[pos-1]];
            
            // assert that thread time is still on track
            assert(tree->nodes[newcoal].age == last_state.time);

            // determine real name of recomb node
            // it may be different due to adding a new branch
            Spr spr2;
            spr2.recomb_node = recombs[irecomb].node;
            spr2.recomb_time = recombs[irecomb].time;
            if (spr2.recomb_node == newleaf)
                spr2.recomb_node = displaced;            
            assert(spr2.recomb_time <= tree->nodes[newcoal].age);

            // determine coal node and time
            //int istate = thread_path[pos];
            if (spr2.recomb_node == -1) {
                // recomb on new branch, coal given thread
                spr2.recomb_node = newleaf;
                spr2.coal_node = state.node;

                // fix coal node due to displacement
                if (spr2.coal_node == newleaf)
                    spr2.coal_node = displaced;

                // rename coal node due to newcoal underneath
                if (state.node == last_state.node &&
                    state.time > last_state.time)
                    spr2.coal_node = newcoal;

            } else {
                // recomb in ARG, coal on new branch
                if (state.time > last_state.time)
                    spr2.coal_node = nodes[newleaf].parent;
                else
                    spr2.coal_node = newleaf;
            }
            spr2.coal_time = state.time;

            
            // determine mapping:
            // all nodes keep their name expect the broken node, which is the
            // parent of recomb
            int *mapping2 = new int [tree->capacity];
            for (int j=0; j<nnodes2; j++)
                mapping2[j] = j;
            mapping2[nodes[spr2.recomb_node].parent] = -1;

            
            // make new local tree and apply SPR operation
            LocalTree *new_tree = new LocalTree(nnodes2, tree->capacity);
            new_tree->copy(*tree);
            apply_spr(new_tree, spr2);

            // calculate block end
            int block_end;
            if (irecomb < recombs.size() - 1)
                // use next recomb in this block to determine block end
                block_end = min(recomb_pos[irecomb+1], end);
            else
                // no more recombs in this block
                block_end = end;
           
            // insert new tree and spr into local trees list
            it->blocklen = pos - start;
            ++it;
            it = trees->trees.insert(it, 
                LocalTreeSpr(new_tree, spr2, block_end - pos, mapping2));
            

            // assert tree and SPR
            assert(assert_tree(new_tree));
            assert(new_tree->nodes[newcoal].age == state.time);
            assert(assert_spr(tree, new_tree, &spr2, mapping2));

            // remember the previous tree for next iteration of loop
            tree = new_tree;
            nodes = tree->nodes;
            start = pos;
        }

        // remember the previous tree for next iteration of loop
        last_tree = tree;
        last_state = states[thread_path[end-1]];
        if (last_state.node == newleaf)
            last_state.node = displaced;
    }

    assert_trees(trees);
}



// Removes a thread from an ARG
// NOTE: if remove_leaf is not last_leaf, nleaves - 1, 
// last_leaf is renamed to remove_leaf
void remove_arg_thread(LocalTrees *trees, int remove_seqid)
{
    int nnodes = trees->nnodes;
    int nleaves = trees->get_num_leaves();
    int displace[nnodes];
    int last_leaf = nleaves - 1;

    // find leaf to remove from seqid
    int remove_leaf = -1;
    for (unsigned int i=0; i<trees->seqids.size(); i++) {
        if (trees->seqids[i] == remove_seqid) {
            remove_leaf = i;
            break;
        }
    }
    assert(remove_leaf != -1);
    
    // special case for trunk genealogy
    if (nnodes == 3) {
        assert(remove_leaf == 0 || remove_leaf == 1);
        trees->make_trunk(trees->start_coord, trees->end_coord,
                          trees->begin()->tree->capacity);
        trees->seqids[0] = trees->seqids[1-remove_leaf];
        trees->seqids.resize(1);
        return;
    }

    
    // remove extra branch
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it) {
        LocalTree *tree = it->tree;
        LocalNode *nodes = tree->nodes;
        
        // get information about removal
        int remove_coal = nodes[remove_leaf].parent;
        int coal_time = nodes[remove_coal].age;
        int *c = nodes[remove_coal].child;
        int coal_child = (c[0] == remove_leaf ? c[1] : c[0]);
        
        // remove branch from tree        
        remove_tree_branch(tree, remove_leaf, displace);
        //assert_tree(tree);
        
        
        // fix this mapping due to displacement
        int *mapping = it->mapping;
        if (mapping) {
            for (int i=0; i<nnodes-2; i++)
                if (mapping[i] != -1)
                    mapping[i] = displace[mapping[i]];
                else
                    mapping[i] = -1;
        }
        
        // get next tree
        LocalTrees::iterator it2 = it;
        ++it2;
        if (it2 == trees->end())
            continue;

        // fix next mapping due to displacement
        mapping = it2->mapping;
        if (displace[last_leaf] != -1)
            mapping[displace[last_leaf]] = mapping[last_leaf];
        if (displace[nnodes-2] != -1)
            mapping[displace[nnodes-2]] = mapping[nnodes-2];
        if (displace[nnodes-1] != -1)
            mapping[displace[nnodes-1]] = mapping[nnodes-1];
        
        
        
        // fix SPR
        Spr *spr = &it2->spr;
        
        // get new name of coal_child
        coal_child = displace[coal_child];

        // if recomb is on branch removed, prune it
        if (spr->recomb_node == remove_leaf) {
            spr->set_null();
            continue;
        }

        // see if recomb node is renamed
        if (spr->recomb_node == remove_coal) {
            spr->recomb_node = coal_child;
        } else {
            // rename recomb_node due to displacement
            spr->recomb_node = displace[spr->recomb_node];
        }
        
        // if recomb is on root branch, prune it
        if (spr->recomb_node == coal_child && nodes[coal_child].parent == -1) {
            spr->set_null();
            continue;
        }

        // rename spr coal_node
        if (spr->coal_node == remove_leaf) {
            // mediated coal
            spr->coal_node = coal_child;
            spr->coal_time = coal_time;

        } else if (spr->coal_node == remove_coal) {
            // move coal down a branch
            spr->coal_node = coal_child;
        } else {
            // rename recomb_node due to displacement
            spr->coal_node = displace[spr->coal_node];
        }


        // check for bubbles
        if (spr->recomb_node == spr->coal_node) {
            spr->set_null();
            continue;
        }
    }
    

    // update trees info
    trees->seqids[remove_leaf] = trees->seqids[last_leaf];
    trees->seqids.resize(nleaves - 1);
    trees->nnodes -= 2;
    
    // remove extra trees
    remove_null_sprs(trees);
    
    assert_trees(trees);
}


//=============================================================================
// internal branch threading operations


// find recoal node, it is the node with no inward mappings
int get_recoal_node(const LocalTree *tree, 
                    const Spr &spr, const int *mapping)
{
    const int nnodes = tree->nnodes;
    bool mapped[nnodes];
    fill(mapped, mapped + nnodes, false);

    for (int i=0; i<nnodes; i++)
        if (mapping[i] != -1)
            mapped[mapping[i]] = true;
    
    for (int i=0; i<nnodes; i++)
        if (!mapped[i])
            return i;

    assert(false);
    return -1;
}


// find the next possible branches in a removal path
void get_next_removal_nodes(const LocalTree *tree1, 
                            const Spr &spr2, const int *mapping2,
                            int node, int next_nodes[2])
{   
    const int recoal = get_recoal_node(tree1, spr2, mapping2);

    // get passive transition
    next_nodes[0] = mapping2[node];
    if (next_nodes[0] == -1) {
        // node is broken by SPR
        // next node is then non-recomb child or recoal
        next_nodes[0] = tree1->get_sibling(spr2.recomb_node);
        if (spr2.coal_node == next_nodes[0])
            next_nodes[0] = recoal;
    }
    
    // get possible active transition
    // if recoal is on this branch (node) then there is a split in the path
    if (spr2.coal_node == node) {
        // find recoal node, its the node with no inward mappings
        next_nodes[1] = recoal;
    } else {
        // no second transition
        next_nodes[1] = -1;
    }
}


void sample_arg_removal_path(LocalTrees *trees, int node, int *path)
{
    int i = 0;
    path[i++] = node;
    LocalTree *last_tree = NULL;

    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it) {
        LocalTree *tree = it->tree;

        if (last_tree) {
            int next_nodes[2];
            get_next_removal_nodes(last_tree, it->spr, it->mapping,
                                   path[i-1], next_nodes);
            int j = (next_nodes[1] != -1 ? irand(2) : 0);
            path[i++] = next_nodes[j];
        }
        
        last_tree = tree;
    }
}


// Removes a thread path from an ARG and returns a partial ARG
void remove_arg_thread_path(LocalTrees *trees, const int *removal_path, 
                            int maxtime)
{
    LocalTree *tree = NULL;
    
    int i=0;
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it, i++) 
    {
        LocalTree *last_tree = tree;
        tree = it->tree;
        LocalNode *nodes = tree->nodes;

        int removal_node = removal_path[i];
        
        // modify local into subtree-maintree format
        int broken_node = nodes[removal_node].parent;
        int coal_time = nodes[broken_node].age;
        int broken_child = tree->get_sibling(removal_node);
        Spr removal_spr(removal_node, nodes[removal_node].age,
                        tree->root, maxtime);

        if (removal_node == tree->root) {
            // removal path has "fallen off the top" there is nothing to edit
            continue;
        } else 
            apply_spr(tree, removal_spr);
        
        // determine subtree and maintree roots
        int subtree_root = removal_node;
        int maintree_root = tree->get_sibling(subtree_root);
        
        // ensure subtree is the first child of the root
        int *c = nodes[tree->root].child;
        if (c[0] == maintree_root) {
            c[0] = subtree_root;
            c[1] = maintree_root;
        }
        
        // fix previous mapping
        if (it->mapping) {
            assert(last_tree);
            if (removal_path[i-1] != last_tree->root)
                it->mapping[last_tree->root] = tree->root;
        }


        // get next tree
        LocalTrees::iterator it2 = it;
        ++it2;
        if (it2 == trees->end())
            continue;



        // fix SPR
        Spr *spr = &it2->spr;
        int *mapping = it2->mapping;

        // if recomb is on branch removed, prune it
        if (spr->recomb_node == removal_node) {
            int p = nodes[spr->recomb_node].parent;
            assert(mapping[p] != -1 || p == tree->root);
            spr->set_null();
            continue;
        }

        // see if recomb node is renamed
        if (spr->recomb_node == broken_node) {
            spr->recomb_node = broken_child;
        }
        
        // detect branch path splits
        int next_nodes[2];
        get_next_removal_nodes(tree, *spr, mapping, 
                               removal_path[i], next_nodes);

        if (spr->coal_node == removal_node) {
            if (removal_path[i+1] == next_nodes[0]) {
                // removal path chooses lower path
                // note: sister_node = broken_child
                
                if (spr->recomb_node == broken_child) {
                    // spr is now bubble, prune it
                    int p = nodes[spr->recomb_node].parent;
                    assert(mapping[p] != -1 || p == tree->root);
                    spr->set_null();
                    continue;
                } else {
                    // recomb is on non-sister branch, therefore it is a
                    // mediated coalescence
                    spr->coal_node = broken_child;
                    spr->coal_time = coal_time;
                }
            } else {
                // removal path chooses upper path
                // keep spr recoal where it is

                // assert that upper path is the new recoal node
                // nobody should map to the new recoal node
                for (int j=0; j<tree->nnodes; j++)
                    assert(mapping[j] != removal_path[i+1]);
            }
        } else if (spr->coal_node == broken_node) {
            // rename spr recoal
            spr->coal_node = broken_child;
        }
        
        // check for bubbles
        if (spr->recomb_node == spr->coal_node) {
            int p = nodes[spr->recomb_node].parent;
            assert(mapping[p] != -1 || p == tree->root);
            spr->set_null();
            continue;
        }

        // ensure broken node maps to -1
        int spr_broken_node = nodes[spr->recomb_node].parent;
        mapping[spr_broken_node] = -1;

        // assert spr
        if (last_tree && !it->spr.is_null())
            assert_spr(last_tree, tree, &it->spr, it->mapping);
    }

    assert_trees(trees);
    
    // remove extra trees
    remove_null_sprs(trees);
    
    assert_trees(trees);
}



//=============================================================================
// C interface

extern "C" {

void arghmm_sample_arg_removal_path(LocalTrees *trees, int node, int *path)
{
    sample_arg_removal_path(trees, node, path);
}

void arghmm_remove_arg_thread_path(LocalTrees *trees, int *removal_path, 
                                   int maxtime)
{
    remove_arg_thread_path(trees, removal_path, maxtime);
}


} // extern C

} // namespace arghmm

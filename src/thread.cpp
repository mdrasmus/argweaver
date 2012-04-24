
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
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it) {
        int start = it->block.start;
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
    }

    return true;
}


void displace_node(LocalTree *tree, int src_node, int dest_node)
{
    // special case
    if (src_node == dest_node)
        return;

    LocalNode *nodes = tree->nodes;

    // copy displaced node
    nodes[dest_node] = nodes[src_node];

    if (nodes[dest_node].parent != -1) {
        int *c = nodes[nodes[dest_node].parent].child;
        if (c[0] == src_node)
            c[0] = dest_node;
        else
            c[1] = dest_node;
    }
    int *c = nodes[dest_node].child;
    nodes[c[0]].parent = dest_node;
    nodes[c[1]].parent = dest_node;
}


void add_tree_branch(LocalTree *tree, State state)
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
    int node2 = (state.node != newleaf ? state.node : displaced);
    int parent = nodes[state.node].parent;
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
    nodes[newcoal].age = state.time;

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
    // assert new branch is where it should be
    assert(tree->nodes[newcoal].age == state.time);   
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
    // however, if it equals newcoal, then recomb branch is 
    // either renamed or we have mediation
    int recoal = nodes[mapping[spr->recomb_node]].parent;
    if (recoal == newcoal) {
        if (mapping[last_state.node] == node2) {
            // recomb is above coal state, we rename spr recomb node
            spr->recomb_node = newcoal;
        } else {
            // this is a mediated coal, rename coal node and time
            spr->coal_node = newleaf;
            spr->coal_time = state.time;
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
            //printf("x=%d, y=%d\n", x, y);
            if (y == spr->coal_node || y == spr->recomb_node)
                y = last_nodes[x].child[1];
            x = y;
            //printf("  x=%d, mapping[x] = %d\n", x, mapping[x]);
            if (mapping[x] != -1)
                break;
        }
        //printf(":: %d -> %d, x=%d\n", 
        //       newcoal, nodes[mapping[x]].parent, x);
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
    
    int newleaf = nleaves;
    int displaced = nnodes;
    int newcoal = nnodes + 1;

    States states;
    State last_state;
    LocalTree *last_tree = NULL;


    // update trees info
    trees->seqids.push_back(seqid);
    //for (int i=0; i<trees->seqids.size(); i++)
    //    printf("seqids[%d] = %d %d\n", i, trees->seqids[i], nleaves);


    // loop through blocks
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it) {
        LocalTree *tree = it->tree;
        Spr *spr = &(it->spr);
        const int start = it->block.start;
        const int end = it->block.end;
        get_coal_states(tree, ntimes, states);

        // DBEUG
        //printf("add %d\n", start);        
        
        // add new branch to local tree
        it->ensure_capacity(nnodes2);
        State state = states[thread_path[start]];
        add_tree_branch(tree, state);
        
        // update mapping and spr
        int *mapping = it->mapping;
        if (mapping) {
            add_spr_branch(tree, last_tree, state, last_state,
                           &it->spr, mapping,
                           newleaf, displaced, newcoal);
            // assert SPR
            if (!assert_spr(last_tree, tree, spr, mapping))
                printf("!!! %d spr fail\n", start);
        }

        // assert new branch is where it should be
        assert(tree->nodes[newcoal].age == states[thread_path[start]].time);


        // break this block for each new recomb within this block
        for (;irecomb < recombs.size() && 
              recomb_pos[irecomb] < end; irecomb++) {
            int pos = recomb_pos[irecomb];
            //printf("start %d pos %d\n", start, pos);

            LocalNode *nodes = tree->nodes;

            assert(tree->nodes[newcoal].age == states[thread_path[pos-1]].time);

            // determine real name of recomb node
            // it may be different due to displacement
            Spr spr2;
            spr2.recomb_node = recombs[irecomb].node;
            spr2.recomb_time = recombs[irecomb].time;
            if (spr2.recomb_node == newleaf)
                spr2.recomb_node = displaced;            
            assert(spr2.recomb_time <= tree->nodes[newcoal].age);

            // determine coal node and time
            int istate = thread_path[pos];
            if (spr2.recomb_node == -1) {
                // recomb on new branch, coal given thread
                spr2.recomb_node = newleaf;
                spr2.coal_node = states[istate].node;

                // fix coal node due to displacement
                if (spr2.coal_node == newleaf)
                    spr2.coal_node = displaced;

                // rename due to newcoal
                if (states[istate].node == states[thread_path[pos-1]].node &&
                    states[istate].time > states[thread_path[pos-1]].time)
                    spr2.coal_node = newcoal;

            } else {
                // recomb in ARG, coal on new branch
                // fix recomb node due to displacement
                if (spr2.recomb_node == newleaf)
                    spr2.recomb_node = displaced;

                if (states[istate].time > states[thread_path[pos-1]].time)
                    spr2.coal_node = nodes[newleaf].parent;
                else
                    spr2.coal_node = newleaf;
            }
            spr2.coal_time = states[istate].time;

            
            // determine mapping
            int *mapping2 = new int [tree->capacity];
            for (int j=0; j<nnodes2; j++)
                mapping2[j] = j;
            mapping2[nodes[spr2.recomb_node].parent] = -1;

            
            // make new local tree and apply SPR operation
            LocalTree *new_tree = new LocalTree(nnodes2, tree->capacity);
            LocalNode *new_nodes = new_tree->nodes;
            for (int j=0; j<nnodes2; j++) {
                new_nodes[j].parent = nodes[j].parent;
                new_nodes[j].child[0] = nodes[j].child[0];
                new_nodes[j].child[1] = nodes[j].child[1];
                new_nodes[j].age = nodes[j].age;
            }
            apply_spr(new_tree, &spr2);

            // calculate block end
            int block_end;
            if (irecomb < recombs.size() - 1)
                block_end = min(recomb_pos[irecomb+1], end);
            else
                block_end = end;
           
            // insert new tree into local trees list
            it->block.end = pos;
            ++it;
            it = trees->trees.insert(it, 
                LocalTreeSpr(pos, block_end, new_tree, spr2, mapping2));


            // assert tree and SPR
            assert(assert_tree(new_tree));
            assert(new_tree->nodes[newcoal].age == 
                   states[thread_path[pos]].time);
            if (!assert_spr(tree, new_tree, &spr2, mapping2))
                printf("!!! %d spr fail\n", pos);

            // remember the previous tree
            tree = new_tree;
            nodes = tree->nodes;
        }

        last_tree = tree;
        last_state = states[thread_path[end-1]];
        if (last_state.node == newleaf)
            last_state.node = displaced;
    }

    
    // update number of nodes
    trees->nnodes = nnodes2;

    //assert_trees(trees);
}


// Removes a thread from an ARG
// NOTE: if remove_leaf is not last_leaf, nleaves - 1, 
// last_leaf is renamed to remove_leaf
void remove_arg_thread(LocalTrees *trees, int remove_leaf)
{
    int nnodes = trees->nnodes;
    int nleaves = trees->get_num_leaves();
    int displace[nnodes];
    int last_leaf = nleaves - 1;


    // TODO: remove leaf must be looked up in seqids
    
    // special case for trunk genealogy
    if (nnodes == 3) {
        assert(remove_leaf == 0 || remove_leaf == 1);
        trees->make_trunk(trees->start_coord, trees->end_coord,
                          trees->begin()->tree->capacity);
        trees->seqids[0] = trees->seqids[1-remove_leaf];
        trees->seqids.resize(1);
        return;
    }

    // update trees info
    trees->seqids[remove_leaf] = trees->seqids[last_leaf];
    trees->seqids.resize(nleaves - 1);
    
    
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it) {
        LocalTree *tree = it->tree;
        LocalNode *nodes = tree->nodes;
        
        // remove coal node
        int remove_coal = nodes[remove_leaf].parent;
        int coal_time = nodes[remove_coal].age;
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
        if (remove_coal != nnodes-2) {
            displace[nnodes-2] = last_leaf;
            displace_node(tree, nnodes-2, last_leaf);
        }
        if (remove_coal != nnodes-1) {
            displace[nnodes-1] = remove_coal;
            displace_node(tree, nnodes-1, remove_coal);
        }

        // set tree data
        tree->nnodes -= 2;
        tree->set_root();

        
        // get new name of coal_child
        coal_child = displace[coal_child];


        // get next tree
        LocalTrees::iterator it2 = it;
        ++it2;
        if (it2 == trees->end())
            continue;
        
        // fix SPR
        Spr *spr = &it2->spr;

        // if recomb is on branch removed, prune it
        if (spr->recomb_node == remove_leaf) {
            it2->clear();
            trees->trees.erase(it2);
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
            it2->clear();
            trees->trees.erase(it2);
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
            it2->clear();
            trees->trees.erase(it2);
            continue;
        }


        // fix this mapping due to displacement
        int *mapping = it->mapping;
        if (mapping) {
            for (int i=0; i<nnodes-2; i++)
                mapping[i] = displace[mapping[i]];
        }

        // fix next mapping due to displacement
        mapping = it2->mapping;
        mapping[displace[last_leaf]] = mapping[last_leaf];
        mapping[displace[nnodes-2]] = mapping[nnodes-2];
        mapping[displace[nnodes-1]] = mapping[nnodes-1];
    }
}


} // namespace arghmm

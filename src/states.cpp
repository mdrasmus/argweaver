

#include "states.h"


namespace arghmm {



// Converts integer-based states to State class
void make_states(intstate *istates, int nstates, States &states) {
    states.clear();
    for (int i=0; i<nstates; i++)
        states.push_back(State(istates[i][0], istates[i][1]));
}


// Converts state class represent to integer-based
void make_intstates(States states, intstate *istates)
{
    const int nstates = states.size();
    for (int i=0; i<nstates; i++) {
        istates[i][0] = states[i].node;
        istates[i][1] = states[i].time;
    }
}



void get_coal_states(const LocalTree *tree, int ntimes, States &states,
                     bool internal)
{
    if (internal)
        get_coal_states_internal(tree, ntimes, states);
    else
        get_coal_states_external(tree, ntimes, states);
}

int get_num_coal_states(const LocalTree *tree, int ntimes, bool internal)
{
    if (internal)
        return get_num_coal_states_internal(tree, ntimes);
    else
        return get_num_coal_states_external(tree, ntimes);
}


// Returns the possible coalescing states for a tree
//
// NOTE: Do not allow coalescing at top time
// states for the same branch are clustered together and ages are given
// in increasing order
void get_coal_states_external(const LocalTree *tree, int ntimes, States &states)
{
    states.clear();
    const LocalNode *nodes = tree->nodes;
    
    // iterate over the branches of the tree
    for (int i=0; i<tree->nnodes; i++) {
        int time = nodes[i].age;
        const int parent = nodes[i].parent;
        
        if (parent == -1) {
            // no parent, allow coalescing up basal branch until ntimes-2
            for (; time<ntimes-1; time++)
                states.push_back(State(i, time));
        } else {
            // allow coalescing up branch until parent
            const int parent_age = nodes[parent].age;
            for (; time<=parent_age; time++)
                states.push_back(State(i, time));
        }
    }
}

// Returns the number of possible coalescing states for a tree
int get_num_coal_states_external(const LocalTree *tree, int ntimes)
{
    int nstates = 0;
    const LocalNode *nodes = tree->nodes;
    
    // iterate over the branches of the tree
    for (int i=0; i<tree->nnodes; i++) {
        int time = nodes[i].age;
        const int parent = nodes[i].parent;
        
        if (parent == -1) {
            // no parent, allow coalescing up basal branch until ntimes-2
            for (; time<ntimes-1; time++)
                nstates++;
        } else {
            // allow coalescing up branch until parent
            const int parent_age = nodes[parent].age;
            for (; time<=parent_age; time++)
                nstates++;
        }
    }

    return nstates;
}


// Returns the possible coalescing states for a tree that is in
// subtree-maintree format for internal branch resampling
//
// NOTE: Do not allow coalescing at top time
// states for the same branch are clustered together and ages are given
// in increasing order
void get_coal_states_internal(const LocalTree *tree, int ntimes, States &states)
{
    states.clear();
    const int nnodes = tree->nnodes;
    const LocalNode *nodes = tree->nodes;

    if (nodes[tree->root].age < ntimes) {
        // we have a fully specified tree, there are no states
        return;
    }

    int subtree_root = nodes[tree->root].child[0];
    int minage = nodes[subtree_root].age;

    // ignore root and whole subtree
    bool ignore[nnodes];
    fill(ignore, ignore + nnodes, false);
    ignore[tree->root] = true;

    // recurse down subtree
    int stack[nnodes];
    int stacki = 0;
    stack[stacki++] = subtree_root;
    while (stacki > 0) {
        // pop stack
        int node = stack[--stacki];
        ignore[node] = true;

        // push children on stack
        if (!nodes[node].is_leaf()) {
            stack[stacki++] = nodes[node].child[0];
            stack[stacki++] = nodes[node].child[1];
        }
    }
    
    
    // iterate over the branches of the tree
    for (int i=0; i<nnodes; i++) {
        int time = max(nodes[i].age, minage);
        const int parent = nodes[i].parent;
        
        // skip subtree and root node
        if (ignore[i])
            continue;

        if (parent == tree->root) {
            // no parent, allow coalescing up basal branch until ntimes-2
            for (; time<ntimes-1; time++)
                states.push_back(State(i, time));
        } else {
            // allow coalescing up branch until parent
            const int parent_age = nodes[parent].age;
            for (; time<=parent_age; time++)
                states.push_back(State(i, time));
        }
    }
}


// Returns the number of possible coalescing states for a tree
int get_num_coal_states_internal(const LocalTree *tree, int ntimes)
{
    int nstates = 0;
    const LocalNode *nodes = tree->nodes;
    int subtree_root = nodes[tree->root].child[0];
    int minage = nodes[subtree_root].age;
    
    // iterate over the branches of the tree
    for (int i=0; i<tree->nnodes; i++) {
        int time = max(nodes[i].age, minage);
        const int parent = nodes[i].parent;
        
        // skip subtree root branch and root node
        if (i == subtree_root || i == tree->root)
            continue;

        if (parent == tree->root) {
            // no parent, allow coalescing up basal branch until ntimes-2
            for (; time<ntimes-1; time++)
                nstates++;
        } else {
            // allow coalescing up branch until parent
            const int parent_age = nodes[parent].age;
            for (; time<=parent_age; time++)
                nstates++;
        }
    }

    return nstates;
}



//=============================================================================
// C interface

extern "C" {


void arghmm_get_nstates(LocalTrees *trees, int ntimes, bool internal,
                        int *nstates)
{
    States states;
    
    // iterate over local trees
    int i = 0;
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); it++) {
        LocalTree *tree = it->tree;

        get_coal_states(tree, ntimes, states, internal);
        int n = states.size();
        for (int j=0; j<it->blocklen; j++)
            nstates[i+j] = n;
        i += it->blocklen;
    }
}


// Returns state-spaces, useful for calling from python
intstate **get_state_spaces(LocalTrees *trees, int ntimes, bool internal)
{
    States states;
    
    // allocate state space
    intstate **all_states = new intstate* [trees->get_num_trees()];

    // iterate over local trees
    int i = 0;
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); it++) {
        LocalTree *tree = it->tree;
        get_coal_states(tree, ntimes, states);
        int nstates = states.size();
        all_states[i] = new intstate [nstates];
        
        for (int j=0; j<nstates; j++) {
            all_states[i][j][0] = states[j].node;
            all_states[i][j][1] = states[j].time;
        }
        i++;
    }

    return all_states;
}


// Deallocate state space memory
void delete_state_spaces(intstate **all_states, int ntrees)
{
    for (int i=0; i<ntrees; i++)
        delete [] all_states[i];
    delete [] all_states;
}


} // extern "C"

} // namespace arghmm

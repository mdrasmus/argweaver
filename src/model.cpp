

#include "model.h"


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



// Retruns the possible coalescing states for a tree
//
// NOTE: Do not allow coalescing at top time
void get_coal_states(LocalTree *tree, int ntimes, States &states)
{
    states.clear();
    LocalNode *nodes = tree->nodes;
    
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


} // namespace arghmm

//=============================================================================
// ArgHmm states


#ifndef ARGHMM_STATES_H
#define ARGHMM_STATES_H

#include "common.h"
#include "local_tree.h"

namespace arghmm {


// A state in the ArgHmm
//
// Each state represents a node and time where coalescing is allowed
class State
{
public:
    State(int node=0, int time=0) :
        node(node), time(time) {}

    inline bool operator==(const State &other) const {
        return (node == other.node) && (time == other.time);
    }

    int node;
    int time;
};

// A state space for a local block
typedef vector<State> States;


// This data structure provides a mapping from (node, time) tuples to
// the corresponding state index.
class NodeStateLookup
{
public:
    NodeStateLookup(const States &states, int nnodes) :
        nstates(states.size()),
        nnodes(nnodes)
    {
        const int nstates = states.size();

        // allocate lookup arrays
        node_offset = new int[nnodes];
        state_lookup = new int[nstates];


        // count number of states per node and mintime per node
        int nstates_per_node[nnodes];
        int node_mintimes[nnodes];

        // initialize arrays
        for (int i=0; i<nnodes; i++) {
            nstates_per_node[i] = 0;
            node_mintimes[i] = nstates;
        }

        for (int i=0; i<nstates; i++) {
            nstates_per_node[states[i].node]++;
            node_mintimes[states[i].node] = min(node_mintimes[states[i].node], 
                                                states[i].time);
        }

        // setup node_offsets
        int offset = 0;
        for (int i=0; i<nnodes; i++) {
            node_offset[i] = offset - node_mintimes[i];
            offset += nstates_per_node[i];
        }

        // set states
        for (int i=0; i<nstates; i++)
            state_lookup[node_offset[states[i].node] + states[i].time] = i;
    }

    ~NodeStateLookup()
    {
        // clean up lookup arrays
        delete [] node_offset;
        delete [] state_lookup;
    }

    // Returns the state index for state (node, time)
    inline int lookup(int node, int time) {
        return state_lookup[node_offset[node] + time];
    }

    int nstates;
    int nnodes;
    int *node_offset;
    int *state_lookup;
};


// A simple representation of a state, useful for passing from python
typedef int intstate[2];


// Converts integer-based states to State class
void make_states(intstate *istates, int nstates, States &states);

// Converts state class represent to integer-based
void make_intstates(States states, intstate *istates);

void get_coal_states(const LocalTree *tree, int ntimes, States &states);
int get_num_coal_states(const LocalTree *tree, int ntimes);

void get_coal_states_internal(const LocalTree *tree, int ntimes,States &states);
int get_num_coal_states_internal(const LocalTree *tree, int ntimes);



} // namespace arghmm


#endif // ARGHMM_STATES_H

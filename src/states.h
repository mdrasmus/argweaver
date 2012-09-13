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

    void set(const int &_node, const int &_time) {
        node = _node;
        time = _time;
    }

    void set_null()
    {
        node = -1;
        time = -1;
    }

    bool is_null() const 
    {
        return node == -1 && time -1;
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
        states(states),
        nstates(states.size()),
        nnodes(nnodes)
    {
        const int nstates = states.size();
        const int MAXTIME = 1000000;

        // allocate lookup arrays
        node_offset = new int[nnodes];
        state_lookup = new int[nstates];
        nstates_per_node = new int[nnodes];
        
        // count number of states per node and mintime per node
        int node_mintimes[nnodes];

        // initialize arrays
        for (int i=0; i<nnodes; i++) {
            nstates_per_node[i] = 0;
            node_mintimes[i] = MAXTIME;
        }

        for (int i=0; i<nstates; i++) {
            state_lookup[i] = -1;
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
        for (int i=0; i<nstates; i++) {
            int j = node_offset[states[i].node] + states[i].time;
            assert(j >=0 && j<nstates);
            state_lookup[j] = i;
        }
    }

    ~NodeStateLookup()
    {
        // clean up lookup arrays
        delete [] node_offset;
        delete [] state_lookup;
        delete [] nstates_per_node;
    }

    // Returns the state index for state (node, time)
    inline int lookup(int node, int time) {
        if (nstates_per_node[node] == 0)
            return -1;
        const int i = node_offset[node] + time;
        if (i < 0 || i >= nstates)
            return -1;
        const int statei = state_lookup[i];
        if (statei == -1)
            return -1;
        if (states[statei].node != node ||
            states[statei].time != time)
            return -1;
        return statei;
    }

protected:
    const States &states;
    int nstates;
    int nnodes;
    int *node_offset;
    int *state_lookup;
    int *nstates_per_node;
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

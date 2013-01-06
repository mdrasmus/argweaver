//=============================================================================
// transitions

#ifndef ARGHMM_TRANS_H
#define ARGHMM_TRANS_H

#include "common.h"
#include "local_tree.h"
#include "model.h"
#include "states.h"

namespace arghmm {


// A compressed representation of the transition matrix
class TransMatrix
{
public:
    TransMatrix(int ntimes, int nstates, bool alloc=true) :
        ntimes(ntimes),
        nstates(nstates),
        own_data(false),
        internal(false)
    {
        if (alloc)
            allocate(ntimes);
    }

    ~TransMatrix()
    {
        if (own_data) {
            delete [] D;
            delete [] E;
            delete [] B_alloc;
            delete [] E2;
            delete [] G1;
            delete [] G2;
            delete [] G3;
            delete [] G4;
            delete [] norecombs;
        }
    }

    // allocate space for transition matrix
    void allocate(int ntimes)
    {
        ntimes = ntimes;
        own_data = true;
        D = new double [ntimes];
        E = new double [ntimes];
        B_alloc = new double [ntimes+1];
        B = &B_alloc[1];
        E2 = new double [ntimes];
        G1 = new double [ntimes];
        G2 = new double [ntimes];
        G3 = new double [ntimes];
        G4 = new double [ntimes];
        norecombs = new double [ntimes];
    }
    
    // Returns the probability of transition from state i to j
    inline double get(
        const LocalTree *tree, const States &states, int i, int j) const
    {
        int minage = 0;
        if (internal) {
            if (nstates == 0)
                return 1.0;
            const int subtree_root = tree->nodes[tree->root].child[0];
            const int subtree_age = tree->nodes[subtree_root].age;
            minage = subtree_age;
        }

        const int node1 = states[i].node;
        const int a = states[i].time;
        const int node2 = states[j].node;
        const int b = states[j].time;
        const int c = tree->nodes[node2].age;

        return get_time(a, b, c, minage, node1 == node2);
    }

    // Returns the probability of transition from state1 with time 'a'
    // to state2 with time 'b'.  The probability also depends on whether
    // the node changes between states ('same_node') or whether there is
    // a minimum age ('minage') allowed for the state.
    inline double get_time(int a, int b, int c, 
                           int minage, bool same_node) const
    {
        if (a < minage || b < minage)
            return 0.0;
        
        if (!same_node) {
            return D[a] * E[b] * 
                (E2[b] * (B[min(a,b-1)] + int(a<b) * G1[a])
                 + int(a>b) * G2[b]
                 + int(a==b) * G3[b]
                 - (minage>0 ? G4[b]*B[minage-1] : 0));
        } else {
            double p = D[a] * E[b] * 
                ((2 * (E2[b] * (B[min(a,b-1)] + int(a<b) * G1[a])
                       + int(a>b) * G2[b]
                       + int(a==b) * G3[b]))
                 - (c>0 ? G4[b]*B[c-1] : 0)
                 - (minage>0 ? G4[b]*B[minage-1] : 0));
            if (a == b)
                p += norecombs[a];
            return p;
        }
    }

    // Returns the log probability of transition from state i to j
    inline double get_log(
        const LocalTree *tree, const States &states, int i, int j) const
    {
        return log(get(tree, states, i, j));
    }


    int ntimes;
    int nstates;
    bool own_data;
    double *D;
    double *E;
    double *B;
    double *B_alloc;
    double *E2;
    double *G1;
    double *G2;
    double *G3;
    double *G4;
    double *norecombs;
    bool internal;
};


// A compressed representation of the switch transition matrix
class TransMatrixSwitch
{
public:
    TransMatrixSwitch(int nstates1, int nstates2, bool alloc=true) :
        nstates1(nstates1),
        nstates2(nstates2),
        own_data(false)
    {
        if (alloc)
            allocate(nstates1, nstates2);
    }

    ~TransMatrixSwitch()
    {
        if (own_data) {
            delete [] determ;
            delete [] determprob;
            delete [] recoalrow;
            delete [] recombrow;
        }
    }

    void allocate(int nstates1, int nstates2)
    {
        nstates1 = nstates1;
        nstates2 = nstates2;
        
        // NOTE: nstates1 and nstates2 might be zero
        // we still calculate transitions for a state space of size zero
        
        own_data = true;
        determ = new int [max(nstates1, 1)];
        determprob = new double [max(nstates1, 1)];
        recoalrow = new double [max(nstates2, 1)];
        recombrow = new double [max(nstates2, 1)];
    }
    
    inline double get_log(int i, int j) const
    {
        return log(get(i, j));
    }
    
    inline double get(int i, int j) const
    {
        if (i == recoalsrc) {
            return recoalrow[j];
        } else if (i == recombsrc) {
            return recombrow[j];
        } else {
            if (determ[i] == j)
                return determprob[i];
            else
                return 0.0;
        }
    }
    

    int nstates1;
    int nstates2;
    int recoalsrc;
    int recombsrc;
    bool own_data;
    int *determ;
    double *determprob;
    double *recoalrow;
    double *recombrow;
};



//=============================================================================

void calc_transition_probs(const LocalTree *tree, const ArgModel *model,
    const States &states, const LineageCounts *lineages, TransMatrix *matrix,
    bool internal=false);
void calc_transition_probs(const LocalTree *tree, const ArgModel *model,
                          const States &states, const LineageCounts *lineages,
                          double **transprob);
void get_transition_probs(const LocalTree *tree, const ArgModel *model,
                           const States &states, const LineageCounts *lineages,
                           const TransMatrix *matrix, double **transprob);
void calc_transition_probs_internal(const LocalTree *tree, 
    const ArgModel *model, const States &states, const LineageCounts *lineages,
    TransMatrix *matrix);


void calc_transition_probs_switch(
    const LocalTree *tree, const LocalTree *last_tree, 
    const Spr &spr, const int *mapping,
    const States &states1, const States &states2,
    const ArgModel *model, const LineageCounts *lineages, 
    TransMatrixSwitch *transmat_switch, bool internal=false);
void calc_transition_probs_switch(
    const LocalTree *tree, const LocalTree *last_tree, 
    const Spr &spr, const int *mapping,
    const States &states1, const States &states2,
    const ArgModel *model, const LineageCounts *lineages, double **transprob);
void get_transition_probs_switch(const TransMatrixSwitch *matrix, 
                                 double **transprob);
void calc_transition_probs_switch_internal(
    const LocalTree *tree, const LocalTree *last_tree, 
    const Spr &spr, const int *mapping,
    const States &states1, const States &states2,
    const ArgModel *model, const LineageCounts *lineages, 
    TransMatrixSwitch *transmat_switch);
void calc_transition_probs_switch_internal2(
    const LocalTree *tree, const LocalTree *last_tree, 
    const Spr &spr, const int *mapping,
    const States &states1, const States &states2,
    const ArgModel *model, const LineageCounts *lineages, 
    TransMatrixSwitch *transmat_switch);


double calc_state_priors(
    int time, const int *nbranches, const int *ncoals, int minage,
    const double *popsizes, const double *coal_time_steps, int ntimes);
void calc_state_priors(const States &states, const LineageCounts *lineages, 
                       const ArgModel *model, double *priors,
                       const int minage=0);

void get_deterministic_transitions(
    const LocalTree *last_tree, const LocalTree *tree, 
    const Spr &spr, const int *mapping,
    const States &states1, const States &states2,
    NodeStateLookup &state2_lookup,
    int ntimes, int *next_states, bool internal=false);


bool assert_transmat(const LocalTree *tree, const ArgModel *model,
                     const TransMatrix *matrix);



//=============================================================================
// misc

void get_recomb_transition_switch(
    const LocalTree *tree, const LocalTree *last_tree, 
    const Spr &spr, const int *mapping,
    const States &states1, const States &states2, int next_states[2]);



} // namespace arghmm

#endif // ARGHMM_TRANS_H

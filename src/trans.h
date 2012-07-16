//=============================================================================
// transitions

#ifndef ARGHMM_TRANS_H
#define ARGHMM_TRANS_H

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
            delete [] B;
            delete [] D;
            delete [] E;
            delete [] G;
            delete [] norecombs;
            //delete [] sums;
        }
    }

    void allocate(int ntimes)
    {
        ntimes = ntimes;
        own_data = true;
        B = new double [ntimes];
        D = new double [ntimes];
        E = new double [ntimes];
        G = new double [ntimes];
        norecombs = new double [ntimes];
        //sums = new double [nstates];
    }


    inline double get(
        const LocalTree *tree, const States &states, int i, int j) const
    {
        double Bq = 0.0;
        int minage = 0;
        if (internal) {
            if (nstates == 0)
                return 0.0;
            const int subtree_root = tree->nodes[tree->root].child[0];
            const int subtree_age = tree->nodes[subtree_root].age;
            minage = subtree_age;
            if (subtree_age > 0)
                Bq = B[subtree_age - 1];
        }

        const int node1 = states[i].node;
        const int a = states[i].time;
        const int node2 = states[j].node;
        const int b = states[j].time;
        const int c = tree->nodes[node2].age;
        const double Bc = (c > 0 ? B[c-1] : 0.0);
        const double I = double(a <= b);

        if (a < minage || b < minage)
            return 0.0;

        if (node1 != node2)
            return D[a] * E[b] * (B[min(a,b)] - Bq - I * G[a]);
        else {
            double p = D[a] * E[b] * (2*B[min(a,b)] - Bq - 2*I*G[a] - Bc);
            if (a == b)
                p += norecombs[a];
            return p;
        }
    }


    inline double get_log(
        const LocalTree *tree, const States &states, int i, int j) const
    {
        return log(get(tree, states, i, j));
    }



    int ntimes;
    int nstates;
    bool own_data;
    double *B;
    double *D;
    double *E;
    double *G;
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
    
    inline double get(int i, int j) const
    {
        return exp(get_log(i, j));
    }


    inline double get_log(int i, int j) const
    {
        if (i == recoalsrc) {
            return recoalrow[j];
        } else if (i == recombsrc) {
            return recombrow[j];
        } else {
            if (determ[i] == j)
                return determprob[i];
            else
                return -INFINITY;
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


void calc_transition_probs(const LocalTree *tree, const ArgModel *model,
    const States &states, const LineageCounts *lineages, TransMatrix *matrix);
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
    TransMatrixSwitch *transmat_switch);
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


void calc_state_priors(const States &states, const LineageCounts *lineages, 
                       const ArgModel *model, double *priors,
                       const int minage=0);

void get_deterministic_transitions(
    const LocalTree *tree, const LocalTree *last_tree, 
    const Spr &spr, const int *mapping,
    const States &states1, const States &states2,
    int ntimes, int *next_states, bool internal=false);


bool assert_transmat(const LocalTree *tree, const ArgModel *model,
                     const TransMatrix *matrix);



} // namespace arghmm

#endif // ARGHMM_TRANS_H

//=============================================================================
// transitions

#ifndef ARGWEAVER_TRANS_H
#define ARGWEAVER_TRANS_H

#include "common.h"
#include "local_tree.h"
#include "model.h"
#include "states.h"

namespace argweaver {


// A compressed representation of the transition matrix.
//
// This transition matrix is used in the chromosome threading HMM within
// one non-recombining alignment block.
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
            delete [] lnB;
            delete [] lnE2;
            delete [] lnNegG1;
            delete [] G2;
            delete [] G3;
            delete [] lnG4;
            delete [] norecombs;
        }
    }

    // allocate space for transition matrix
    void allocate(int _ntimes)
    {
        ntimes = _ntimes;
        own_data = true;
        D = new double [ntimes];
        E = new double [ntimes];
        lnB = new double [ntimes];
        lnE2 = new double [ntimes];
        lnNegG1 = new double [ntimes];
        G2 = new double [ntimes];
        G3 = new double [ntimes];
        lnG4 = new double [ntimes];
        norecombs = new double [ntimes];
    }

    // Probability of transition from state i to state j.
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
        double term1 = D[a] * E[b];
        double minage_term = 0.0;
        if (minage > 0) {
          minage_term = exp(lnG4[b] + lnB[minage-1]);
        }

        if (!same_node) {
          if (a < b) {
            return term1 * (exp(lnE2[b] + lnB[a]) -
                            exp(lnE2[b] + lnNegG1[a]) 
                            - minage_term);
          } else if (a == b) {
            return term1 * ( (b > 0 ? exp(lnE2[b] + lnB[b-1]) : 0.0) +
                             G3[b] - minage_term);
          } else { // b < a
            return term1 * ( (b > 0 ? exp(lnE2[b] + lnB[b-1]) : 0.0)
                             + G2[b] - minage_term);
          }
        } else {
            double c_term = 0.0;
            if (c > 0) c_term = exp(lnG4[b] + lnB[c-1]);

            if (a < b) {
              return term1 * (2 * (exp(lnE2[b] + lnB[a]) -
                                   exp(lnE2[b] + lnNegG1[a]))
                              - c_term - minage_term);
            } else if (a == b) {
              return term1 * ((2 * ( (b > 0 ? exp(lnE2[b] + lnB[b-1]) : 0.0) + G3[b]))
                              - c_term - minage_term)
                + norecombs[a];
            } else { // b < a
              return term1 * ((2 * ( (b > 0 ? exp(lnE2[b] + lnB[b-1]) : 0.0) + G2[b]))
                              - c_term - minage_term);
            }
        }
    }

    // Log probability of transition from state i to state j.
    inline double get_log(
        const LocalTree *tree, const States &states, int i, int j) const
    {
        return log(get(tree, states, i, j));
    }


    int ntimes;     // Number of time steps in model.
    int nstates;    // Number of states in HMM.
    bool own_data;  // If true, delete matrix data when object is deleted

    double *D;      // Intermediate terms in calculating entries in the full
    double *E;      // transition matrix.
    double *lnB;
    double *lnE2;
    double *lnNegG1;
    double *G2;
    double *G3;
    double *lnG4;
    double *norecombs;

    bool internal;  // If true, this matrix is for threading an internal branch.
    int minage;     // Minimum age of a state we can consider (due to threading
                    // an internal branch).
};


// A compressed representation of the switch transition matrix.
//
// This transition matrix is used in the chromosome threading HMM to go between
// one non-recombining alignment block to the next (i.e. switching blocks).
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

    // Allocate matrix with dimensions (nstates1, nstates2).
    void allocate(int _nstates1, int _nstates2)
    {
        nstates1 = _nstates1;
        nstates2 = _nstates2;

        // NOTE: nstates1 and nstates2 might be zero
        // we still calculate transitions for a state space of size zero

        own_data = true;
        determ = new int [max(nstates1, 1)];
        determprob = new double [max(nstates1, 1)];
        recoalrow = new double [max(nstates2, 1)];
        recombrow = new double [max(nstates2, 1)];
    }

    // Log probability of transition from state i to state j.
    inline double get_log(int i, int j) const
    {
        return log(get(i, j));
    }

    // Probability of transition from state i to state j.
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


    int nstates1;   // Number of states in beginning block
    int nstates2;   // Number of states in the ending block
    int recoalsrc;  // Row for the state of the re-coalescence point
    int recombsrc;  // Row for the state of the recombination point
    bool own_data;  // If true, delete matrix data when object is deleted
    int *determ;    // Compressed representation of deterministic transitions
                    // determ[i] = j indicates (i --> j) is a deterministic
                    // transition.  determ[i] = -1, means row i has many
                    // transitions.
    double *determprob;  // Unnormalized probability determprob[i] of
                         // transition (i -> determ[i])
    double *recoalrow;   // Transition probabilities for row recoalsrc
    double *recombrow;   // Transition probabilities for row recombsrc
};



//=============================================================================

void calc_transition_probs(const LocalTree *tree, const ArgModel *model,
    const States &states, const LineageCounts *lineages, TransMatrix *matrix,
    bool internal=false, int minage=0);
void calc_transition_probs(const LocalTree *tree, const ArgModel *model,
                          const States &states, const LineageCounts *lineages,
                          double **transprob);
void get_transition_probs(const LocalTree *tree, const ArgModel *model,
                           const States &states, const LineageCounts *lineages,
                           const TransMatrix *matrix, double **transprob);


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



} // namespace argweaver

#endif // ARGWEAVER_TRANS_H

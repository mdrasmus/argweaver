//=============================================================================
// transitions

#ifndef ARGHMM_TRANS_H
#define ARGHMM_TRANS_H

#include "local_tree.h"
#include "model.h"
#include "states.h"

namespace arghmm {

// A compressed representation of the transition matrix
class TransMatrixCompress
{
public:
    TransMatrixCompress(int ntimes, int nstates, bool alloc=true) :
        ntimes(ntimes),
        nstates(nstates),
        own_data(false)
    {
        if (alloc)
            allocate(ntimes);
    }

    ~TransMatrixCompress()
    {
        if (own_data) {
            delete [] B;
            delete [] D;
            delete [] E;
            delete [] norecombs;
            delete [] sums;
        }
    }

    void allocate(int ntimes)
    {
        ntimes = ntimes;
        own_data = true;
        B = new double [ntimes];
        D = new double [ntimes];
        E = new double [ntimes];
        norecombs = new double [ntimes];
        sums = new double [nstates];
    }


    inline double get_transition_prob(LocalTree *tree, const States &states, 
                                      int i, int j)
    {
        const int node1 = states[i].node;
        const int a = states[i].time;
        const int c = tree->nodes[node1].age;
        const int node2 = states[j].node;
        const int b = states[j].time;
            
        if (node1 != node2)
            return D[a] * E[b] * B[min(a,b)] / sums[i];
        else {
            double p = D[a] * E[b] * (2 * B[min(a,b)] - B[min(c,b)]);
            if (a == b)
                p += norecombs[a];
            return p / sums[i];
        }
    }


    int ntimes;
    int nstates;
    bool own_data;
    double *B;
    double *D;
    double *E;
    double *norecombs;
    double *sums;
};


void calc_transition_probs(LocalTree *tree, ArgModel *model,
                           const States &states, LineageCounts *lineages,
                           double **transprob);

void calc_transition_probs_compress(LocalTree *tree, ArgModel *model,
    const States &states, LineageCounts *lineages, TransMatrixCompress *matrix);
void calc_transition_probs(LocalTree *tree, ArgModel *model,
                           const States &states, LineageCounts *lineages,
                           TransMatrixCompress *matrix, double **transprob);

void get_deterministic_transitions(
    LocalTree *tree, LocalTree *last_tree, const Spr &spr, int *mapping,
    const States &states1, const States &states2,
    int ntimes, int *next_states);

void calc_transition_probs_switch(
    LocalTree *tree, LocalTree *last_tree, const Spr &spr, int *mapping,
    const States &states1, const States &states2,
    ArgModel *model, LineageCounts *lineages, double **transprob);

void calc_state_priors(const States &states, LineageCounts *lineages, 
                       ArgModel *model, double *priors);


} // namespace arghmm

#endif // ARGHMM_TRANS_H

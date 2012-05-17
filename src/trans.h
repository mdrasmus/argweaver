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
    TransMatrixCompress(int ntimes, bool alloc=true) :
        ntimes(ntimes),
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
    }

    int ntimes;
    bool own_data;
    double *B;
    double *D;
    double *E;
    double *norecombs;
};


void calc_transition_probs(LocalTree *tree, ArgModel *model,
                           const States &states, LineageCounts *lineages,
                           double **transprob);

void calc_transition_probs_compress(LocalTree *tree, ArgModel *model,
    LineageCounts *lineages, TransMatrixCompress *matrix);

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

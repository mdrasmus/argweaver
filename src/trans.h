//=============================================================================
// transitions

#ifndef ARGHMM_TRANS_H
#define ARGHMM_TRANS_H

#include "local_tree.h"
#include "model.h"
#include "states.h"

namespace arghmm {

void calc_transition_probs(LocalTree *tree, ArgModel *model,
                           const States &states, LineageCounts *lineages,
                           double **transprob);

void calc_transition_probs_compressed(
    LocalTree *tree, ArgModel *model, LineageCounts *lineages,
    double **transprob);

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

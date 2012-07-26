#ifndef ARGHMM_TOTAL_PROB_H
#define ARGHMM_TOTAL_PROB_H

#include "local_tree.h"
#include "model.h"

namespace arghmm {

double calc_spr_prob(const ArgModel *model, const LocalTree *tree, 
                     const Spr &spr, LineageCounts &lineages);

double calc_arg_likelihood(const ArgModel *model, const Sequences *sequences, 
                           LocalTrees *trees);
double calc_arg_prior(const ArgModel *model, LocalTrees *trees);
double calc_arg_joint_prob(const ArgModel *model, const Sequences *sequences, 
                           LocalTrees *trees);



} // namespace arghmm

#endif // ARGHMM_TOTAL_PROB_H

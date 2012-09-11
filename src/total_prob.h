#ifndef ARGHMM_TOTAL_PROB_H
#define ARGHMM_TOTAL_PROB_H

#include "local_tree.h"
#include "model.h"

namespace arghmm {

double calc_spr_prob(const ArgModel *model, const LocalTree *tree, 
                     const Spr &spr, LineageCounts &lineages);

double calc_arg_likelihood(const ArgModel *model, const Sequences *sequences, 
                           const LocalTrees *trees);

// NOTE: trees should be uncompressed and sequences compressed
double calc_arg_likelihood(const ArgModel *model, const Sequences *sequences, 
                           const LocalTrees *trees, 
                           const SitesMapping* sites_mapping);

double calc_arg_prior(const ArgModel *model, const LocalTrees *trees);
double calc_arg_joint_prob(const ArgModel *model, const Sequences *sequences, 
                           const LocalTrees *trees);



} // namespace arghmm

#endif // ARGHMM_TOTAL_PROB_H

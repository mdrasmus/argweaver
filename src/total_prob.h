#ifndef ARGHMM_TOTAL_PROB_H
#define ARGHMM_TOTAL_PROB_H

#include "local_tree.h"
#include "model.h"

namespace arghmm {

double calc_spr_prob(const ArgModel *model, const LocalTree *tree, 
                     const Spr &spr, LineageCounts &lineages);

} // namespace arghmm

#endif // ARGHMM_TOTAL_PROB_H

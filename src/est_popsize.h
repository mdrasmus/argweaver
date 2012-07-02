#ifndef ARGHMM_EST_POPSIZE_H
#define ARGHMM_EST_POPSIZE_H

#include "local_tree.h"
#include "model.h"

namespace arghmm {

void est_popsize_local_trees(const ArgModel *model, const LocalTrees *trees, 
                             double *popsizes);

} // namespace arghmm

#endif // ARGHMM_EST_POPSIZE_H

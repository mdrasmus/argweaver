#ifndef ARGHMM_EST_POPSIZE_H
#define ARGHMM_EST_POPSIZE_H

#include "local_tree.h"
#include "model.h"

namespace argweaver {

void est_popsize_local_trees(const ArgModel *model, const LocalTrees *trees,
                             double *popsizes);

} // namespace argweaver

#endif // ARGHMM_EST_POPSIZE_H

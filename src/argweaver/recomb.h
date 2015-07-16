//=============================================================================
// Local trees

#ifndef ARGWEAVER_RECOMB_H
#define ARGWEAVER_RECOMB_H

// c++ includes
#include <assert.h>
#include <list>
#include <vector>
#include <string.h>
#include <stdio.h>

// argweaver includes
#include "local_tree.h"
#include "matrices.h"
#include "states.h"

namespace argweaver {

using namespace std;


double recomb_prob_unnormalized(const ArgModel *model, const LocalTree *tree,
                                const LineageCounts &lineages,
                                const State &last_state,
                                const State &state,
                                const NodePoint &recomb, bool internal);

void get_possible_recomb(const LocalTree *tree,
                         const State last_state, const State state,
                         bool internal, vector<NodePoint> &candidates);

void sample_recombinations(
    const LocalTrees *trees, const ArgModel *model,
    ArgHmmMatrixIter *matrix_list,
    int *thread_path, vector<int> &recomb_pos, vector<NodePoint> &recombs,
    bool internal=false);


} // namespace argweaver

#endif // ARGWEAVER_RECOMB_H

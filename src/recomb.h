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


namespace argweaver {

using namespace std;


void sample_recombinations(
    const LocalTrees *trees, const ArgModel *model,
    ArgHmmMatrixIter *matrix_list,
    int *thread_path, vector<int> &recomb_pos, vector<NodePoint> &recombs,
    bool internal=false);


} // namespace argweaver

#endif // ARGWEAVER_RECOMB_H

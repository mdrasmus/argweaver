//=============================================================================
// Local trees

#ifndef ARGHMM_RECOMB_H
#define ARGHMM_RECOMB_H

// c++ includes
#include <assert.h>
#include <list>
#include <vector>
#include <string.h>
#include <stdio.h>


namespace arghmm {

using namespace std;


void sample_recombinations(
    LocalTrees *trees, ArgModel *model, ArgHmmMatrixIter *matrix_list,
    int *thread_path, vector<int> &recomb_pos, vector<NodePoint> &recombs);

void max_recombinations(
    LocalTrees *trees, ArgModel *model, ArgHmmMatrixIter *matrix_list,
    int *thread_path, vector<int> &recomb_pos, vector<NodePoint> &recombs);


} // namespace arghmm

#endif // ARGHMM_RECOMB_H

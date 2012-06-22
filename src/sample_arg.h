//=============================================================================
// sample full ARGs
//

#ifndef ARGHMM_SAMPLE_ARG_H
#define ARGHMM_SAMPLE_ARG_H

// c++ includes
#include <vector>

// arghmm includes
#include "local_tree.h"
#include "model.h"
#include "sequences.h"


namespace arghmm {

using namespace std;

void sample_arg_seq(ArgModel *model, Sequences *sequences, LocalTrees *trees);

void resample_arg(ArgModel *model, Sequences *sequences, LocalTrees *trees,
                  int nremove=1);

void remax_arg(ArgModel *model, Sequences *sequences, LocalTrees *trees,
               int nremove=1);

void cond_sample_arg_seq(ArgModel *model, Sequences *sequences, 
                         LocalTrees *trees, 
                         LocalTree *start_tree, LocalTree *end_tree,
                         const vector<int> &full_seqids);

void sample_arg_seq_region(ArgModel *model, Sequences *sequences, 
                           LocalTrees *trees, int region_start, int region_end);

} // namespace arghmm

#endif // ARGHMM_SAMPLE_ARG_H

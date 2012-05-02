// c++ includes
#include <list>
#include <vector>
#include <string.h>

// arghmm includes
#include "common.h"
#include "emit.h"
#include "local_tree.h"
#include "seq.h"


using namespace std;

namespace arghmm {

double calc_arg_likelihood(ArgModel *model, Sequences *sequences, 
                           LocalTrees *trees)
{
    double lnl = 0.0;
    int nseqs = sequences->nseqs;
    double minlen = model->times[1];

    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it) {
        const int start = it->block.start;
        const int end = it->block.end;
        LocalTree *tree = it->tree;
        LocalNode *nodes = tree->nodes;
        int blocklen = end - start;

        // initialize per branch substitution counts
        int counts[tree->nnodes];
        for (int j=0; j<nseqs; j++)
            counts[j] = 0;

        // get subsequences
        char *seqs[nseqs];
        for (int j=0; j<nseqs; j++)
            seqs[j] = &sequences->seqs[j][start];

        // process block
        for (int i=start; i<end; i++) {
            char ancestral[tree->nnodes];
            parsimony_ancestral_seq(tree, seqs, nseqs, blocklen, i, ancestral);
            
            // count substituions per branch
            for (int j=0; j<tree->nnodes; j++) {
                if (j != tree->root && 
                    ancestral[j] != ancestral[nodes[j].parent])
                    counts[j]++;
            }
        }

        // compute likelihood of seeing each of those substituions

        // prior probability of root sequence
        lnl += log(.25) * blocklen;
        for (int j=0; j<tree->nnodes; j++) {
            double t = max(tree->get_dist(j), minlen);
            lnl += (counts[j] - blocklen) * model->mu * t +
                counts[j] * log(1/3. - 1/3. * exp(-model->mu * t));
        }
    }

    return lnl;
}


//=============================================================================
// C interface

extern "C" {

double arghmm_likelihood(
    int **ptrees, int **ages, int **sprs, int *blocklens,
    int ntrees, int nnodes, 
    double *times, int ntimes,
    double mu,
    char **seqs, int nseqs, int seqlen)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, NULL, 0.0, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    LocalTrees trees(ptrees, ages, sprs, blocklens, ntrees, nnodes);

    return calc_arg_likelihood(&model, &sequences, &trees);
}


} // extern "C"


} // namespace arghmm

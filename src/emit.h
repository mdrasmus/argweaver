#ifndef ARGWEAVER_EMIT_H
#define ARGWEAVER_EMIT_H

#include "local_tree.h"
#include "model.h"
#include "states.h"
#include "sequences.h"

namespace argweaver {

void find_masked_sites(const char *const *seqs, int nseqs, int seqlen,
                       bool *masked, bool *invariant=NULL);

void parsimony_ancestral_seq(const LocalTree *tree, const char *const *seqs,
                             int nseqs, int pos, char *ancestral,
                             int *postorder=NULL);
void parsimony_ancestral_set(const LocalTree *tree, const char * const *seqs,
                             int pos, int *postorder, int npostorder,
                             char *ancestral);
int parsimony_cost_seq(const LocalTree *tree, const char * const *seqs,
                       int nseqs, int pos, int *postorder);
void calc_emissions_external(const States &states, const LocalTree *tree,
                             const char * const*seqs, int nseqs, int seqlen,
                             const ArgModel *model, double **emit,
                             PhaseProbs *phase_pr);
void calc_emissions_internal(const States &states, const LocalTree *tree,
                             const char *const *seqs, int nseqs, int seqlen,
                             const ArgModel *model, double **emit);

double likelihood_tree(const LocalTree *tree, const ArgModel *model,
                       const char *const *seqs, const int nseqs,
                       const int start, const int end);

int count_noncompat(const LocalTrees *trees, const char * const *seqs,
                    int nseqs, int seqlen);


//=============================================================================
// C interface
extern "C" {

double **new_emissions(intstate *istates, int nstates,
                       int *ptree, int nnodes, int *ages_index,
                       char **seqs, int nseqs, int seqlen,
                       double *times, int ntimes,
                       double mu);
void delete_emissions(double **emit, int seqlen);

} // extern "C"

} // namesapce arghmm

#endif // ARGWEAVER_EMIT_H

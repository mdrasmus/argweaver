#ifndef ARGHMM_EMIT_H
#define ARGHMM_EMIT_H

#include "local_tree.h"
#include "model.h"
#include "states.h"

namespace arghmm {


void parsimony_ancestral_seq(const LocalTree *tree, const char *const *seqs, 
                             int nseqs, int pos, char *ancestral);
void calc_emissions(const States &states, const LocalTree *tree,
                    const char * const*seqs, int nseqs, int seqlen, 
                    const ArgModel *model, double **emit);
void calc_emissions_internal(const States &states, const LocalTree *tree,
                             const char *const *seqs, int nseqs, int seqlen, 
                             const ArgModel *model, double **emit);


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

#endif // ARGHMM_EMIT_H

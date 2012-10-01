

#include "matrices.h"

namespace arghmm {


// calculate transition and emission matrices for current block
void ArgHmmMatrixIter::calc_matrices(ArgHmmMatrices *matrices)
{
    if (internal)
        return calc_matrices_internal(matrices);

    // get block information
    const int blocklen = tree_iter->blocklen;
    matrices->blocklen = blocklen;
    const LocalTree *tree = tree_iter->tree;
    get_coal_states(tree, model->ntimes, *states);
    const int nstates = states->size();
    model->get_local_model(pos, local_model);
    
    // calculate emissions
    if (seqs) {
        const int nleaves = trees->get_num_leaves();
        char *subseqs[seqs->get_num_seqs()];
        for (int i=0; i<nleaves; i++)
            subseqs[i] = &seqs->seqs[trees->seqids[i]][pos];
        subseqs[nleaves] = &seqs->seqs[new_chrom][pos];
        matrices->emit = new_matrix<double>(blocklen, nstates);
        calc_emissions(*states, tree, subseqs, nleaves + 1, blocklen, 
                       &local_model, matrices->emit);
    } else {
        matrices->emit = NULL;
    }
    
    
    // calculate switch transition matrix if we are starting a new block
    if (!last_states) {
        // no previous states (first block)
        // therefore we have no switch transition matrix
        matrices->transmat_switch = NULL;
        matrices->nstates1 = matrices->nstates2 = nstates;
        
    } else {
        matrices->nstates1 = last_states->size();
        matrices->nstates2 = nstates;
        lineages.count(last_tree, internal);
        
        // calculate transmat_switch
        matrices->transmat_switch = new TransMatrixSwitch(
            matrices->nstates1, matrices->nstates2);
        calc_transition_probs_switch(tree, last_tree, 
                                     tree_iter->spr, tree_iter->mapping,
                                     *last_states, *states, &local_model, 
                                     &lineages, matrices->transmat_switch);
    }
    
    // update lineages to current tree
    lineages.count(tree, internal);
        
    // calculate transmat and use it for rest of block
    matrices->transmat = new TransMatrix(local_model.ntimes, nstates);
    calc_transition_probs(tree, &local_model, *states, &lineages, 
                          matrices->transmat);
}



// calculate transition and emission matrices for current block and
// for the internal branch resampling case
void ArgHmmMatrixIter::calc_matrices_internal(ArgHmmMatrices *matrices)
{
    // get block information
    const int blocklen = tree_iter->blocklen;
    matrices->blocklen = blocklen;
    const LocalTree *tree = tree_iter->tree;
    get_coal_states_internal(tree, model->ntimes, *states);
    const int nstates = states->size();
    model->get_local_model(pos, local_model);
    
    // calculate emissions
    if (seqs) {
        const int nleaves = trees->get_num_leaves();
        char *subseqs[nleaves];
        for (int i=0; i<nleaves; i++)
            subseqs[i] = &seqs->seqs[trees->seqids[i]][pos];
        matrices->emit = new_matrix<double>(blocklen, max(nstates, 1));
        calc_emissions_internal(*states, tree, subseqs, nleaves, 
                                blocklen, &local_model, matrices->emit);
    } else {
        matrices->emit = NULL;
    }
    
    
    // calculate switch transition matrix if we are starting a new block
    if (!last_states) {
        // no previous states (first block)
        // therefore we have no switch transition matrix
        matrices->transmat_switch = NULL;
        matrices->nstates1 = matrices->nstates2 = nstates;
        
    } else {
        matrices->nstates1 = last_states->size();
        matrices->nstates2 = nstates;
        lineages.count(last_tree, internal);
        
        // calculate transmat_switch
        matrices->transmat_switch = new TransMatrixSwitch(
            matrices->nstates1, matrices->nstates2);
        calc_transition_probs_switch_internal(tree, last_tree, 
            tree_iter->spr, tree_iter->mapping,
            *last_states, *states, &local_model, &lineages, 
            matrices->transmat_switch);
    }
    
    // update lineages to current tree
    lineages.count(tree, internal);
    
    // calculate transmat and use it for rest of block
    matrices->transmat = new TransMatrix(local_model.ntimes, nstates);
    calc_transition_probs_internal(tree, &local_model, *states, &lineages, 
                                   matrices->transmat);
}


} // namespace arghmm



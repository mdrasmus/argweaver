

#include "matrices.h"

namespace arghmm {


// calculate transition and emission matrices for current block
void ArgHmmMatrixIter::calc_matrices(ArgHmmMatrices *matrices)
{
    if (internal)
        return calc_matrices_internal(matrices);

    int blocklen = tree_iter->blocklen;
    LocalTree *tree = tree_iter->tree;
    get_coal_states(tree, model->ntimes, *states);
    int nstates = states->size();
    matrices->blocklen = blocklen;
    int nleaves = trees->get_num_leaves();
        
    // calculate emissions
    if (seqs) {
        char *subseqs[seqs->get_num_seqs()];
        for (int i=0; i<nleaves; i++)
            subseqs[i] = &seqs->seqs[trees->seqids[i]][pos];
        subseqs[nleaves] = &seqs->seqs[new_chrom][pos];
        matrices->emit = new_matrix<double>(blocklen, nstates);
        calc_emissions(*states, tree, subseqs, nleaves + 1, blocklen, 
                       model, matrices->emit);
    } else {
        matrices->emit = NULL;
    }
    
    
    // use switch matrix for first column of forward table
    // if we have a previous state space (i.e. not first block)
    if (!last_states) {
        matrices->transmat_switch = NULL;
        matrices->transprobs_switch = NULL;
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
                                     *last_states, *states, model, &lineages, 
                                     matrices->transmat_switch);

        if (calc_full) {
            matrices->transprobs_switch = new_matrix<double>(
                matrices->nstates1, matrices->nstates2);
            get_transition_probs_switch(matrices->transmat_switch,
                                        matrices->transprobs_switch);
        } else {
            matrices->transprobs_switch = NULL;
        }
    }
    
    // update lineages to current tree
    lineages.count(tree, internal);
        
    // calculate transmat and use it for rest of block
    matrices->transmat = new TransMatrix(model->ntimes, nstates);
    calc_transition_probs(tree, model, *states, &lineages, 
                          matrices->transmat);
    if (calc_full) {
        matrices->transprobs = new_matrix<double>(nstates, nstates);
        get_transition_probs(tree, model, *states, &lineages,
                             matrices->transmat, matrices->transprobs);
    } else {
        matrices->transprobs = NULL;
    }
}



// calculate transition and emission matrices for current block
void ArgHmmMatrixIter::calc_matrices_internal(ArgHmmMatrices *matrices)
{
    int blocklen = tree_iter->blocklen;
    LocalTree *tree = tree_iter->tree;
    get_coal_states_internal(tree, model->ntimes, *states);
    int nstates = states->size();
    matrices->blocklen = blocklen;
    int nleaves = trees->get_num_leaves();
    
    // calculate emissions
    if (seqs) {
        char *subseqs[nleaves];
        for (int i=0; i<nleaves; i++)
            subseqs[i] = &seqs->seqs[trees->seqids[i]][pos];
        matrices->emit = new_matrix<double>(blocklen, max(nstates, 1));
        calc_emissions_internal(*states, tree, subseqs, nleaves, 
                                blocklen, model, matrices->emit);
    } else {
        matrices->emit = NULL;
    }
    
    
    // use switch matrix for first column of forward table
    // if we have a previous state space (i.e. not first block)
    if (!last_states) {
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
            *last_states, *states, model, &lineages, 
            matrices->transmat_switch);
    }
    
    // update lineages to current tree
    lineages.count(tree, internal);
    
    // calculate transmat and use it for rest of block
    matrices->transmat = new TransMatrix(model->ntimes, nstates);
    calc_transition_probs_internal(tree, model, *states, &lineages, 
                                   matrices->transmat);
    
    // NOTE: full matrix calculation is not implemented for internal branch
    // sampling
    assert(!calc_full);
    matrices->transprobs = NULL;
    matrices->transprobs_switch = NULL;
}


} // namespace arghmm



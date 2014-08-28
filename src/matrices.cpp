

#include "matrices.h"

namespace argweaver {

// calculate transition and emission matrices for current block
void calc_arghmm_matrices_internal(
    const ArgModel *model, const Sequences *seqs, const LocalTrees *trees,
    const LocalTreeSpr *last_tree_spr, const LocalTreeSpr *tree_spr,
    const int start, const int end, int minage,
    ArgHmmMatrices *matrices, PhaseProbs *phase_pr)
{
    const bool internal = true;

    assert(phase_pr == NULL);  // for now focus on external case

    // get block information
    const int blocklen = end - start;
    matrices->blocklen = blocklen;
    const LocalTree *tree = tree_spr->tree;

    LineageCounts lineages(model->ntimes);
    States last_states;
    States states;
    matrices->states_model.set(model->ntimes, internal, minage);
    matrices->states_model.get_coal_states(tree, states);
    const int nstates = states.size();

    // calculate emissions
    if (seqs) {
        const int nleaves = trees->get_num_leaves();
        char *subseqs[nleaves];
	//	int phase_nodes[2]={-1,-1};
        for (int i=0; i<nleaves; i++)
            subseqs[i] = &seqs->seqs[trees->seqids[i]][start];
        matrices->emit = new_matrix<double>(blocklen, max(nstates, 1));
	/*	if (model->unphased) {
	    phase_nodes[0] = tree->nodes[tree->root].child[0];
	    while (!tree->nodes[phase_nodes[0]].is_leaf())
		phase_nodes[0] = tree->nodes[phase_nodes[0]].child[irand(2)];
	    for (unsigned int i=0; i < trees->seqids.size(); i++) {
	      if (trees->seqids[i] == seqs->get_pair(phase_nodes[0])) {
		phase_nodes[1] = i;
		break;
		}}}*/
	calc_emissions_internal(states, tree, subseqs, nleaves,
                                blocklen, model, matrices->emit);
    } else {
        matrices->emit = NULL;
    }


    // calculate switch transition matrix if we are starting a new block
    if (!last_tree_spr || last_tree_spr == tree_spr) {
        // no switch transition matrix
        matrices->transmat_switch = NULL;
        matrices->nstates1 = matrices->nstates2 = nstates;

    } else {
        LocalTree *last_tree = last_tree_spr->tree;
        matrices->states_model.get_coal_states(last_tree, last_states);
        matrices->nstates1 = last_states.size();
        matrices->nstates2 = nstates;
        lineages.count(last_tree, internal);

        // calculate transmat_switch
        matrices->transmat_switch = new TransMatrixSwitch(
            matrices->nstates1, matrices->nstates2);
        calc_transition_probs_switch_internal(tree, last_tree,
            tree_spr->spr, tree_spr->mapping,
            last_states, states, model, &lineages,
            matrices->transmat_switch);
    }

    // update lineages to current tree
    lineages.count(tree, internal);

    // calculate transmat and use it for rest of block
    matrices->transmat = new TransMatrix(model->ntimes, nstates);
    calc_transition_probs(tree, model, states, &lineages, matrices->transmat,
                          internal, matrices->states_model.minage);
}



// calculate transition and emission matrices for current block
void calc_arghmm_matrices_external(
    const ArgModel *model, const Sequences *seqs, const LocalTrees *trees,
    const LocalTreeSpr *last_tree_spr, const LocalTreeSpr *tree_spr,
    const int start, const int end, const int new_chrom,
    ArgHmmMatrices *matrices, PhaseProbs *phase_pr)
{
    // get block information
    const int blocklen = end - start;
    matrices->blocklen = blocklen;
    const LocalTree *tree = tree_spr->tree;

    LineageCounts lineages(model->ntimes);
    States last_states;
    States states;
    matrices->states_model.set(model->ntimes, false, 0);
    matrices->states_model.get_coal_states(tree, states);
    const int nstates = states.size();

    // calculate emissions
    if (seqs) {
        const int nleaves = trees->get_num_leaves();
        char *subseqs[seqs->get_num_seqs()];
        for (int i=0; i<nleaves; i++)
            subseqs[i] = &seqs->seqs[trees->seqids[i]][start];
        subseqs[nleaves] = &seqs->seqs[new_chrom][start];
        matrices->emit = new_matrix<double>(blocklen, nstates);
	if (model->unphased) {
	  //	    phase_pr->treemap1 = nleaves;
	    //	    phase_pr->updateTreeMap2(trees);
	    phase_pr->offset = start;
	    //	    printf("treemap = %i %i\n", phase_pr->treemap1, phase_pr->treemap2);
	}
        calc_emissions_external(states, tree, subseqs, nleaves + 1, blocklen,
                                model, matrices->emit, phase_pr);
    } else {
        matrices->emit = NULL;
    }


    // calculate switch transition matrix if we are starting a new block
    if (!last_tree_spr || last_tree_spr == tree_spr) {
        // no switch transition matrix
        matrices->transmat_switch = NULL;
        matrices->nstates1 = matrices->nstates2 = nstates;

    } else {
        LocalTree *last_tree = last_tree_spr->tree;
        matrices->states_model.get_coal_states(last_tree, last_states);
        matrices->nstates1 = last_states.size();
        matrices->nstates2 = nstates;
        lineages.count(last_tree);

        // calculate transmat_switch
        matrices->transmat_switch = new TransMatrixSwitch(
            matrices->nstates1, matrices->nstates2);
        calc_transition_probs_switch(tree, last_tree,
                                     tree_spr->spr, tree_spr->mapping,
                                     last_states, states, model,
                                     &lineages, matrices->transmat_switch);
    }

    // update lineages to current tree
    lineages.count(tree);

    // calculate transmat and use it for rest of block
    matrices->transmat = new TransMatrix(model->ntimes, nstates);
    calc_transition_probs(tree, model, states, &lineages,
                          matrices->transmat, false, matrices->states_model.minage);
}


void calc_arghmm_matrices(
    const ArgModel *model, const Sequences *seqs, const LocalTrees *trees,
    const LocalTreeSpr *last_tree_spr, const LocalTreeSpr *tree_spr,
    const int start, const int end, const int new_chrom,
    const StatesModel &states_model, ArgHmmMatrices *matrices,
    PhaseProbs *phase_pr)
{
    if (states_model.internal)
        calc_arghmm_matrices_internal(
            model, seqs, trees, last_tree_spr, tree_spr,
            start, end, states_model.minage, matrices,
            phase_pr);
    else
        calc_arghmm_matrices_external(
            model, seqs, trees, last_tree_spr,  tree_spr,
            start, end, new_chrom, matrices, phase_pr);
}


} // namespace argweaver


//=============================================================================
// sample full ARGs
//

// c++ includes
#include <vector>

// arghmm includes
#include "common.h"
#include "local_tree.h"
#include "logging.h"
#include "model.h"
#include "sample_arg.h"
#include "sample_thread.h"
#include "sequences.h"



namespace arghmm {

using namespace std;


// sequentially sample an ARG from scratch
// sequences are sampled in the order given
void sample_arg_seq(const ArgModel *model, const Sequences *sequences, 
                    LocalTrees *trees)
{
    const int nseqs = sequences->get_num_seqs();
    const int seqlen = sequences->length();

    if (trees->get_num_leaves() == 0) {
        // initialize ARG as trunk
        const int capacity = 2 * sequences->get_num_seqs() - 1;
        int start = trees->start_coord;
        int end = trees->end_coord;
        if (end != seqlen) {
            start = 0;
            end = seqlen;
        }
        trees->make_trunk(start, end, capacity);
    }
    
    // record which sequences are already in the tree
    bool has_sequence[nseqs];
    fill(has_sequence, has_sequence+nseqs, false);
    for (int i=0; i<trees->get_num_leaves(); i++)
        has_sequence[trees->seqids[i]] = true;

    // add more chromosomes one by one
    for (int new_chrom=0; new_chrom<nseqs; new_chrom++) {
        if (!has_sequence[new_chrom])
            sample_arg_thread(model, sequences, trees, new_chrom);
    }
}


/*
// sequentially sample an ARG from scratch
// sequences are sampled in the order given
void sample_arg_seq(const ArgModel *model, const Sequences *sequences, 
                    LocalTrees *trees)
{
    const int seqlen = sequences->length();

    // initialize ARG as trunk
    const int capacity = 2 * sequences->get_num_seqs() - 1;
    int start = trees->start_coord;
    int end = trees->end_coord;
    if (end != seqlen) {
        start = 0;
        end = seqlen;
    }
    trees->make_trunk(start, end, capacity);
    
    // add more chromosomes one by one
    for (int nchroms=2; nchroms<=sequences->get_num_seqs(); nchroms++) {
        // use first nchroms sequences
        //Sequences sequences2(sequences, nchroms);
        int new_chrom = nchroms - 1;
        sample_arg_thread(model, sequences, trees, new_chrom);
    }
}
 */

// resample the threading of all the chromosomes
void resample_arg(const ArgModel *model, const Sequences *sequences, 
                  LocalTrees *trees, int nremove)
{
    const int nleaves = trees->get_num_leaves();

    if (nremove == 1) {
        // cycle through chromosomes

        for (int chrom=0; chrom<nleaves; chrom++) {
            // remove chromosome from ARG and resample its thread
            remove_arg_thread(trees, chrom);
            sample_arg_thread(model, sequences, trees, chrom);
        }
    } else {
        // randomly choose which chromosomes to remove

        // clamp nremove
        if (nremove <= 0)
            return;
        if (nremove > nleaves - 1)
            nremove = nleaves - 1;

        // randomly choose which chromosomes to remove
        int chroms_avail[nleaves];
        for (int i=0; i<nleaves; i++)
            chroms_avail[i] = i;
        shuffle(chroms_avail, nleaves);

        // remove chromosomes from ARG
        for (int i=0; i<nremove; i++)
            remove_arg_thread(trees, chroms_avail[i]);

        // resample chromosomes
        for (int i=0; i<nremove; i++)
            sample_arg_thread(model, sequences, trees, chroms_avail[i]);
    }
}


// resample the threading of an internal branch
void resample_arg_all(const ArgModel *model, const Sequences *sequences, 
                      LocalTrees *trees, double prob_path_switch=.1)
{
    const int maxtime = model->ntimes + 1;
    int *removal_path = new int [trees->get_num_trees()];
    
    // ramdomly choose a removal path
    int node = irand(trees->nnodes);
    int pos = irand(trees->start_coord, trees->end_coord);
    sample_arg_removal_path(trees, node, pos, removal_path, prob_path_switch);
    
    remove_arg_thread_path(trees, removal_path, maxtime);
    sample_arg_thread_internal(model, sequences, trees);
    
    delete [] removal_path;
}


// resample the threading of a leaf of an ARG
void resample_arg_leaf(const ArgModel *model, const Sequences *sequences, 
                       LocalTrees *trees)
{
    const int maxtime = model->ntimes + 1;
    int *removal_path = new int [trees->get_num_trees()];
    
    // ramdomly choose a removal path
    int node = irand(trees->get_num_leaves());
    sample_arg_removal_leaf_path(trees, node, removal_path);
    
    remove_arg_thread_path(trees, removal_path, maxtime);
    sample_arg_thread_internal(model, sequences, trees);
    
    delete [] removal_path;
}



// resample the threading of an internal branch using MCMC
bool resample_arg_mcmc(const ArgModel *model, const Sequences *sequences, 
                       LocalTrees *trees)
{
    const int maxtime = model->ntimes + 1;
    int *removal_path = new int [trees->get_num_trees()];
    
    // save a copy of the local trees
    LocalTrees trees2;
    trees2.copy(*trees);

    // ramdomly choose a removal path
    double npaths = sample_arg_removal_path_uniform(trees, removal_path);
    remove_arg_thread_path(trees, removal_path, maxtime);
    sample_arg_thread_internal(model, sequences, trees);
    double npaths2 = count_total_arg_removal_paths(trees);

    // perform reject if needed
    double accept_prob = exp(npaths - npaths2);
    bool accept = (frand() < accept_prob);
    if (!accept)
        trees->copy(trees2);

    // logging
    printLog(LOG_LOW, "accept_prob = exp(%lf - %lf) = %f, accept = %d\n", 
             npaths, npaths2, accept_prob, (int) accept);
    

    // clean up
    delete [] removal_path;

    return accept;
}

// resample the threading of an internal branch using MCMC
void resample_arg_mcmc_all(const ArgModel *model, const Sequences *sequences, 
                           LocalTrees *trees, double frac_leaf,
                           int window, int step, int niters)
{
    if (frand() < frac_leaf) {
        resample_arg_leaf(model, sequences, trees);
        printLog(LOG_LOW, "resample_arg_leaf: accept=%f\n", 1.0);
    } else {
        double accept_rate = resample_arg_regions(
            model, sequences, trees, window, step, niters);
        printLog(LOG_LOW, "resample_arg_regions: accept=%f\n", accept_rate);
    }
}



// resample the threading of an internal branch with preference for recombs
void resample_arg_recomb(const ArgModel *model, const Sequences *sequences, 
                         LocalTrees *trees, double recomb_preference)
{
    const int maxtime = model->ntimes + 1;
    int *removal_path = new int [trees->get_num_trees()];
    
    // ramdomly choose a removal path weighted by recombinations
    sample_arg_removal_path_recomb(trees, recomb_preference, removal_path);
    remove_arg_thread_path(trees, removal_path, maxtime);
    sample_arg_thread_internal(model, sequences, trees);
    
    delete [] removal_path;
}



// resample an ARG heuristically and aggressively to high joint probability
void resample_arg_climb(const ArgModel *model, const Sequences *sequences, 
                        LocalTrees *trees, double recomb_preference)
{
    resample_arg_recomb(model, sequences, trees, recomb_preference);
}




/*
// sample an ARG with both sequential and gibbs iterations
void sample_arg_seq_gibbs(const ArgModel *model, const Sequences *sequences, 
                          LocalTrees *trees, int seqiters, int gibbsiters)
{
    const int nseqs = sequences->get_num_seqs();
    const int seqlen = sequences->length();
    const int minseqs = 3;

    // initialize ARG as trunk
    const int capacity = 2 * nseqs - 1;
    trees->make_trunk(0, seqlen, capacity);
    
    int nleaves = 1;
    Timer time;
    while (nleaves < nseqs) {
        time.start();

        // sequential stage
        int nleaves2 = max(min(nleaves + seqiters, nseqs), minseqs);

        // add more chromosomes one by one
        for (int nchroms=nleaves+1; nchroms<=nleaves2; nchroms++) {
            // use first nchroms sequences
            Sequences sequences2(sequences, nchroms);
            int new_chrom = nchroms - 1;
            sample_arg_thread(model, &sequences2, trees, new_chrom);
        }
        nleaves = nleaves2;


        // gibbs stage 
        // randomly choose gibbsiters chromosomes
        int chroms[nleaves];
        for (int i=0; i<nleaves; i++)
            chroms[i] = i;
        shuffle(chroms, nleaves);

        for (int i=0; i<min(gibbsiters, nleaves); i++)
            resample_arg_all(model, sequences, trees);

        printTimerLog(time, LOG_QUIET, 
                      "seq_gibbs stage (%3d leaves):       ", nleaves);
    }
}
*/


//=============================================================================
// sub-region resampling


State find_state_sub_tree(
    const LocalTree *full_tree, const vector<int> &full_seqids,
    const LocalTree *partial_tree, const vector<int> &partial_seqids, 
    int new_chrom)
{
    // reconcile full tree to partial tree
    int recon[full_tree->nnodes];
    map_congruent_trees(full_tree, &full_seqids[0],
                        partial_tree, &partial_seqids[0], recon);
    
    // find new chrom in full_tree
    int ptr = find_array(&full_seqids[0], full_seqids.size(), new_chrom);
    assert(ptr != -1);

    // walk up from new chrom until we hit a reconciled portion of full tree
    while (recon[ptr] == -1)
        ptr = full_tree->nodes[ptr].parent;

    return State(recon[ptr], full_tree->nodes[ptr].age);
}



// sequentially sample an ARG from scratch
// sequences are sampled in the order given
void cond_sample_arg_seq(const ArgModel *model, const Sequences *sequences, 
                         LocalTrees *trees, 
                         LocalTree *start_tree, LocalTree *end_tree,
                         const vector<int> &full_seqids)
{
    // initialize ARG as trunk
    const int capacity = 2 * sequences->get_num_seqs() - 1;
    trees->make_trunk(trees->start_coord, trees->end_coord, capacity);
    
    // add more chromosomes one by one
    for (int nchroms=2; nchroms<=sequences->get_num_seqs(); nchroms++) {
        // use first nchroms sequences
        Sequences sequences2(sequences, nchroms);
        int new_chrom = nchroms - 1;

        // determine start and end states from given trees
        LocalTree *first_tree = trees->front().tree;
        LocalTree *last_tree = trees->back().tree;
        State start_state = find_state_sub_tree(
            start_tree, full_seqids, first_tree, trees->seqids, new_chrom);
        State end_state = find_state_sub_tree(
            end_tree, full_seqids, last_tree, trees->seqids, new_chrom);

        cond_sample_arg_thread(model, &sequences2, trees, new_chrom,
                               start_state, end_state);
        
        assert_trees(trees);
    }
}



// sequentially sample an ARG only for a given region
// sequences are sampled in the order given
void sample_arg_seq_region(const ArgModel *model, const Sequences *sequences, 
                           LocalTrees *trees, int region_start, int region_end)
{
    // ensure region is within ARG
    assert(region_start > trees->start_coord);
    assert(region_end < trees->end_coord);
    assert(region_start < region_end);

    // partion trees into three segments
    LocalTrees *trees2 = partition_local_trees(trees, region_start);
    LocalTrees *trees3 = partition_local_trees(trees2, region_end);
    assert(trees2->length() == region_end - region_start);
    
    // resample region conditioning on starting and ending trees
    cond_sample_arg_seq(model, sequences, 
                        trees2, trees->back().tree, trees3->front().tree,
                        trees->seqids);

    // rejoin trees
    append_local_trees(trees, trees2);
    append_local_trees(trees, trees3);
    
    // clean up
    delete trees2;
    delete trees3;
}


State find_state_sub_tree_internal(
    const LocalTree *full_tree, const LocalTree *partial_tree, int maxtime)
{
    if (partial_tree->nodes[partial_tree->root].age < maxtime) {
        // fully specified tree
        return State(-1, -1);
    }

    // NOTE: do not assume internal nodes have same naming scheme between
    // trees

    int subtree_root = partial_tree->nodes[partial_tree->root].child[0];

    // identify node by path length from left most leaf
    int count = 0;
    int leaf = subtree_root;
    while (!partial_tree->nodes[leaf].is_leaf()) {
        leaf = partial_tree->nodes[leaf].child[0];
        count++;
    }

    // find equivalent node in full tree
    int ptr = leaf;
    while (count > 0) {
        ptr = full_tree->nodes[ptr].parent;
        count--;
    }

    // find sibling and age of coalescence
    int sib = full_tree->get_sibling(ptr);
    assert(sib != -1);
    int parent = full_tree->nodes[ptr].parent;
    int coal_time = full_tree->nodes[parent].age;

    // identify sibling by leaf and path length
    count = 0;
    leaf = sib;
    while (!full_tree->nodes[leaf].is_leaf()) {
        leaf = full_tree->nodes[leaf].child[0];
        count++;
    }

    // map sib back to partial tree
    ptr = leaf;
    while (count > 0) {
        ptr = partial_tree->nodes[ptr].parent;
        count--;
    }

    return State(ptr, coal_time);
}


// resample an ARG only for a given region
// all branches are possible to resample
// open_ended -- If true and region touches start or end of local trees do not
//               conditioned on state.
double resample_arg_region(
    const ArgModel *model, const Sequences *sequences, 
    LocalTrees *trees, int region_start, int region_end, int niters,
    bool open_ended)
{
    const int maxtime = model->ntimes + 1;

    // special case: zero length region
    if (region_start == region_end)
        return 1.0;
    
    // assert region is within trees
    assert(region_start >= trees->start_coord);
    assert(region_end <= trees->end_coord);
    assert(region_start < region_end);

    // partion trees into three segments
    LocalTrees *trees2 = partition_local_trees(trees, region_start);
    LocalTrees *trees3 = partition_local_trees(trees2, region_end);
    assert(trees2->length() == region_end - region_start);

    // TODO: refactor
    // extend stub (zero length block) if it happens to exist
    bool stub = (trees2->trees.back().blocklen == 0);
    if (stub) {
        trees2->trees.back().blocklen += 1;
        trees2->end_coord++;
    }

    // perform several iterations of resampling
    int accepts = 0;
    for (int i=0; i<niters; i++) {
        printLog(LOG_LOW, "region sample: iter=%d, region=(%d, %d)\n", 
                 i, region_start, region_end);

        // save a copy of the local trees
        LocalTrees old_trees2;
        old_trees2.copy(*trees2);

        // get starting and ending trees
        LocalTree start_tree(*trees2->front().tree);
        LocalTree end_tree(*trees2->back().tree);

        // remove internal branch from trees2
        int *removal_path = new int [trees2->get_num_trees()];
        double npaths = sample_arg_removal_path_uniform(trees2, removal_path);
        remove_arg_thread_path(trees2, removal_path, maxtime);
        delete [] removal_path;
        
        // determine start and end states from start and end trees
        LocalTree *start_tree_partial = trees2->front().tree;
        LocalTree *end_tree_partial = trees2->back().tree;
        State start_state = find_state_sub_tree_internal(
            &start_tree, start_tree_partial, maxtime);
        State end_state = find_state_sub_tree_internal(
            &end_tree, end_tree_partial, maxtime);
        
        // set start/end state to null if open ended is requested
        if (open_ended) {
            if (region_start == trees->start_coord)
                start_state.set_null();
            if (region_end == trees3->end_coord)
                end_state.set_null();
        }

        // sample new ARG conditional on start and end states
        decLogLevel();
        cond_sample_arg_thread_internal(model, sequences, trees2,
                                        start_state, end_state);
        incLogLevel();
        assert_trees(trees2);
        double npaths2 = count_total_arg_removal_paths(trees2);
        
        // perform reject if needed
        double accept_prob = exp(npaths - npaths2);
        bool accept = (frand() < accept_prob);
        if (!accept)
            trees2->copy(old_trees2);
        else
            accepts++;
        
        // logging
        printLog(LOG_LOW, "accept_prob = exp(%lf - %lf) = %f, accept = %d\n", 
                 npaths, npaths2, accept_prob, (int) accept);
    }

    // remove stub if it exists
    if (stub) {
        trees2->trees.back().blocklen -= 1;
        trees2->end_coord--;
    }
    
    // rejoin trees
    append_local_trees(trees, trees2);
    append_local_trees(trees, trees3);

    // clean up
    delete trees2;
    delete trees3;

    return accepts / double(niters);
}


// resample an ARG a region at a time in a sliding window
double resample_arg_regions(
    const ArgModel *model, const Sequences *sequences, 
    LocalTrees *trees, int window, int step, int niters)
{
    decLogLevel();
    double accept_rate = 0.0;
    int nwindows = 0;
    for (int start=trees->start_coord; 
         start == trees->start_coord || start+window/2 <trees->end_coord; 
         start+=step)
    {
        nwindows++;
        int end = min(start + window, trees->end_coord);
        accept_rate += resample_arg_region(
            model, sequences, trees, start, end, niters);
    }
    incLogLevel();

    accept_rate /= nwindows;
    return accept_rate;
}




// resample the threading of all the chromosomes
void remax_arg(const ArgModel *model, const Sequences *sequences, 
               LocalTrees *trees, int nremove)
{
    const int nleaves = trees->get_num_leaves();

    if (nremove == 1) {
        // cycle through chromosomes

        for (int chrom=0; chrom<nleaves; chrom++) {
            // remove chromosome from ARG and resample its thread
            remove_arg_thread(trees, chrom);
            max_arg_thread(model, sequences, trees, chrom);
        }
    } else {
        // clamp nremove
        if (nremove <= 0)
            return;
        if (nremove > nleaves)
            nremove = nleaves;

        // randomly choose which chromosomes to remove
        int chroms_avail[nleaves];
        for (int i=0; i<nleaves; i++)
            chroms_avail[i] = i;
        shuffle(chroms_avail, nleaves);

        // remove chromosomes from ARG
        for (int i=0; i<nremove; i++) {
            remove_arg_thread(trees, chroms_avail[i]);
        }

        // resample chromosomes
        for (int i=0; i<nremove; i++)
            max_arg_thread(model, sequences, trees, chroms_avail[i]);
    }
}



//=============================================================================
// C interface
extern "C" {

// sequentially sample until all chromosomes are present
LocalTrees *arghmm_complete_arg(
    LocalTrees *trees, ArgModel *model, Sequences *sequences)
{
    const int nseqs = sequences->get_num_seqs();
    for (int new_chrom=trees->get_num_leaves(); new_chrom<nseqs; new_chrom++)
        sample_arg_thread(model, sequences, trees, new_chrom);
    return trees;
}


// sequentially sample an ARG
LocalTrees *arghmm_sample_arg_seq(
    double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    LocalTrees *trees = new LocalTrees();    

    sample_arg_seq(&model, &sequences, trees);

    return trees;
}


// sequentially sample an ARG and then refine with gibbs
LocalTrees *arghmm_sample_arg_refine(
    double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, int niters, int nremove)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    LocalTrees *trees = new LocalTrees();    

    sample_arg_seq(&model, &sequences, trees);
    for (int i=0; i<niters; i++)
        resample_arg(&model, &sequences, trees, nremove);
    
    return trees;
}


// resample an ARG with gibbs
LocalTrees *arghmm_resample_arg(
    LocalTrees *trees, double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, int niters, int nremove)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    
    // sequentially sample until all chromosomes are present
    arghmm_complete_arg(trees, &model, &sequences);

    // gibbs sample
    for (int i=0; i<niters; i++)
        resample_arg(&model, &sequences, trees, nremove);
    
    return trees;
}


// resample all branches in an ARG with gibbs
LocalTrees *arghmm_resample_all_arg(
    LocalTrees *trees, double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, int niters, double prob_path_switch)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    
    // sequentially sample until all chromosomes are present
    arghmm_complete_arg(trees, &model, &sequences);
    
    // gibbs sample
    for (int i=0; i<niters; i++)
        resample_arg_all(&model, &sequences, trees, prob_path_switch);
    
    return trees;
}


// resample all branches in an ARG with gibbs
LocalTrees *arghmm_resample_mcmc_arg(
    LocalTrees *trees, double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, int niters, int niters2, int window)
{
    // setup model, local trees, sequences
    double frac_leaf = 0.5;
    int step = window / 2;
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    
    // sequentially sample until all chromosomes are present
    arghmm_complete_arg(trees, &model, &sequences);
    
    // gibbs sample
    for (int i=0; i<niters; i++) {
        printLog(LOG_LOW, "sample %d\n", i);
        resample_arg_mcmc_all(&model, &sequences, trees, frac_leaf,
                              window, step, niters2);
    }

    return trees;
}


LocalTrees *arghmm_resample_arg_leaf(
    LocalTrees *trees, double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, int niters)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    
    // sequentially sample until all chromosomes are present
    arghmm_complete_arg(trees, &model, &sequences);
    
    // gibbs sample
    for (int i=0; i<niters; i++)
        resample_arg_leaf(&model, &sequences, trees);
    
    return trees;    
}


// resample ARG focused on recombinations
LocalTrees *arghmm_resample_climb_arg(
    LocalTrees *trees, double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, int niters, double recomb_preference)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    
    // sequentially sample until all chromosomes are present
    arghmm_complete_arg(trees, &model, &sequences);
    
    // gibbs sample
    for (int i=0; i<niters; i++)
        resample_arg_climb(&model, &sequences, trees, recomb_preference);
    
    return trees;
}


LocalTrees *arghmm_resample_arg_region(
    LocalTrees *trees, double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, 
    int region_start, int region_end, int niters)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);

    resample_arg_region(&model, &sequences, trees, 
                        region_start, region_end, niters);
    
    return trees;
}


// remax an ARG with viterbi
LocalTrees *arghmm_remax_arg(
    LocalTrees *trees, double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, int niters, int nremove)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    
    // sequentially sample until all chromosomes are present
    arghmm_complete_arg(trees, &model, &sequences);
    
    // gibbs sample
    for (int i=0; i<niters; i++)
        remax_arg(&model, &sequences, trees, nremove);
    
    return trees;
}


// DISABLED
LocalTrees *arghmm_sample_arg_seq_gibbs(double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, int seqiters, int gibbsiters)
{
    // setup model, local trees, and sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    LocalTrees *trees = new LocalTrees();    

    //sample_arg_seq_gibbs(&model, &sequences, trees, seqiters, gibbsiters);

    return trees;
}


} // extern C

} // namespace arghmm

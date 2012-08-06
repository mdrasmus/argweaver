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
        Sequences sequences2(sequences, nchroms);
        int new_chrom = nchroms - 1;
        sample_arg_thread(model, &sequences2, trees, new_chrom);
    }
}


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
        for (int i=0; i<nremove; i++) {
            remove_arg_thread(trees, chroms_avail[i]);
        }

        // resample chromosomes
        for (int i=0; i<nremove; i++)
            sample_arg_thread(model, sequences, trees, chroms_avail[i]);
    }
}


// resample the threading of an internal branch
void resample_arg_all(const ArgModel *model, const Sequences *sequences, 
                      LocalTrees *trees)
{
    const int maxtime = model->ntimes + 1;
    int *removal_path = new int [trees->get_num_trees()];
    
    // ramdomly choose a removal path
    int node = irand(trees->nnodes);
    int pos = irand(trees->start_coord, trees->end_coord);
    sample_arg_removal_path(trees, node, pos, removal_path);
    
    remove_arg_thread_path(trees, removal_path, maxtime);
    sample_arg_thread_internal(model, sequences, trees);
    
    delete [] removal_path;
}



// resample the threading of an internal branch
void resample_arg_climb(const ArgModel *model, const Sequences *sequences, 
                        LocalTrees *trees, double recomb_preference)
{
    const int maxtime = model->ntimes + 1;
    int *removal_path = new int [trees->get_num_trees()];
    
    if (true) { //frand() < .5) {
        // ramdomly choose a removal path weighted by recombinations
        sample_arg_removal_path_recomb(trees, recomb_preference, removal_path);

    } else {
        int node = irand(trees->get_num_leaves());
        sample_arg_removal_leaf_path(trees, node, removal_path);        
    }

    remove_arg_thread_path(trees, removal_path, maxtime);
    sample_arg_thread_internal(model, sequences, trees);
    
    delete [] removal_path;
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

        for (int i=0; i<min(gibbsiters, nleaves); i++) {
            resample_arg_all(model, sequences, trees);

            /*
            remove_arg_thread(trees, chroms[i]);
            sample_arg_thread(model, sequences, trees, chroms[i]);
            */
        }

        printTimerLog(time, LOG_QUIET, 
                      "seq_gibbs stage (%3d leaves):       ", nleaves);
    }
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

    int subtree_root = partial_tree->nodes[partial_tree->root].child[0];
    int sib = full_tree->get_sibling(subtree_root);
    int parent = full_tree->nodes[subtree_root].parent;

    return State(sib, full_tree->nodes[parent].age);
}


// resample an ARG only for a given region
// all branches are possible to resample
void resample_arg_all_region(
    const ArgModel *model, const Sequences *sequences, 
    LocalTrees *trees, int region_start, int region_end, int niters)
{
    const int maxtime = model->ntimes + 1;
    const double recomb_preference = .7;

    assert(region_start > trees->start_coord);
    assert(region_end < trees->end_coord);
    assert(region_start < region_end);

    // partion trees into three segments
    LocalTrees *trees2 = partition_local_trees(trees, region_start);
    LocalTrees *trees3 = partition_local_trees(trees2, region_end);

    assert(trees2->length() == region_end - region_start);
    
    LocalTree *start_tree = trees->back().tree;
    LocalTree *end_tree = trees3->front().tree;

    // perform several iterations of resampling
    int *removal_path = new int [trees2->get_num_trees()];
    for (int i=0; i<niters; i++) {
        // remove internal branch from trees2
        // ramdomly choose a removal path weighted by recombinations

        sample_arg_removal_path_recomb(trees2, recomb_preference, removal_path);
        remove_arg_thread_path(trees2, removal_path, maxtime);

        // determine start and end states from given trees
        LocalTree *first_tree = trees2->front().tree;
        LocalTree *last_tree = trees2->back().tree;
        State start_state = find_state_sub_tree_internal(
            start_tree, first_tree, maxtime);
        State end_state = find_state_sub_tree_internal(
            end_tree, last_tree, maxtime);

        cond_sample_arg_thread_internal(model, sequences, trees2,
                                        start_state, end_state);

        assert_trees(trees2);
    }
    
    // rejoin trees
    append_local_trees(trees, trees2);
    append_local_trees(trees, trees3);
    
    // clean up
    delete [] removal_path;
    delete trees2;
    delete trees3;
}


//=============================================================================
// C interface
extern "C" {


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
    for (int new_chrom=trees->get_num_leaves(); new_chrom<nseqs; new_chrom++)
        sample_arg_thread(&model, &sequences, trees, new_chrom);

    // gibbs sample
    for (int i=0; i<niters; i++)
        resample_arg(&model, &sequences, trees, nremove);
    
    return trees;
}


// resample all branches in an ARG with gibbs
LocalTrees *arghmm_resample_all_arg(
    LocalTrees *trees, double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, int niters)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    
    // sequentially sample until all chromosomes are present
    for (int new_chrom=trees->get_num_leaves(); new_chrom<nseqs; new_chrom++)
        sample_arg_thread(&model, &sequences, trees, new_chrom);

    // gibbs sample
    for (int i=0; i<niters; i++)
        resample_arg_all(&model, &sequences, trees);
    
    return trees;
}



// resample all branches in an ARG with gibbs
LocalTrees *arghmm_resample_climb_arg(
    LocalTrees *trees, double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, int niters, double recomb_preference)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    
    // sequentially sample until all chromosomes are present
    for (int new_chrom=trees->get_num_leaves(); new_chrom<nseqs; new_chrom++)
        sample_arg_thread(&model, &sequences, trees, new_chrom);

    // gibbs sample
    for (int i=0; i<niters; i++)
        resample_arg_climb(&model, &sequences, trees, recomb_preference);
    
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
    for (int new_chrom=trees->get_num_leaves(); new_chrom<nseqs; new_chrom++) {
        max_arg_thread(&model, &sequences, trees, new_chrom);
    }

    // gibbs sample
    for (int i=0; i<niters; i++)
        remax_arg(&model, &sequences, trees, nremove);
    
    return trees;
}


LocalTrees *arghmm_resample_arg_region(
    LocalTrees *trees, double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, int region_start, int region_end)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    
    sample_arg_seq_region(&model, &sequences, trees, region_start, region_end);
    
    return trees;
}

LocalTrees *arghmm_resample_arg_all_region(
    LocalTrees *trees, double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, 
    int region_start, int region_end, int niters)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);

    resample_arg_all_region(&model, &sequences, trees, 
                            region_start, region_end, niters);
    
    return trees;
}



LocalTrees *arghmm_sample_arg_seq_gibbs(double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, int seqiters, int gibbsiters)
{
    // setup model, local trees, and sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    LocalTrees *trees = new LocalTrees();    

    sample_arg_seq_gibbs(&model, &sequences, trees, seqiters, gibbsiters);

    return trees;
}


} // extern C

} // namespace arghmm

// c++ includes
#include <list>
#include <vector>
#include <string.h>

// arghmm includes
#include "common.h"
#include "emit.h"
#include "hmm.h"
#include "local_tree.h"
#include "logging.h"
#include "matrices.h"
#include "model.h"
#include "recomb.h"
#include "sample_thread.h"
#include "sequences.h"
#include "sequences.h"
#include "states.h"
#include "thread.h"
#include "trans.h"



namespace arghmm {

using namespace std;


//=============================================================================
// Forward algorithm for thread path

inline double logsum_fast(const double *vals, int nvals)
{
    return logsum(vals, nvals, -10);
}


// compute one block of forward algorithm with compressed transition matrices
void arghmm_forward_alg_block(const LocalTree *tree, const ArgModel *model,
                              const int blocklen, const States &states, 
                              const LineageCounts &lineages,
                              const TransMatrix *matrix,
                              const double* const *emit, double **fw,
                              bool internal=false)
{
    const int nstates = states.size();
    const int ntimes = model->ntimes;
    const LocalNode *nodes = tree->nodes;

    // get aliases for various precomputed terms
    const double *B = matrix->B;
    const double *D = matrix->D;
    const double *E = matrix->E;
    const double *G = matrix->G;
    const double *norecombs = matrix->norecombs;

    double Bq = 0.0;
    int minage = 0;
    int subtree_root = 0, maintree_root = 0;
    if (internal) {
        subtree_root = nodes[tree->root].child[0];
        maintree_root = nodes[tree->root].child[1];
        const int subtree_age = nodes[subtree_root].age;
        if (subtree_age > 0) {
            Bq = B[subtree_age - 1];
            minage = subtree_age;
        }

        if (nstates == 0) {
            // handle fully given case
            for (int i=1; i<blocklen; i++)
                fw[i][0] = fw[i-1][0];
            return;
        }
    }

    // compute ntimes*ntimes and ntime*nstates temp matrices
    double lnorecombs[ntimes];
    double tmatrix[ntimes][ntimes];
    double tmatrix2[ntimes][nstates];
    for (int a=0; a<ntimes-1; a++) {
        lnorecombs[a] = log(norecombs[a]);

        for (int b=0; b<ntimes-1; b++) {
            const double I = double(a <= b);
            if (a < minage || b < minage)
                tmatrix[a][b] = -INFINITY;
            else
                tmatrix[a][b] = log(D[a] * E[b] * (B[min(a,b)] -Bq - I*G[a]));
            assert(!isnan(tmatrix[a][b]));
        }

        for (int k=0; k<nstates; k++) {
            const int b = states[k].time;
            const int node2 = states[k].node;
            const int c = nodes[node2].age;
            const double Bc = (c > 0 ? B[c-1] : 0.0);
            const double I = double(a <= b);
            assert(b >= minage);
            if (a < minage)
                tmatrix2[a][k] = -INFINITY;
            else
                tmatrix2[a][k] = log(D[a] * E[b] * (B[min(a,b)] - I*G[a] - Bc));
        }
    }

    // get max time
    int maxtime = 0;
    for (int k=0; k<nstates; k++)
        if (maxtime < states[k].time)
            maxtime = states[k].time;

    // get branch ages
    int ages1[tree->nnodes];
    int ages2[tree->nnodes];
    for (int i=0; i<tree->nnodes; i++) {
        ages1[i] = max(nodes[i].age, minage);
        if (internal)
            ages2[i] = (i == maintree_root) ? 
                maxtime : nodes[nodes[i].parent].age;
        else
            ages2[i] = (i == tree->root) ? maxtime : nodes[nodes[i].parent].age;
    }
    

    // group states by age
    int nstates_per_time[ntimes];
    fill(nstates_per_time, nstates_per_time+ntimes, 0);
    for (int j=0; j<nstates; j++)
        nstates_per_time[states[j].time]++;
    int offset[ntimes], offset2[ntimes];
    offset[0] = offset2[0] = 0;
    for (int a=1; a<ntimes-1; a++)
        offset[a] = offset2[a] = offset[a-1] + nstates_per_time[a-1];
    int state_map[nstates];
    for (int j=0; j<nstates; j++) {
        const int a = states[j].time;
        state_map[j] = offset2[a];
        offset2[a]++;
    }


    NodeStateLookup state_lookup(states, tree->nnodes);

    double vec[nstates + ntimes];
    double tmatrix_fgroups[ntimes];
    double fgroups[ntimes];
    for (int i=1; i<blocklen; i++) {
        const double *col1 = fw[i-1];
        double *col2 = fw[i];
        const double *emit2 = emit[i];
        
        // precompute the fgroup sums
        for (int j=0; j<nstates; j++)
            vec[state_map[j]] = col1[j];
        for (int a=0; a<ntimes-1; a++)
            fgroups[a] = logsum_fast(&vec[offset[a]], nstates_per_time[a]);
        
        // multiply tmatrix and fgroups together
        for (int b=0; b<ntimes-1; b++) {
            for (int a=0; a<ntimes-1; a++)
                vec[a] = tmatrix[a][b] + fgroups[a];
            tmatrix_fgroups[b] = logsum_fast(vec, ntimes-1);
        }
        
        // fill in one column of forward table
        for (int k=0; k<nstates; k++) {
            const int b = states[k].time;
            const int node2 = states[k].node;
            const int age1 = ages1[node2];
            const int age2 = ages2[node2];
            
            // same branch case (extra terms substracted, no 2*B[min(a,b)])
            vec[0] = tmatrix_fgroups[b];
            int m = 1;
            const int j1 = state_lookup.lookup(node2, age1);
            const int j2 = state_lookup.lookup(node2, age2);
            assert(j2 != -1 && j1 != -1);
            for (int j=j1, a=age1; j<=j2; j++, a++)
                vec[m++] = tmatrix2[a][k] + col1[j];
            
            // same state case (add possibility of no recomb)
            vec[m++] = lnorecombs[b] + col1[k];
            
            col2[k] = logsum_fast(vec, m) + emit2[k];
        }
    }
}


// run forward algorithm for one column of the table
// use switch matrix
void arghmm_forward_switch(const double *col1, double* col2, 
                           const TransMatrixSwitch *matrix,
                           const double *emit)
{        
    // if state space is size zero, we still treat it as size 1
    const int nstates1 = max(matrix->nstates1, 1);
    const int nstates2 = max(matrix->nstates2, 1);
    
    // initialize all entries in col2 to log(0)
    for (int k=0; k<nstates2; k++)
        col2[k] = -INFINITY;

    // add deterministic transitions
    for (int j=0; j<nstates1; j++) {
        int k = matrix->determ[j];
        if (j != matrix->recombsrc && j != matrix->recoalsrc && k != -1) {
            col2[k] = logadd(col2[k], col1[j] + matrix->determprob[j]);
        }
    }
    
    for (int k=0; k<nstates2; k++) {
        if (matrix->recombsrc != -1 && matrix->recombrow[k] > -INFINITY)
            col2[k] = logadd(col2[k], col1[matrix->recombsrc] + 
                             matrix->recombrow[k]);
        if (matrix->recoalsrc != -1 && matrix->recoalrow[k] > -INFINITY)
            col2[k] = logadd(col2[k], col1[matrix->recoalsrc] + 
                             matrix->recoalrow[k]);
        col2[k] += emit[k];
        assert(!isnan(col2[k]));
    }
    
    assert(max_array(col2, nstates2) > -INFINITY);
}



// run forward algorithm with matrices precomputed
void arghmm_forward_alg_fast(LocalTrees *trees, ArgModel *model,
    Sequences *sequences, ArgHmmMatrixIter *matrix_iter, 
    ArgHmmForwardTable *forward, bool prior_given, bool internal)
{
    LineageCounts lineages(model->ntimes);
    States states;
    ArgHmmMatrices matrices;
    
    double **fw = forward->get_table();

    // forward algorithm over local trees
    LocalTree *last_tree = NULL;
    for (matrix_iter->begin(); matrix_iter->more(); matrix_iter->next()) {
        LocalTrees::iterator it = matrix_iter->get_tree_iter();
        matrix_iter->get_matrices(&matrices);
        int pos = matrix_iter->get_position();
        LocalTree *tree = it->tree;

        // allocate the forward table
        if (pos > trees->start_coord || !prior_given)
            forward->new_block(pos, pos+matrices.blocklen, matrices.nstates2);
        
        if (internal)
            get_coal_states_internal(it->tree, model->ntimes, states);
        else
            get_coal_states(it->tree, model->ntimes, states);
        lineages.count(it->tree, internal);
        
        // use switch matrix for first column of forward table
        // if we have a previous state space (i.e. not first block)
        if (pos == trees->start_coord) {
            // calculate prior of first state
            if (!prior_given) {
                int minage = 0;
                if (internal) {
                    int subtree_root = it->tree->nodes[it->tree->root].child[0];
                    minage = it->tree->nodes[subtree_root].age;
                }
                calc_state_priors(states, &lineages, model, fw[pos], minage);
            }
        } else {
            // perform one column of forward algorithm with transmat_switch
            Spr &spr = it->spr;
            arghmm_forward_switch(fw[pos-1], fw[pos], 
                matrices.transmat_switch, matrices.emit[0]);
        }
        
        // calculate rest of block
        arghmm_forward_alg_block(it->tree, model, matrices.blocklen, 
                                 states, lineages, 
                                 matrices.transmat,
                                 matrices.emit, &fw[pos], internal);

        // safety check
        int nstates = max(matrices.transmat->nstates, 1);
        double top = max_array(fw[pos + matrices.blocklen - 1], nstates);
        assert(top > -INFINITY);

        last_tree = tree;
    }
}



// run forward algorithm with matrices precomputed
void arghmm_forward_alg(LocalTrees *trees, ArgModel *model,
    Sequences *sequences, ArgHmmMatrixIter *matrix_iter, 
    ArgHmmForwardTable *forward, bool prior_given)
{
    LineageCounts lineages(model->ntimes);
    States states;
    ArgHmmMatrices matrices;

    double **fw = forward->get_table();

    // forward algorithm over local trees
    for (matrix_iter->begin(); matrix_iter->more(); matrix_iter->next()) {
        LocalTrees::iterator it = matrix_iter->get_tree_iter();
        int pos = matrix_iter->get_position();
        matrix_iter->get_matrices(&matrices);

        assert(matrices.transmat != NULL &&
               matrices.transmat_switch != NULL);

        // allocate the forward table
        if (pos > trees->start_coord || !prior_given)
            forward->new_block(pos, pos+matrices.blocklen, matrices.nstates2);
        
        get_coal_states(it->tree, model->ntimes, states);
        lineages.count(it->tree);
        
        // use switch matrix for first column of forward table
        // if we have a previous state space (i.e. not first block)
        if (pos == trees->start_coord) {
            // calculate prior of first state
            if (!prior_given)
                calc_state_priors(states, &lineages, model, fw[pos]);
        } else {
            // perform one column of forward algorithm with transmat_switch
            forward_step(fw[pos-1], fw[pos], matrices.nstates1,
                         matrices.nstates2,
                         matrices.transprobs_switch, matrices.emit[0]);
        }

        // calculate rest of block
        forward_alg(matrices.blocklen, states.size(), matrices.transprobs,
                    matrices.emit, &fw[pos]);
    }
}



//=============================================================================
// Sample thread paths


void max_hmm_posterior(int n, const LocalTree *tree, const States &states,
                       const TransMatrix *matrix, const double *const *emit, 
                       const double *const *fw, int *path)
{
    // NOTE: path[n-1] must already be sampled
    
    const int nstates = states.size();
    double trans[nstates];
    int last_k = -1;

    // recurse
    for (int i=n-2; i>=0; i--) {
        int k = path[i+1];
        
        // recompute transition probabilities if state (k) changes
        if (k != last_k) {
            for (int j=0; j<nstates; j++)
                trans[j] = matrix->get_transition_prob(tree, states, j, k);
            last_k = k;
        }

        // find max transition
        int maxj = 0;
        double maxprob = fw[i][0] + trans[0];
        for (int j=1; j<nstates; j++) {
            double prob = fw[i][j] + trans[j];
            if (prob > maxprob) {
                maxj = j;
                maxprob = prob;
            }
        }
        path[i] = maxj;
    }
}


int max_hmm_posterior_step(const TransMatrixSwitch *matrix, 
                           const double *col1, int state2)
{
    const int nstates1 = matrix->nstates1;
    
    int maxj = 0;
    double maxprob = col1[0] + matrix->get_transition_prob(0, state2);
    for (int j=1; j<nstates1; j++) {
        double prob = col1[j] + matrix->get_transition_prob(j, state2);
        if (prob > maxprob) {
            maxj = j;
            maxprob = prob;
        }
    }
    return maxj;
}




void sample_hmm_posterior(int n, const LocalTree *tree, const States &states,
                          const TransMatrix *matrix, const double *const *emit, 
                          const double *const *fw, int *path)
{
    // NOTE: path[n-1] must already be sampled
    
    const int nstates = max(states.size(), (size_t)1);
    double A[nstates];
    double trans[nstates];
    int last_k = -1;

    // recurse
    for (int i=n-2; i>=0; i--) {
        int k = path[i+1];
        
        // recompute transition probabilities if state (k) changes
        if (k != last_k) {
            for (int j=0; j<nstates; j++)
                trans[j] = matrix->get_transition_prob(tree, states, j, k);
            last_k = k;
        }

        for (int j=0; j<nstates; j++)
            A[j] = fw[i][j] + trans[j];
        double total = logsum_fast(A, nstates);
        assert(total != -INFINITY);
        for (int j=0; j<nstates; j++)
            A[j] = t2exp(A[j] - total);
        path[i] = sample(A, nstates);

        // DEBUG
        assert(trans[path[i]] != -INFINITY);
    }
}


int sample_hmm_posterior_step(const TransMatrixSwitch *matrix, 
                              const double *col1, int state2)
{
    const int nstates1 = max(matrix->nstates1, 1);
    double A[nstates1];
    
    for (int j=0; j<nstates1; j++)
        A[j] = col1[j] + matrix->get_transition_prob(j, state2);
    const double total = logsum_fast(A, nstates1);
    assert(total != -INFINITY);
    for (int j=0; j<nstates1; j++)
        A[j] = t2exp(A[j] - total);
    int k = sample(A, nstates1);

    // DEBUG
    assert(matrix->get_transition_prob(k, state2) != -INFINITY);
    return k;
}


void stochastic_traceback_fast(LocalTrees *trees, ArgModel *model, 
                               ArgHmmMatrixIter *matrix_iter, 
                               double **fw, int *path, 
                               bool last_state_given, bool internal)
{
    ArgHmmMatrices mat;
    States states;

    // choose last column first
    matrix_iter->rbegin();    
    int pos = trees->end_coord;

    if (!last_state_given) {
        matrix_iter->get_matrices(&mat);
        int nstates = max(mat.nstates2, 1);
        double total = logsum_fast(fw[pos - 1], nstates);
        double vec[nstates];
        for (int i=0; i<nstates; i++)
            vec[i] = exp(fw[pos - 1][i] - total);
        path[pos - 1] = sample(vec, nstates);
    }
    
    // iterate backward through blocks
    for (; matrix_iter->more(); matrix_iter->prev()) {
        matrix_iter->get_matrices(&mat);
        LocalTree *tree = matrix_iter->get_tree_iter()->tree;
        if (internal)
            get_coal_states_internal(tree, model->ntimes, states);
        else
            get_coal_states(tree, model->ntimes, states);
        pos -= mat.blocklen;
        
        sample_hmm_posterior(mat.blocklen, tree, states,
                             mat.transmat, mat.emit, 
                             &fw[pos], &path[pos]);

        // use switch matrix for last col of next block
        if (pos > trees->start_coord) {
            int i = pos - 1;
            path[i] = sample_hmm_posterior_step(
                mat.transmat_switch, fw[i], path[i+1]);
        }
    }
}


void stochastic_traceback(ArgHmmMatrixIter *matrix_iter, 
                          double **fw, int *path, 
                          bool last_state_given)
{
    ArgHmmMatrices mat;

    // choose last column first
    matrix_iter->rbegin();
    matrix_iter->get_matrices(&mat);
    int start_coord = matrix_iter->get_position();
    int pos = start_coord + mat.blocklen;

    if (!last_state_given) {
        int nstates = mat.nstates2;
        double total = logsum_fast(fw[pos - 1], nstates);
        double vec[nstates];
        for (int i=0; i<nstates; i++)
            vec[i] = exp(fw[pos - 1][i] - total);
        path[pos - 1] = sample(vec, nstates);
    }
    
    // iterate backward through blocks
    for (; matrix_iter->more(); matrix_iter->prev()) {        
        matrix_iter->get_matrices(&mat);
        pos -= mat.blocklen;
        
        sample_hmm_posterior(mat.blocklen, mat.nstates2, mat.transprobs, 
                             &fw[pos], &path[pos]);
        
        // use switch matrix for last col of next block
        if (pos > start_coord) {
            int i = pos - 1;
            path[i] = sample_hmm_posterior_step(mat.nstates1, 
                                                mat.transprobs_switch, 
                                                fw[i], path[i+1]);
        }
    }
}


void max_traceback_fast(LocalTrees *trees, ArgModel *model, 
                        ArgHmmMatrixIter *matrix_iter, 
                        double **fw, int *path, 
                        bool last_state_given)
{
    ArgHmmMatrices mat;
    States states;

    // choose last column first
    matrix_iter->rbegin();
    int pos = trees->end_coord;

    if (!last_state_given) {
        matrix_iter->get_matrices(&mat);
        int nstates = mat.nstates2;
        int maxi = 0;
        double maxprob = fw[pos - 1][0];
        for (int i=1; i<nstates; i++) {
            if (fw[pos - 1][i] > maxprob) {
                maxi = i;
                maxprob = fw[pos - 1][i];
            }
        }
        path[pos - 1] = maxi;
    }
    
    // iterate backward through blocks
    for (; matrix_iter->more(); matrix_iter->prev()) {
        matrix_iter->get_matrices(&mat);
        LocalTree *tree = matrix_iter->get_tree_iter()->tree;
        get_coal_states(tree, model->ntimes, states);
        pos -= mat.blocklen;
        
        max_hmm_posterior(mat.blocklen, tree, states,
                          mat.transmat, mat.emit, 
                          &fw[pos], &path[pos]);

        // use switch matrix for last col of next block
        if (pos > trees->start_coord) {
            int i = pos - 1;
            path[i] = max_hmm_posterior_step(
                mat.transmat_switch, fw[i], path[i+1]);
        }
    }
}



//=============================================================================
// ARG sampling


// sample the thread of the last chromosome
void sample_arg_thread(ArgModel *model, Sequences *sequences, 
                       LocalTrees *trees, int new_chrom)
{
    // allocate temp variables
    ArgHmmForwardTable forward(trees->start_coord, trees->length());
    int *thread_path_alloc = new int [trees->length()];
    int *thread_path = &thread_path_alloc[-trees->start_coord];

    // build matrices
    ArgHmmMatrixIter matrix_iter(model, sequences, trees, new_chrom);
    
    // compute forward table
    Timer time;
    arghmm_forward_alg_fast(trees, model, sequences, &matrix_iter, &forward);
    int nstates = get_num_coal_states(trees->front().tree, model->ntimes);
    printTimerLog(time, LOG_QUIET, 
                  "forward (%3d states, %6d blocks):", 
                  nstates, trees->get_num_trees());

    // traceback
    time.start();
    double **fw = forward.get_table();
    ArgHmmMatrixIter matrix_iter2(model, NULL, trees, new_chrom);
    stochastic_traceback_fast(trees, model, &matrix_iter2, fw, thread_path);
    printTimerLog(time, LOG_QUIET, 
                  "trace:                              ");

    time.start();

    // sample recombination points
    vector<int> recomb_pos;
    vector<NodePoint> recombs;
    sample_recombinations(trees, model, &matrix_iter2,
                          thread_path, recomb_pos, recombs);

    // add thread to ARG
    add_arg_thread(trees, model->ntimes, thread_path, new_chrom, 
                   recomb_pos, recombs);

    //printf("add thread:  %e s\n", time.time());
    printTimerLog(time, LOG_QUIET, 
                  "add thread:                         ");

    // clean up
    delete [] thread_path_alloc;
}


// sample the thread of the last chromosome
void sample_arg_thread_internal(ArgModel *model, Sequences *sequences, 
                                LocalTrees *trees)
{
    const bool internal = true;

    // allocate temp variables
    ArgHmmForwardTable forward(trees->start_coord, trees->length());
    int *thread_path_alloc = new int [trees->length()];
    int *thread_path = &thread_path_alloc[-trees->start_coord];

    // build matrices
    ArgHmmMatrixIter matrix_iter(model, sequences, trees);
    matrix_iter.set_internal(internal);
    
    
    // compute forward table
    Timer time;
    arghmm_forward_alg_fast(trees, model, sequences, &matrix_iter, &forward,
                            false, internal);
    int nstates = get_num_coal_states_internal(
        trees->front().tree, model->ntimes);
    printTimerLog(time, LOG_QUIET, 
                  "forward (%3d states, %6d blocks):", 
                  nstates, trees->get_num_trees());

    // traceback
    time.start();
    double **fw = forward.get_table();
    ArgHmmMatrixIter matrix_iter2(model, NULL, trees);
    matrix_iter2.set_internal(internal);
    stochastic_traceback_fast(trees, model, &matrix_iter2, fw, thread_path,
                              false, internal);
    printTimerLog(time, LOG_QUIET, 
                  "trace:                              ");

    time.start();


    // sample recombination points
    vector<int> recomb_pos;
    vector<NodePoint> recombs;
    sample_recombinations(trees, model, &matrix_iter2,
                          thread_path, recomb_pos, recombs, internal);

    // add thread to ARG
    add_arg_thread_path(trees, model->ntimes, thread_path,
                        recomb_pos, recombs);

    //printf("add thread:  %e s\n", time.time());
    printTimerLog(time, LOG_QUIET, 
                  "add thread:                         ");

    // clean up
    delete [] thread_path_alloc;
}


// sample the thread of the last chromosome
void max_arg_thread(ArgModel *model, Sequences *sequences, 
                    LocalTrees *trees, int new_chrom)
{
    // allocate temp variables
    ArgHmmForwardTable forward(trees->start_coord, trees->length());
    double **fw = forward.get_table();
    int *thread_path_alloc = new int [trees->length()];
    int *thread_path = &thread_path_alloc[-trees->start_coord];

    // build matrices
    Timer time;
    ArgHmmMatrixList matrix_list(model, sequences, trees, new_chrom);
    matrix_list.setup();
    printf("matrix calc: %e s\n", time.time());
    
    // compute forward table
    time.start();
    arghmm_forward_alg_fast(trees, model, sequences, &matrix_list, &forward);
    int nstates = get_num_coal_states(trees->front().tree, model->ntimes);
    printf("forward:     %e s  (%d states, %d blocks)\n", time.time(),
           nstates, trees->get_num_trees());

    // traceback
    time.start();
    //stochastic_traceback(&matrix_list, fw, thread_path, sequences->seqlen);
    max_traceback_fast(trees, model, &matrix_list, fw, 
                       thread_path, sequences->length());
    printf("trace:       %e s\n", time.time());


    time.start();

    // sample recombination points
    vector<int> recomb_pos;
    vector<NodePoint> recombs;
    max_recombinations(trees, model, &matrix_list,
                       thread_path, recomb_pos, recombs);

    // add thread to ARG
    add_arg_thread(trees, model->ntimes, thread_path, new_chrom, 
                   recomb_pos, recombs);

    printf("add thread:  %e s\n", time.time());
    
    // clean up
    delete [] thread_path_alloc;
}


// sample the thread of the last chromosome, conditioned on a given
// start and end state
void cond_sample_arg_thread(ArgModel *model, Sequences *sequences, 
                            LocalTrees *trees, int new_chrom,
                            State start_state, State end_state)
{
    // allocate temp variables
    ArgHmmForwardTable forward(trees->start_coord, trees->length());
    States states;
    LocalTree *tree;
    double **fw = forward.get_table();
    int *thread_path_alloc = new int [trees->length()];
    int *thread_path = &thread_path_alloc[-trees->start_coord];

    // build matrices
    Timer time;
    ArgHmmMatrixList matrix_list(model, sequences, trees, new_chrom);
    matrix_list.setup();
    printf("matrix calc: %e s\n", time.time());
    
    // fill in first column of forward table
    matrix_list.begin();
    tree = matrix_list.get_tree_iter()->tree;
    get_coal_states(tree, model->ntimes, states);
    forward.new_block(trees->start_coord, trees->start_coord + 
                      trees->begin()->blocklen, states.size());
    bool found = false;
    for (unsigned int j=0; j<states.size(); j++) {
        if (states[j] == start_state) {
            fw[trees->start_coord][j] = 0.0;
            found = true;
        } else {
            fw[trees->start_coord][j] = -INFINITY;
        }
    }
    assert(found);

    // compute forward table
    time.start();
    arghmm_forward_alg_fast(trees, model, sequences, &matrix_list, &forward,
                            true);
    int nstates = get_num_coal_states(trees->front().tree, model->ntimes);
    printf("forward:     %e s  (%d states, %d blocks)\n", time.time(),
           nstates, trees->get_num_trees());

    // fill in last state of traceback
    matrix_list.rbegin();
    tree = matrix_list.get_tree_iter()->tree;
    get_coal_states(tree, model->ntimes, states);
    thread_path[trees->end_coord-1] = find_vector(states, end_state);
    assert(thread_path[trees->end_coord-1] != -1);

    // traceback
    time.start();
    stochastic_traceback_fast(trees, model, &matrix_list, fw, 
                              thread_path, true);
    printf("trace:       %e s\n", time.time());
    assert(fw[trees->start_coord][thread_path[trees->start_coord]] == 0.0);




    // sample recombination points
    time.start();
    vector<int> recomb_pos;
    vector<NodePoint> recombs;
    sample_recombinations(trees, model, &matrix_list,
                          thread_path, recomb_pos, recombs);

    // add thread to ARG
    add_arg_thread(trees, model->ntimes, thread_path, new_chrom, 
                   recomb_pos, recombs);

    printf("add thread:  %e s\n", time.time());

    // clean up
    delete [] thread_path_alloc;
}



// resample the threading of one chromosome
void resample_arg_thread(ArgModel *model, Sequences *sequences, 
                         LocalTrees *trees, int chrom)
{
    // remove chromosome from ARG and resample its thread
    remove_arg_thread(trees, chrom);
    sample_arg_thread(model, sequences, trees, chrom);
}



//=============================================================================
// C interface

extern "C" {


// perform forward algorithm
double **arghmm_forward_alg(
    LocalTrees *trees, double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, bool prior_given, double *prior,
    bool internal)
{    
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);

    // build matrices
    ArgHmmMatrixList matrix_list(&model, &sequences, trees);
    matrix_list.set_internal(internal);
    matrix_list.setup();

    ArgHmmForwardTableOld forward(0, sequences.length());

    // setup prior
    if (prior_given) {
        LocalTree *tree = trees->begin()->tree;
        int blocklen = trees->begin()->blocklen;
        LineageCounts lineages(ntimes);
        States states;
        if (internal)
            get_coal_states_internal(tree, model.ntimes, states);
        else
            get_coal_states(tree, model.ntimes, states);

        forward.new_block(0, blocklen, states.size());
        double **fw = forward.get_table();
        for (unsigned int i=0; i<states.size(); i++)
            fw[0][i] = prior[i];
    }

    arghmm_forward_alg_fast(trees, &model, &sequences, &matrix_list,
                            &forward, prior_given, internal);

    // steal pointer
    double **fw = forward.detach_table();

    return fw;
}


// perform forward algorithm and sample threading path from posterior
intstate *arghmm_sample_posterior(
    int **ptrees, int **ages, int **sprs, int *blocklens,
    int ntrees, int nnodes, double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, intstate *path=NULL)
{    
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    LocalTrees trees(ptrees, ages, sprs, blocklens, ntrees, nnodes);
    Sequences sequences(seqs, nseqs, seqlen);
    
    // build matrices
    ArgHmmMatrixList matrix_list(&model, &sequences, &trees, -1, true);
    matrix_list.setup();
    
    // compute forward table
    ArgHmmForwardTable forward(0, seqlen);
    arghmm_forward_alg_fast(&trees, &model, &sequences, &matrix_list, &forward);

    // traceback
    int *ipath = new int [seqlen];
    stochastic_traceback(&matrix_list, forward.get_table(), ipath, seqlen);

    
    // convert path
    if (path == NULL)
        path = new intstate [seqlen];

    States states;
    int end = trees.start_coord;
    for (LocalTrees::iterator it=trees.begin(); it != trees.end(); ++it) {
        int start = end;
        int end = start + it->blocklen;
        get_coal_states(it->tree, ntimes, states);

        for (int i=start; i<end; i++) {
            int istate = ipath[i];
            path[i][0] = states[istate].node;
            path[i][1] = states[istate].time;
        }
    }


    // clean up
    delete [] ipath;

    return path;
}


// sample the thread of an internal branch
void arghmm_sample_arg_thread_internal(LocalTrees *trees,
    double *times, int ntimes, double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, int *thread_path)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);    
    const bool internal = true;

    // allocate temp variables
    ArgHmmForwardTable forward(trees->start_coord, trees->length());
    thread_path = &thread_path[-trees->start_coord];

    // forward algorithm
    ArgHmmMatrixIter matrix_iter(&model, &sequences, trees);
    matrix_iter.set_internal(internal);
    arghmm_forward_alg_fast(trees, &model, &sequences, &matrix_iter, &forward,
                            false, internal);

    // traceback
    double **fw = forward.get_table();
    ArgHmmMatrixIter matrix_iter2(&model, NULL, trees);
    matrix_iter2.set_internal(internal);
    stochastic_traceback_fast(trees, &model, &matrix_iter2, fw, thread_path,
                              false, internal);
}


// add one chromosome to an ARG
LocalTrees *arghmm_sample_thread(
    LocalTrees *trees, double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    int new_chrom = nseqs -  1;
    
    sample_arg_thread(&model, &sequences, trees, new_chrom);
    
    return trees;
}


// add one chromosome to an ARG using maximization
LocalTrees *arghmm_max_thread(
    LocalTrees *trees, double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    int new_chrom = nseqs -  1;
    
    max_arg_thread(&model, &sequences, trees, new_chrom);

    return trees;
}


void delete_path(int *path)
{
    delete [] path;
}


void delete_double_matrix(double **mat, int nrows)
{
    delete_matrix<double>(mat, nrows);
}

void delete_forward_matrix(double **mat, int nrows)
{
    for (int i=0; i<nrows; i++)
        delete [] mat[i];
    delete [] mat;
}



} // extern "C"

} // namespace arghmm

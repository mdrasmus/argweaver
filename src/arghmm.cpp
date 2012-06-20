// c++ includes
#include <list>
#include <vector>
#include <string.h>

// arghmm includes
#include "common.h"
#include "emit.h"
#include "hmm.h"
#include "itree.h"
#include "local_tree.h"
#include "logging.h"
#include "matrices.h"
#include "model.h"
#include "ptree.h"
#include "recomb.h"
#include "seq.h"
#include "states.h"
#include "thread.h"
#include "trans.h"


using namespace std;

using namespace spidir;
using namespace dlcoal;


namespace arghmm {



//=============================================================================
// Forward tables


class ArgHmmForwardTable
{
public:
    ArgHmmForwardTable(int seqlen) :
        seqlen(seqlen)
    {
        fw = new double *[seqlen];
    }

    virtual ~ArgHmmForwardTable()
    {
        delete_blocks();
        if (fw) {
            delete [] fw;
            fw = NULL;
        }
    }


    // allocate another block of the forward table
    virtual void new_block(int start, int end, int nstates)
    {
        // allocate block
        int blocklen = end - start;
        double *block = new double [blocklen * nstates];
        blocks.push_back(block);

        // link block to fw table
        for (int i=start; i<end; i++)
            fw[i] = &block[(i-start)*nstates];
    }

    // delete all blocks
    virtual void delete_blocks()
    {
        for (unsigned int i=0; i<blocks.size(); i++)
            delete [] blocks[i];
        blocks.clear();
    }

    int seqlen;
    double **fw;

protected:
    vector<double*> blocks;
};


// older style allocation for testing with python
class ArgHmmForwardTableOld : public ArgHmmForwardTable
{
public:
    ArgHmmForwardTableOld(int seqlen) :
        ArgHmmForwardTable(seqlen)
    {
        fw = new double *[seqlen];
    }

    virtual ~ArgHmmForwardTableOld() {}

    // allocate another block of the forward table
    virtual void new_block(int start, int end, int nstates)
    {
        // allocate block
        for (int i=start; i<end; i++)
            fw[i] = new double [nstates];
    }
};


//=============================================================================
// Forward algorithm for thread path


// compute one block of forward algorithm with compressed transition matrices
void arghmm_forward_alg_block(const LocalTree *tree, const ArgModel *model,
                              int blocklen, const States &states, 
                              const LineageCounts &lineages,
                              const TransMatrixCompress *matrix,
                              const double* const *emit, double **fw)
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

    // compute ntimes*ntimes and ntime*nstates temp matrices
    double lnorecombs[ntimes];
    double tmatrix[ntimes][ntimes];
    double tmatrix2[ntimes][nstates];
    for (int a=0; a<ntimes-1; a++) {
        lnorecombs[a] = log(norecombs[a]);

        for (int b=0; b<ntimes-1; b++) {
            const double I = double(a <= b);
            tmatrix[a][b] = log(D[a] * E[b] * (B[min(a,b)] - I * G[a]));
        }

        for (int k=0; k<nstates; k++) {
            const int b = states[k].time;
            const int node2 = states[k].node;
            const int c = nodes[node2].age;
            const double Bc = (c > 0 ? B[c-1] : 0.0);
            const double I = float(a <= b);
            tmatrix2[a][k] = log(D[a] * E[b] * (B[min(a,b)] - I * G[a] - Bc));
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
        ages1[i] = nodes[i].age;
        ages2[i] = (tree->root == i) ? maxtime : nodes[nodes[i].parent].age;
    }
    

    // group states by age
    int nstates_per_time[ntimes];
    memset(nstates_per_time, 0, ntimes*sizeof(int));
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

    double vec[nstates];
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
            fgroups[a] = logsum(&vec[offset[a]], nstates_per_time[a]);
        
        // multiply tmatrix and fgroups together
        for (int b=0; b<ntimes-1; b++) {
            for (int a=0; a<ntimes-1; a++)
                vec[a] = tmatrix[a][b] + fgroups[a];
            tmatrix_fgroups[b] = logsum(vec, ntimes-1);
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
            int j2 = state_lookup.lookup(node2, age2);
            for (int j=state_lookup.lookup(node2, age1), a=age1; j<=j2; j++,a++)
                vec[m++] = tmatrix2[a][k] + col1[j];
            
            // same state case (add possibility of no recomb)
            vec[m++] = lnorecombs[b] + col1[k];
            
            col2[k] = logsum(vec, m) + emit2[k];
        }
    }
}


// run forward algorithm for one column of the table
// use switch matrix
void arghmm_forward_switch(const double *col1, double* col2, 
                           const TransMatrixSwitchCompress *matrix,
                           const double *emit)
{        
    // initialize all entries in col2 to log(0)
    for (int k=0; k<matrix->nstates2; k++)
        col2[k] = -INFINITY;

    // add deterministic transitions
    for (int j=0; j<matrix->nstates1; j++) {
        if (j != matrix->recombsrc && j != matrix->recoalsrc) {
            int k = matrix->determ[j];
            col2[k] = logadd(col2[k], col1[j] + matrix->determprob[j]);
        }
    }
    
    for (int k=0; k<matrix->nstates2; k++) {
        if (matrix->recombrow[k] > -INFINITY)
            col2[k] = logadd(col2[k], col1[matrix->recombsrc] + 
                             matrix->recombrow[k]);
        if (matrix->recoalrow[k] > -INFINITY)
            col2[k] = logadd(col2[k], col1[matrix->recoalsrc] + 
                             matrix->recoalrow[k]);
        col2[k] += emit[k];
    }
}



// run forward algorithm with matrices precomputed
void arghmm_forward_alg_fast(LocalTrees *trees, ArgModel *model,
    Sequences *sequences, ArgHmmMatrixIter *matrix_iter, 
       ArgHmmForwardTable *forward, bool prior_given=false)
{
    LineageCounts lineages(model->ntimes);
    States states;
    ArgHmmMatrices matrices;

    double **fw = forward->fw;

    // forward algorithm over local trees
    for (matrix_iter->begin(); matrix_iter->more(); matrix_iter->next()) {
        LocalTrees::iterator it = matrix_iter->get_tree_iter();
        matrix_iter->get_matrices(&matrices);
        int pos = matrix_iter->get_position();

        // allocate the forward table
        forward->new_block(pos, pos+matrices.blocklen, matrices.nstates2);
        
        get_coal_states(it->tree, model->ntimes, states);
        lineages.count(it->tree);
        
        // use switch matrix for first column of forward table
        // if we have a previous state space (i.e. not first block)
        if (pos == 0) {
            // calculate prior of first state
            if (!prior_given)
                calc_state_priors(states, &lineages, model, fw[0]);
        } else {
            // perform one column of forward algorithm with transmat_switch
            arghmm_forward_switch(fw[pos-1], fw[pos], 
                matrices.transmat_switch_compress, matrices.emit[0]);
        }

        // calculate rest of block
        arghmm_forward_alg_block(it->tree, model, matrices.blocklen, 
                                 states, lineages, 
                                 matrices.transmat_compress,
                                 matrices.emit, &fw[pos]);
    }
}



// run forward algorithm with matrices precomputed
void arghmm_forward_alg(LocalTrees *trees, ArgModel *model,
    Sequences *sequences, ArgHmmMatrixIter *matrix_iter, 
    ArgHmmForwardTable *forward, bool prior_given=false)
{
    LineageCounts lineages(model->ntimes);
    States states;
    ArgHmmMatrices matrices;

    double **fw = forward->fw;

    // forward algorithm over local trees
    for (matrix_iter->begin(); matrix_iter->more(); matrix_iter->next()) {
        LocalTrees::iterator it = matrix_iter->get_tree_iter();
        int pos = matrix_iter->get_position();
        matrix_iter->get_matrices(&matrices);

        assert(matrices.transmat != NULL &&
               matrices.transmat_switch != NULL);

        // allocate the forward table
        forward->new_block(pos, pos+matrices.blocklen, matrices.nstates2);
        
        get_coal_states(it->tree, model->ntimes, states);
        lineages.count(it->tree);
        
        // use switch matrix for first column of forward table
        // if we have a previous state space (i.e. not first block)
        if (pos == 0) {
            // calculate prior of first state
            if (!prior_given)
                calc_state_priors(states, &lineages, model, fw[0]);
        } else {
            // perform one column of forward algorithm with transmat_switch
            forward_step(fw[pos-1], fw[pos], matrices.nstates1,
                         matrices.nstates2,
                         matrices.transmat_switch, matrices.emit[0]);
        }

        // calculate rest of block
        forward_alg(matrices.blocklen, states.size(), matrices.transmat,
                    matrices.emit, &fw[pos]);
    }
}



//=============================================================================
// Sample thread paths


void max_hmm_posterior(int n, LocalTree *tree, const States &states,
                       TransMatrixCompress *matrix, double **emit, 
                       double **fw, int *path)
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


int max_hmm_posterior_step(TransMatrixSwitchCompress *matrix, 
                           double *col1, int state2)
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




void sample_hmm_posterior(int n, LocalTree *tree, const States &states,
                          TransMatrixCompress *matrix, double **emit, 
                          double **fw, int *path)
{
    // NOTE: path[n-1] must already be sampled
    
    const int nstates = states.size();
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
        double total = logsum(A, nstates);
        for (int j=0; j<nstates; j++)
            A[j] = exp(A[j] - total);
        path[i] = sample(A, nstates);
    }
}


int sample_hmm_posterior_step(TransMatrixSwitchCompress *matrix, 
                              double *col1, int state2)
{
    const int nstates1 = matrix->nstates1;
    double A[nstates1];
    
    for (int j=0; j<nstates1; j++)
        A[j] = col1[j] + matrix->get_transition_prob(j, state2);
    double total = logsum(A, nstates1);
    for (int j=0; j<nstates1; j++)
        A[j] = exp(A[j] - total);
    return sample(A, nstates1);
}


void stochastic_traceback_fast(LocalTrees *trees, ArgModel *model, 
                               ArgHmmMatrixIter *matrix_iter, 
                               double **fw, int *path, int seqlen,
                               bool last_state_given=false)
{
    ArgHmmMatrices mat;
    States states;

    // choose last column first
    matrix_iter->rbegin();    
    
    if (!last_state_given) {
        matrix_iter->get_matrices(&mat);
        int nstates = mat.nstates2;
        double total = logsum(fw[seqlen - 1], nstates);
        double vec[nstates];
        for (int i=0; i<nstates; i++)
            vec[i] = exp(fw[seqlen - 1][i] - total);
        path[seqlen - 1] = sample(vec, nstates);
    }
    
    // iterate backward through blocks
    int pos = seqlen;
    for (; matrix_iter->more(); matrix_iter->prev()) {
        matrix_iter->get_matrices(&mat);
        LocalTree *tree = matrix_iter->get_tree_iter()->tree;
        get_coal_states(tree, model->ntimes, states);
        pos -= mat.blocklen;
        
        sample_hmm_posterior(mat.blocklen, tree, states,
                             mat.transmat_compress, mat.emit, 
                             &fw[pos], &path[pos]);

        // use switch matrix for last col of next block
        if (pos > 0) {
            int i = pos - 1;
            path[i] = sample_hmm_posterior_step(
                mat.transmat_switch_compress, fw[i], path[i+1]);
        }
    }
}


void stochastic_traceback(ArgHmmMatrixIter *matrix_iter, 
                          double **fw, int *path, int seqlen,
                          bool last_state_given=false)
{
    ArgHmmMatrices mat;

    // choose last column first
    matrix_iter->rbegin();

    if (!last_state_given) {
        matrix_iter->get_matrices(&mat);
        int nstates = mat.nstates2;
        double total = logsum(fw[seqlen - 1], nstates);
        double vec[nstates];
        for (int i=0; i<nstates; i++)
            vec[i] = exp(fw[seqlen - 1][i] - total);
        path[seqlen - 1] = sample(vec, nstates);
    }
    
    // iterate backward through blocks
    int pos = seqlen;
    for (; matrix_iter->more(); matrix_iter->prev()) {        
        matrix_iter->get_matrices(&mat);
        pos -= mat.blocklen;
        
        sample_hmm_posterior(mat.blocklen, mat.nstates2, mat.transmat, 
                             &fw[pos], &path[pos]);
        
        // use switch matrix for last col of next block
        if (pos > 0) {
            int i = pos - 1;
            path[i] = sample_hmm_posterior_step(mat.nstates1, 
                                                mat.transmat_switch, 
                                                fw[i], path[i+1]);
        }
    }
}


void max_traceback_fast(LocalTrees *trees, ArgModel *model, 
                        ArgHmmMatrixIter *matrix_iter, 
                        double **fw, int *path, int seqlen,
                        bool last_state_given=false)
{
    ArgHmmMatrices mat;
    States states;

    // choose last column first
    matrix_iter->rbegin();

    if (!last_state_given) {
        matrix_iter->get_matrices(&mat);
        int nstates = mat.nstates2;
        int maxi = 0;
        double maxprob = fw[seqlen - 1][0];
        for (int i=1; i<nstates; i++) {
            if (fw[seqlen - 1][i] > maxprob) {
                maxi = i;
                maxprob = fw[seqlen - 1][i];
            }
        }
        path[seqlen - 1] = maxi;
    }
    
    // iterate backward through blocks
    int pos = seqlen;
    for (; matrix_iter->more(); matrix_iter->prev()) {
        matrix_iter->get_matrices(&mat);
        LocalTree *tree = matrix_iter->get_tree_iter()->tree;
        get_coal_states(tree, model->ntimes, states);
        pos -= mat.blocklen;
        
        max_hmm_posterior(mat.blocklen, tree, states,
                          mat.transmat_compress, mat.emit, 
                          &fw[pos], &path[pos]);

        // use switch matrix for last col of next block
        if (pos > 0) {
            int i = pos - 1;
            path[i] = max_hmm_posterior_step(
                mat.transmat_switch_compress, fw[i], path[i+1]);
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
    ArgHmmForwardTable forward(sequences->seqlen);
    double **fw = forward.fw;
    int *thread_path = new int [sequences->seqlen];

    // build matrices
    Timer time;
    ArgHmmMatrixList matrix_list(model, sequences, trees, new_chrom, false);
    matrix_list.setup();
    printf("matrix calc: %e s\n", time.time());
    
    // compute forward table
    time.start();
    arghmm_forward_alg_fast(trees, model, sequences, &matrix_list, &forward);
    printf("forward:     %e s  (%d states, %d blocks)\n", time.time(),
           matrix_list.matrices[0].nstates2,
           (int) matrix_list.matrices.size());


    // traceback
    time.start();
    stochastic_traceback_fast(trees, model, &matrix_list, fw, 
                              thread_path, sequences->seqlen);
    printf("trace:       %e s\n", time.time());


    time.start();

    // sample recombination points
    vector<int> recomb_pos;
    vector<NodePoint> recombs;
    sample_recombinations(trees, model, &matrix_list,
                          thread_path, recomb_pos, recombs);

    // add thread to ARG
    add_arg_thread(trees, model->ntimes, thread_path, new_chrom, 
                   recomb_pos, recombs);

    printf("add thread:  %e s\n", time.time());

    // clean up
    delete [] thread_path;
}


// resample the threading of one chromosome
void resample_arg_thread(ArgModel *model, Sequences *sequences, 
                         LocalTrees *trees, int chrom)
{
    // remove chromosome from ARG and resample its thread
    remove_arg_thread(trees, chrom);
    sample_arg_thread(model, sequences, trees, chrom);
}


// sample the thread of the last chromosome
void max_arg_thread(ArgModel *model, Sequences *sequences, 
                    LocalTrees *trees, int new_chrom)
{
    // allocate temp variables
    ArgHmmForwardTable forward(sequences->seqlen);
    double **fw = forward.fw;
    int *thread_path = new int [sequences->seqlen];

    // build matrices
    Timer time;
    ArgHmmMatrixList matrix_list(model, sequences, trees, new_chrom, false);
    matrix_list.setup();
    printf("matrix calc: %e s\n", time.time());
    
    // compute forward table
    time.start();
    arghmm_forward_alg_fast(trees, model, sequences, &matrix_list, &forward);
    printf("forward:     %e s  (%d states, %d blocks)\n", time.time(),
           matrix_list.matrices[0].nstates2,
           (int) matrix_list.matrices.size());


    // traceback
    time.start();
    //stochastic_traceback(&matrix_list, fw, thread_path, sequences->seqlen);
    max_traceback_fast(trees, model, &matrix_list, fw, 
                       thread_path, sequences->seqlen);
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
    delete [] thread_path;
}


// sample the thread of the last chromosome, conditioned on a given
// start and end state
void cond_sample_arg_thread(ArgModel *model, Sequences *sequences, 
                            LocalTrees *trees, int new_chrom,
                            State start_state, State end_state)
{
    // allocate temp variables
    ArgHmmForwardTable forward(sequences->seqlen);
    States states;
    LocalTree *tree;
    double **fw = forward.fw;
    int *thread_path = new int [sequences->seqlen];

    // build matrices
    Timer time;
    ArgHmmMatrixList matrix_list(model, sequences, trees, new_chrom, false);
    matrix_list.setup();
    printf("matrix calc: %e s\n", time.time());
    
    // fill in first column of forward table
    matrix_list.begin();
    tree = matrix_list.get_tree_iter()->tree;
    get_coal_states(tree, model->ntimes, states);
    for (unsigned int j=0; j<states.size(); j++) {
        if (states[j] == start_state)
            fw[0][j] = 0.0;
        else
            fw[0][j] = -INFINITY;
    }

    // compute forward table
    time.start();
    arghmm_forward_alg_fast(trees, model, sequences, &matrix_list, &forward,
                            true);
    printf("forward:     %e s  (%d states, %d blocks)\n", time.time(),
           matrix_list.matrices[0].nstates2,
           (int) matrix_list.matrices.size());

    // fill in last state of traceback
    matrix_list.begin();
    tree = matrix_list.get_tree_iter()->tree;
    get_coal_states(tree, model->ntimes, states);
    thread_path[sequences->seqlen-1] = -1;
    for (unsigned int j=0; j<states.size(); j++) {
        if (states[j] == end_state) {
            thread_path[sequences->seqlen-1] = j;
            break;
        }
    }
    assert(thread_path[sequences->seqlen-1] != -1);

    // traceback
    time.start();
    stochastic_traceback_fast(trees, model, &matrix_list, fw, 
                              thread_path, sequences->seqlen, true);
    printf("trace:       %e s\n", time.time());


    time.start();

    // sample recombination points
    vector<int> recomb_pos;
    vector<NodePoint> recombs;
    sample_recombinations(trees, model, &matrix_list,
                          thread_path, recomb_pos, recombs);

    // add thread to ARG
    add_arg_thread(trees, model->ntimes, thread_path, new_chrom, 
                   recomb_pos, recombs);

    printf("add thread:  %e s\n", time.time());

    // clean up
    delete [] thread_path;
}


//=============================================================================
// sample full ARGs

// sequentially sample an ARG from scratch
// sequences are sampled in the order given
void sample_arg_seq(ArgModel *model, Sequences *sequences, LocalTrees *trees)
{
    const int seqlen = sequences->seqlen;

    // initialize ARG as trunk
    const int capacity = 2 * sequences->nseqs - 1;
    trees->make_trunk(0, seqlen, capacity);
    
    // add more chromosomes one by one
    for (int nchroms=2; nchroms<=sequences->nseqs; nchroms++) {
        // use first nchroms sequences
        Sequences sequences2(sequences->seqs, nchroms, seqlen);
        int new_chrom = nchroms - 1;
        sample_arg_thread(model, &sequences2, trees, new_chrom);
    }
}


// resample the threading of all the chromosomes
void resample_arg(ArgModel *model, Sequences *sequences, LocalTrees *trees,
                  int nremove=1)
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


// resample the threading of all the chromosomes
void remax_arg(ArgModel *model, Sequences *sequences, LocalTrees *trees,
               int nremove=1)
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


State find_state_sub_tree(LocalTree *full_tree, vector<int> full_seqids,
                          LocalTree *partial_tree, vector<int> partial_seqids,
                          int new_chrom)
{
    // reconcile partial_tree to full tree
    int recon[full_tree->nnodes];
    
    // reconcile leaves
    for (unsigned int i=0; i<full_seqids.size(); i++) {
        int seqid = full_seqids[i];
        recon[i] = -1;        
        for (unsigned int j=0; i<partial_seqids.size(); i++) {
            if (partial_seqids[j] == seqid) {
                recon[i] = j;
                break;
            }
        }
    }

    // postorder iterate over full tree to reconcile internal nodes
    int order[full_tree->nnodes];
    full_tree->get_postorder(order);
    LocalNode *nodes = full_tree->nodes;
    for (int i=0; i<full_tree->nnodes; i++) {
        int j = order[i];
        const int *child = nodes[j].child;
        
        if (!nodes[j].is_leaf()) {
            if (recon[child[0]] != -1) {
                if (recon[child[1]] != -1) {
                    // both children recon, so we recon to their LCA
                    recon[j] = partial_tree->nodes[recon[child[0]]].parent;
                    assert(partial_tree->nodes[recon[child[0]]].parent ==
                           partial_tree->nodes[recon[child[1]]].parent);
                } else {
                    // single child recons, copy its mapping
                    recon[j] = recon[child[0]];
                }
            } else {
                if (recon[child[1]] != -1) {
                    // single child recons, copy its mapping
                    recon[j] = recon[child[0]];
                } else {
                    // neither child recons, so neither do we
                    recon[j] = -1;
                }
            }
        }
    }

    // find new chrom in full_tree
    int ptr = -1;
    for (unsigned int i=0; i<full_seqids.size(); i++) {
        if (full_seqids[i] == new_chrom) {
            ptr = i;
            break;
        }
    }
    assert(ptr != -1);

    // walk up from new chrom until we hit a reconciled portion of full tree
    while (recon[ptr] == -1)
        ptr = nodes[ptr].parent;

    return State(recon[ptr], nodes[ptr].age);
}


// sequentially sample an ARG from scratch
// sequences are sampled in the order given
void cond_sample_arg_seq(ArgModel *model, Sequences *sequences, 
                         LocalTrees *trees, 
                         LocalTree *start_tree, LocalTree *end_tree,
                         vector<int> full_seqids)
{
    const int seqlen = sequences->seqlen;

    // initialize ARG as trunk
    const int capacity = 2 * sequences->nseqs - 1;
    trees->make_trunk(0, seqlen, capacity);
    
    // add more chromosomes one by one
    for (int nchroms=2; nchroms<=sequences->nseqs; nchroms++) {
        // use first nchroms sequences
        Sequences sequences2(sequences->seqs, nchroms, seqlen);
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
    }
}




//=============================================================================
// C interface
extern "C" {


double **arghmm_forward_alg(
    int **ptrees, int **ages, int **sprs, int *blocklens,
    int ntrees, int nnodes, double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen)
{    
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    LocalTrees trees(ptrees, ages, sprs, blocklens, ntrees, nnodes);
    Sequences sequences(seqs, nseqs, seqlen);

    // build matrices
    ArgHmmMatrixList matrix_list(&model, &sequences, &trees);
    matrix_list.setup();

    ArgHmmForwardTableOld forward(sequences.seqlen);
    arghmm_forward_alg_fast(&trees, &model, &sequences, &matrix_list,
                            &forward);

    // steal pointer
    double **fw = forward.fw;
    forward.fw = NULL;

    return fw;
}



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
    ArgHmmMatrixList matrix_list(&model, &sequences, &trees);
    matrix_list.setup();
    
    // compute forward table
    ArgHmmForwardTable forward(seqlen);
    arghmm_forward_alg_fast(&trees, &model, &sequences, &matrix_list, &forward);

    // traceback
    int *ipath = new int [seqlen];
    stochastic_traceback(&matrix_list, forward.fw, ipath, seqlen);

    
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

//----------------------
// sample ARGs


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




void delete_path(int *path)
{
    delete [] path;
}


void delete_double_matrix(double **mat, int nrows)
{
    delete_matrix<double>(mat, nrows);
}





} // extern C

}

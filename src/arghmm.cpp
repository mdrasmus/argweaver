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
#include "model.h"
#include "ptree.h"
#include "seq.h"
#include "states.h"
#include "thread.h"
#include "trans.h"


using namespace std;

using namespace spidir;
using namespace dlcoal;


namespace arghmm {


class ArgHmmMatrices
{
public:
    ArgHmmMatrices() :
        nstates1(0),
        nstates2(0),
        blocklen(0),
        transmat_compress(NULL),
        transmat(NULL),
        transmat_switch(NULL),
        emit(NULL)
    {}

    ArgHmmMatrices(int nstates1, int nstates2, int blocklen,
                   double **transmat, double **transmat_switch, double **emit):
        nstates1(nstates1),
        nstates2(nstates2),
        blocklen(blocklen),
        transmat_compress(NULL),
        transmat(transmat),
        transmat_switch(transmat_switch),
        emit(emit)
    {}

    ArgHmmMatrices(int nstates1, int nstates2, int blocklen,
                   TransMatrixCompress *transmat_compress,
                   double **transmat, double **transmat_switch, double **emit):
        nstates1(nstates1),
        nstates2(nstates2),
        blocklen(blocklen),
        transmat_compress(transmat_compress),
        transmat(transmat),
        transmat_switch(transmat_switch),
        emit(emit)
    {}

    ~ArgHmmMatrices() 
    {}
    
    // delete all matrices
    void clear()
    {
        if (transmat_compress) {
            delete transmat_compress;
            transmat_compress = NULL;
        }
        if (transmat) {
            delete_matrix<double>(transmat, nstates2);
            transmat = NULL;
        }
        if (transmat_switch) {
            delete_matrix<double>(transmat_switch, nstates1);
            transmat_switch = NULL;
        }
        if (emit) {
            delete_matrix<double>(emit, blocklen);
            emit = NULL;
        }
    }

    int nstates1;
    int nstates2;
    int blocklen;
    TransMatrixCompress* transmat_compress;
    double **transmat;
    double **transmat_switch;
    double **emit;
};


// iterates through matricies for the ArgHmm
class ArgHmmMatrixIter
{
public:
    ArgHmmMatrixIter(ArgModel *model, Sequences *seqs, LocalTrees *trees, 
                     int _new_chrom=-1, bool calc_transmat=true) :
        model(model),
        seqs(seqs),
        trees(trees),
        new_chrom(_new_chrom),
        calc_transmat(calc_transmat),
        lineages(model->ntimes)
    {
        if (new_chrom == -1)
            new_chrom = trees->get_num_leaves();
    }

    virtual ~ArgHmmMatrixIter() {}

    // initializes iterator
    virtual void begin()
    {
        begin(trees->begin());
    }
    
    virtual void begin(LocalTrees::iterator start)
    {
        tree_iter = start;
        
        // setup last_tree information
        if (start == trees->begin()) {
            // no last tree
            last_tree = NULL;
            last_states = NULL;
            states = &states1;
            
        } else {
            LocalTrees::iterator tree_iter2 = tree_iter;
            --tree_iter2;
            last_states = &states2;
            last_tree = tree_iter2->tree;
            states = &states1;
            get_coal_states(last_tree, model->ntimes, *last_states);
        }
    }
    
    
    virtual void rbegin()
    {
        LocalTrees::iterator it = --trees->end();
        rbegin(it);
    }
    
    virtual void rbegin(LocalTrees::iterator start)
    {
        begin(start);
        last_states = &states2;
    }


    virtual bool next()
    {
        // update pointers
        last_tree = tree_iter->tree;
        last_states = states;
        states = ((states == &states1) ? &states2 : &states1);

        ++tree_iter;
        return tree_iter != trees->end();
    }


    // moves iterator to previous block
    virtual bool prev()
    {
        --tree_iter;
        LocalTrees::iterator it = tree_iter;
        --it;

        if (it != trees->end()) {
            last_tree = it->tree;
            get_coal_states(last_tree, model->ntimes, *last_states);
        } else {
            last_tree = NULL;
        }
        
        return tree_iter != trees->end();
    }
    

    // returns true if there are more blocks
    virtual bool more()
    {
        return tree_iter != trees->end();
    }


    virtual void get_matrices(ArgHmmMatrices *matrices)
    {
        calc_matrices(matrices);
    }

    LocalTrees::iterator get_tree_iter()
    {
        return tree_iter;
    }

    
    // calculate transition and emission matrices for current block
    void calc_matrices(ArgHmmMatrices *matrices)
    {
        int pos = tree_iter->block.start;
        int blocklen = tree_iter->block.end - tree_iter->block.start;
        LocalTree *tree = tree_iter->tree;
        get_coal_states(tree, model->ntimes, *states);
        int nstates = states->size();
        matrices->blocklen = blocklen;
                
        // calculate emissions
        if (seqs) {
            char *subseqs[seqs->nseqs];
            for (int i=0; i<seqs->nseqs-1; i++)
                subseqs[i] = &seqs->seqs[trees->seqids[i]][pos];
            subseqs[seqs->nseqs-1] = &seqs->seqs[new_chrom][pos];
            matrices->emit = new_matrix<double>(blocklen, nstates);
            calc_emissions(*states, tree, subseqs, seqs->nseqs, blocklen, 
                           model, matrices->emit);
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
            lineages.count(last_tree);
                
            // calculate transmat_switch
            matrices->transmat_switch = new_matrix<double>(
                matrices->nstates1, matrices->nstates2);
            calc_transition_probs_switch(
                tree, last_tree, tree_iter->spr, tree_iter->mapping,
                *last_states, *states,
                model, &lineages, matrices->transmat_switch);
        }

        // update lineages to current tree
        lineages.count(tree);
        
        // calculate transmat and use it for rest of block
        matrices->transmat_compress = new TransMatrixCompress(
            model->ntimes, nstates);
        calc_transition_probs_compress(tree, model, *states, &lineages, 
                                       matrices->transmat_compress);
        if (calc_transmat) {
            matrices->transmat = new_matrix<double>(nstates, nstates);
            calc_transition_probs(tree, model, *states, &lineages,
                matrices->transmat_compress, matrices->transmat);
        } else {
            matrices->transmat = NULL;
        }
        
        // non-compressed version
        //calc_transition_probs(tree, model, *states, &lineages, transmat);
    }


    

protected:
    ArgModel *model;
    Sequences *seqs;
    LocalTrees *trees;
    int new_chrom;
    bool calc_transmat;
    
    LocalTrees::iterator tree_iter;
    LineageCounts lineages;
    States states1;
    States states2;
    States *states;
    States *last_states;
    LocalTree *last_tree;
};


class ArgHmmMatrixList : public ArgHmmMatrixIter
{
public:
    ArgHmmMatrixList(ArgModel *model, Sequences *seqs, LocalTrees *trees, 
                     int new_chrom=-1, bool calc_transmat=true) :
        ArgHmmMatrixIter(model, seqs, trees, 
                         new_chrom, calc_transmat)
    {}
    ~ArgHmmMatrixList() 
    {
        clear();
    }

    // precompute all matrices
    void setup() 
    {
        setup(trees->begin(), trees->end());
    }

    // precompute all matrices
    void setup(LocalTrees::iterator start, LocalTrees::iterator end)
    {
        for (begin(); more(); next()) {
            matrices.push_back(ArgHmmMatrices());
            calc_matrices(&matrices.back());
        }
    }
    
    // free all computed matrices
    void clear()
    {
        for (unsigned int i=0; i<matrices.size(); i++)
            matrices[i].clear();
        matrices.clear();
    }


    virtual void begin()
    {
        ArgHmmMatrixIter::begin();
    }
    
    // initializes iterator
    virtual void begin(LocalTrees::iterator start)
    {
        ArgHmmMatrixIter::begin(start);
        matrix_index = 0;
    }

    virtual void rbegin()
    {
        ArgHmmMatrixIter::rbegin();
        matrix_index = matrices.size() - 1;
    }
    
    // initializes iterator
    virtual void rbegin(LocalTrees::iterator start)
    {
        ArgHmmMatrixIter::rbegin(start);
        matrix_index = matrices.size() - 1;
    }


    // moves iterator to next block
    virtual bool next()
    {
        matrix_index++;
        return ArgHmmMatrixIter::next();
    }

    // moves iterator to previous block
    virtual bool prev()
    {
        matrix_index--;
        return ArgHmmMatrixIter::prev();
    }


    virtual void get_matrices(ArgHmmMatrices *mat)
    {
        *mat = matrices[matrix_index];
    }

    vector<ArgHmmMatrices> matrices;

protected:
    int matrix_index;
};


//=============================================================================
// HMM algorithms


// compute one block of forward algorithm with compressed transition matrices
void arghmm_forward_alg_block_fast(LocalTree *tree, ArgModel *model, 
                                   int blocklen, const States &states, 
                                   LineageCounts &lineages,
                                   TransMatrixCompress *matrix,
                                   double **emit, double **fw)
{
    const int nstates = states.size();
    const int ntimes = model->ntimes;
    LocalNode *nodes = tree->nodes;

    // get aliases for various precomputed terms
    const double *B = matrix->B;
    const double *D = matrix->D;
    const double *E = matrix->E;
    const double *G = matrix->G;
    const double *norecombs = matrix->norecombs;
    const double *sums = matrix->sums;

    // compute ntimes*ntimes and ntime*nstates temp matrices
    double lnorecombs[ntimes];
    double tmatrix[ntimes][ntimes];
    double tmatrix2[ntimes][nstates];
    for (int a=0; a<ntimes-1; a++) {
        lnorecombs[a] = log(norecombs[a]);

        for (int b=0; b<ntimes-1; b++) {
            const double I = float(a <= b);
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

    // get normalization factors specific to the local tree
    double lsums[nstates];
    for (int j=0; j<nstates; j++)
        lsums[j] = log(sums[j]);


    // get max time
    int maxtime = 0;
    for (int k=0; k<nstates; k++)
        if (maxtime < states[k].time)
            maxtime = states[k].time;

    // get branch ages
    double ages1[tree->nnodes];
    double ages2[tree->nnodes];
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
        int a = states[j].time;
        state_map[j] = offset2[a];
        offset2[a]++;
    }


    NodeStateLookup state_lookup(states, tree->nnodes);

    double vec[nstates];
    double tmatrix_fgroups[ntimes];
    double fgroups[ntimes];
    for (int i=1; i<blocklen; i++) {
        double *col1 = fw[i-1];
        double *col2 = fw[i];
        double *emit2 = emit[i];
        
        // precompute the fgroup sums
        for (int j=0; j<nstates; j++)
            vec[state_map[j]] = col1[j] - lsums[j];
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
                vec[m++] = tmatrix2[a][k] + col1[j] - lsums[j];
            
            // same state case (add possibility of no recomb)
            vec[m++] = lnorecombs[b] + col1[k] - lsums[k];
            
            col2[k] = logsum(vec, m) + emit2[k];
        }
    }
}



// run forward algorithm with matrices precomputed
double **arghmm_forward_alg(LocalTrees *trees, ArgModel *model,
    Sequences *sequences, ArgHmmMatrixList *matrix_list, double **fw=NULL)
{
    LineageCounts lineages(model->ntimes);
    States states;

    // allocate the forward table if necessary
    if (fw == NULL)
        fw = new double* [sequences->seqlen];    

    // forward algorithm over local trees
    int blocki = 0;
    for (LocalTrees::iterator it=trees->begin(); 
         it != trees->end(); ++it, blocki++) {
        int pos = it->block.start;
        ArgHmmMatrices matrices = matrix_list->matrices[blocki];

        // allocate the forward table column if necessary
        for (int i=pos; i<pos+matrices.blocklen; i++)
            fw[i] = new double [matrices.nstates2];
        
        // use switch matrix for first column of forward table
        // if we have a previous state space (i.e. not first block)
        if (pos == 0) {
            // calculate prior of first state
            LocalTree *tree = it->tree;
            get_coal_states(tree, model->ntimes, states);
            lineages.count(tree);
            calc_state_priors(states, &lineages, model, fw[0]);
        } else {
            // perform one column of forward algorithm with transmat_switch
            forward_step(fw[pos-1], fw[pos], matrices.nstates1, 
                matrices.nstates2, matrices.transmat_switch, matrices.emit[0]);
        }

        // calculate rest of block
        forward_alg(matrices.blocklen, matrices.nstates2, 
                    matrices.transmat, matrices.emit, &fw[pos]);
    }
    
    return fw;
}


// run forward algorithm with matrices precomputed
double **arghmm_forward_alg_fast(LocalTrees *trees, ArgModel *model,
    Sequences *sequences, ArgHmmMatrixIter *matrix_list, double **fw=NULL)
{
    LineageCounts lineages(model->ntimes);
    States states;
    ArgHmmMatrices matrices;
    
    // allocate the forward table if necessary
    if (fw == NULL)
        fw = new double* [sequences->seqlen];

    // forward algorithm over local trees
    for (matrix_list->begin(); matrix_list->more(); matrix_list->next()) {
        LocalTrees::iterator it = matrix_list->get_tree_iter();
        int pos = it->block.start;
        matrix_list->get_matrices(&matrices);

        // allocate the forward table column if necessary
        for (int i=pos; i<pos+matrices.blocklen; i++)
            fw[i] = new double [matrices.nstates2];

        
        get_coal_states(it->tree, model->ntimes, states);
        lineages.count(it->tree);
        
        // use switch matrix for first column of forward table
        // if we have a previous state space (i.e. not first block)
        if (pos == 0) {
            // calculate prior of first state
            calc_state_priors(states, &lineages, model, fw[0]);
        } else {
            // perform one column of forward algorithm with transmat_switch
            forward_step(fw[pos-1], fw[pos], matrices.nstates1, 
                matrices.nstates2, matrices.transmat_switch, matrices.emit[0]);
        }

        // calculate rest of block
        arghmm_forward_alg_block_fast(it->tree, model, matrices.blocklen, 
                                      states, lineages, 
                                      matrices.transmat_compress,
                                      matrices.emit, &fw[pos]);
    }
    
    return fw;
}



void sample_hmm_posterior_fast(int n, LocalTree *tree, const States &states,
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
                trans[j] = log(matrix->get_transition_prob(tree, states, j, k));
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



void stochastic_traceback(ArgHmmMatrixList *matrix_list, 
                          double **fw, int *path, int seqlen)
{
    const int ntrees = matrix_list->matrices.size();

    // choose last column first
    int nstates = matrix_list->matrices[ntrees-1].nstates2;
    double total = logsum(fw[seqlen - 1], nstates);
    double vec[nstates];
    for (int i=0; i<nstates; i++)
        vec[i] = exp(fw[seqlen - 1][i] - total);
    path[seqlen - 1] = sample(vec, nstates);
    

    // iterate backward through blocks
    int pos = seqlen;
    for (int blocki = ntrees - 1; blocki >= 0; blocki--) {
        ArgHmmMatrices mat = matrix_list->matrices[blocki];
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


void stochastic_traceback_fast(LocalTrees *trees, ArgModel *model, 
                               ArgHmmMatrixIter *matrix_list, 
                               double **fw, int *path, int seqlen)
{
    ArgHmmMatrices mat;
    States states;

    // choose last column first
    matrix_list->rbegin();
    matrix_list->get_matrices(&mat);
    int nstates = mat.nstates2;
    double total = logsum(fw[seqlen - 1], nstates);
    double vec[nstates];
    for (int i=0; i<nstates; i++)
        vec[i] = exp(fw[seqlen - 1][i] - total);
    path[seqlen - 1] = sample(vec, nstates);
    
    // iterate backward through blocks
    int pos = seqlen;
    for (; matrix_list->more(); matrix_list->prev())
    {        
        matrix_list->get_matrices(&mat);
        LocalTree *tree = matrix_list->get_tree_iter()->tree;
        get_coal_states(tree, model->ntimes, states);
        pos -= mat.blocklen;
        
        sample_hmm_posterior_fast(mat.blocklen, tree, states,
                                  mat.transmat_compress, mat.emit, 
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



void sample_recombinations(
    LocalTrees *trees, ArgModel *model, ArgHmmMatrixList *matrix_list,
    int *thread_path, vector<int> &recomb_pos, vector<NodePoint> &recombs)
{
    States states;
    LineageCounts lineages(model->ntimes);
    const int new_node = -1;
    vector <NodePoint> candidates;
    vector <double> probs;


    // loop through local blocks
    int blocki = 0;
    for (LocalTrees::iterator it=trees->begin(); 
         it != trees->end(); ++it, blocki++) {

        // get local block information
        int start = it->block.start + 1;  // don't allow new recomb at start
        int end = it->block.end;
        ArgHmmMatrices matrices = matrix_list->matrices[blocki];
        LocalTree *tree = it->tree;
        double treelen_b = get_treelen(tree, model->times, model->ntimes);
        double treelen = treelen_b - get_basal_branch(
            tree, model->times, model->ntimes, -1, -1);
        double treelen2_b, treelen2;
        lineages.count(tree);
        get_coal_states(tree, model->ntimes, states);
        int statei = thread_path[start];
        int next_recomb = -1;

        
        // loop through positions in block
        for (int i=start; i<end; i++) {
            
            if (thread_path[i] == thread_path[i-1]) {
                // no change in state, recombination is optional
                if (i > next_recomb) {
                    // sample the next recomb pos
                    int last_state = thread_path[i-1];
                    treelen2_b = get_treelen_branch(
                        tree, model->times, model->ntimes,
                        states[last_state].node,
                        states[last_state].time, treelen_b);
                    treelen2 = treelen2_b - get_basal_branch(
                        tree, model->times, model->ntimes, 
                        states[last_state].node, states[last_state].time);

                    double self_trans;
                    if (matrices.transmat)
                        self_trans = matrices.transmat[last_state][last_state];
                    else
                        self_trans = log(matrices.transmat_compress->get_transition_prob(tree, states, last_state, last_state));
                    double rate = max(
                        1.0 - exp(-model->rho * (treelen2 - treelen)
                              - self_trans), model->rho);

                    // NOTE: the min prevents large floats from overflowing
                    // when cast to int
                    next_recomb = int(min(float(end), i + expovariate(rate)));
                }

                if (i < next_recomb)
                    continue;
            }

            // sample recombination
            next_recomb = -1;
            statei = thread_path[i];
            State state = states[statei];
            State last_state = states[thread_path[i-1]];
            treelen2_b = get_treelen_branch(tree, model->times, model->ntimes,
                                            last_state.node,
                                            last_state.time, treelen_b);
            treelen2 = treelen2_b - get_basal_branch(
                tree, model->times, model->ntimes, 
                last_state.node, last_state.time);
            

            // there must be a recombination
            // either because state changed or we choose to recombine
            // find candidates
            candidates.clear();
            int end_time = min(state.time, last_state.time);
            if (state.node == last_state.node) {
                // y = v, k in [0, min(timei, last_timei)]
                // y = node, k in Sr(node)
                for (int k=tree->nodes[state.node].age; k<=end_time; k++)
                    candidates.push_back(NodePoint(state.node, k));
            }
            
            for (int k=0; k<=end_time; k++)
                candidates.push_back(NodePoint(new_node, k));


            // compute probability of each candidate
            double C[model->ntimes];
            C[0] = 0.0;
            for (int b=1; b<model->ntimes; b++) {
                const int l = b - 1;
                C[b] = C[l] + model->time_steps[l] * lineages.nbranches[l] / 
                    (2.0 * model->popsizes[l]);
            }

            probs.clear();
            int j = state.time;
            for (vector<NodePoint>::iterator it=candidates.begin(); 
                 it != candidates.end(); ++it) {
                int k = it->time;
                double p = 
                    (lineages.nbranches[k] + 1) * model->time_steps[k] /
                    (lineages.ncoals[j] * (lineages.nrecombs[k] + 1.0) * 
                     treelen2_b) * 
                    (1.0 - exp(-model->time_steps[j] * lineages.nbranches[j] /
                               (2.0 * model->popsizes[j-1]))) *
                    (1.0 - exp(-model->rho * treelen2)) *
                    exp(-C[k] + C[k]);
                probs.push_back(p);
            }

            // sample recombination
            recomb_pos.push_back(i);
            recombs.push_back(candidates[sample(&probs[0], probs.size())]);

            assert(recombs[recombs.size()-1].time <= min(state.time,
                                                         last_state.time));
        }
    }
}


//=============================================================================
// ARG sampling


// sample the thread of the last chromosome
void sample_arg_thread(ArgModel *model, Sequences *sequences, 
                       LocalTrees *trees, int new_chrom,
                       double **fw=NULL, int *thread_path=NULL)
{
    bool tmp_given = (fw != NULL);
    //Timer time;

    // allocate temp variables
    if (fw == NULL) {
        fw = new double* [sequences->seqlen];
        thread_path = new int [sequences->seqlen];
    }
    assert(thread_path != NULL);    

    // build matrices
    Timer time;
    ArgHmmMatrixList matrix_list(model, sequences, trees, new_chrom, false);
    matrix_list.setup();
    printf("matrix calc: %e s\n", time.time());
    
    // compute forward table
    time.start();
    //arghmm_forward_alg(trees, model, sequences, &matrix_list, fw);
    arghmm_forward_alg_fast(trees, model, sequences, &matrix_list, fw);
    printf("forward:     %e s  (%d states, %d blocks)\n", time.time(),
           matrix_list.matrices[0].nstates2,
           (int) matrix_list.matrices.size());


    // traceback
    time.start();
    //stochastic_traceback(&matrix_list, fw, thread_path, sequences->seqlen);
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

    // cleanup
    for (int i=0; i<sequences->seqlen; i++)
        delete [] fw[i];

    // clean up
    if (!tmp_given) {
        delete [] fw;
        delete [] thread_path;
    }
}


// sequentially sample an ARG from scratch
// sequences are sampled in the order given
void sample_arg_seq(ArgModel *model, Sequences *sequences, LocalTrees *trees)
{
    const int seqlen = sequences->seqlen;

    // initialize ARG as trunk
    const int capacity = 2 * sequences->nseqs - 1;
    trees->make_trunk(0, seqlen, capacity);

    // allocate temp variables
    double **fw = new double* [seqlen];
    int *thread_path = new int [seqlen];

    // add more chromosomes one by one
    for (int nchroms=2; nchroms<=sequences->nseqs; nchroms++) {
        // use first nchroms sequences
        Sequences sequences2(sequences->seqs, nchroms, seqlen);
        int new_chrom = nchroms - 1;
        sample_arg_thread(model, &sequences2, trees, new_chrom,fw,thread_path);
    }

    // clean up    
    delete [] fw;
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


// resample the threading of all the chromosomes
void resample_arg(ArgModel *model, Sequences *sequences, LocalTrees *trees)
{
    const int nleaves = trees->get_num_leaves();
    for (int chrom=0; chrom<nleaves; chrom++) {
        // remove chromosome from ARG and resample its thread
        remove_arg_thread(trees, chrom);
        sample_arg_thread(model, sequences, trees, chrom);
    }
}



//=============================================================================
// C interface
extern "C" {

double **arghmm_forward_alg(
    int **ptrees, int **ages, int **sprs, int *blocklens,
    int ntrees, int nnodes, double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, double **fw=NULL)
{    
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    LocalTrees trees(ptrees, ages, sprs, blocklens, ntrees, nnodes);
    Sequences sequences(seqs, nseqs, seqlen);

    // build matrices
    ArgHmmMatrixList matrix_list(&model, &sequences, &trees);
    matrix_list.setup();

    return arghmm_forward_alg_fast(&trees, &model, &sequences, &matrix_list);
    //return arghmm_forward_alg(&trees, &model, &sequences, fw);
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
    double **fw = new double* [seqlen];
    arghmm_forward_alg(&trees, &model, &sequences, &matrix_list, fw);

    // traceback
    int *ipath = new int [seqlen];
    stochastic_traceback(&matrix_list, fw, ipath, seqlen);

    
    // convert path
    if (path == NULL)
        path = new intstate [seqlen];

    States states;
    for (LocalTrees::iterator it=trees.begin(); it != trees.end(); ++it) {
        int start = it->block.start;
        int end = it->block.end;
        get_coal_states(it->tree, ntimes, states);

        for (int i=start; i<end; i++) {
            int istate = ipath[i];
            path[i][0] = states[istate].node;
            path[i][1] = states[istate].time;
        }
    }


    // clean up
    delete_matrix<double>(fw, seqlen);
    delete [] ipath;

    return path;
}


LocalTrees *arghmm_sample_thread(
    int **ptrees, int **ages, int **sprs, int *blocklens,
    int ntrees, int nnodes, double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    LocalTrees *trees = new LocalTrees(ptrees, ages, sprs, blocklens, 
                                       ntrees, nnodes);
    Sequences sequences(seqs, nseqs, seqlen);
    
    // build matrices
    ArgHmmMatrixList matrix_list(&model, &sequences, trees);
    matrix_list.setup();
    
    // compute forward table
    double **fw = new double* [seqlen];
    arghmm_forward_alg(trees, &model, &sequences, &matrix_list, fw);

    // traceback
    int *thread_path = new int [seqlen];
    stochastic_traceback(&matrix_list, fw, thread_path, seqlen);

    // sample recombination points
    vector<int> recomb_pos;
    vector<NodePoint> recombs;
    sample_recombinations(trees, &model, &matrix_list,
                          thread_path, recomb_pos, recombs);

    // add thread to ARG
    //assert_trees_thread(trees, thread_path, ntimes);
    int new_chrom = trees->get_num_leaves();
    add_arg_thread(trees, model.ntimes, thread_path, new_chrom, 
                   recomb_pos, recombs);

    // clean up
    delete_matrix<double>(fw, seqlen);
    delete [] thread_path;

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
    char **seqs, int nseqs, int seqlen, int niters)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    LocalTrees *trees = new LocalTrees();    

    sample_arg_seq(&model, &sequences, trees);
    for (int i=0; i<niters; i++)
        resample_arg(&model, &sequences, trees);
    
    return trees;
}


// resample an ARG with gibbs
LocalTrees *arghmm_resample_arg(
    int **ptrees, int **ages, int **sprs, int *blocklens,
    int ntrees, int nnodes, 
    double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, int niters)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    LocalTrees *trees = new LocalTrees(ptrees, ages, sprs, blocklens, 
                                       ntrees, nnodes);

    // sequentially sample until all chromosomes are present
    // then gibbs
    for (int new_chrom=trees->get_num_leaves(); new_chrom<nseqs; new_chrom++) {
        sample_arg_thread(&model, &sequences, trees, new_chrom);
    }

    for (int i=0; i<niters; i++)
        resample_arg(&model, &sequences, trees);
    
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

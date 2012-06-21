//=============================================================================
// transition and emission matrices

#ifndef ARGHMM_MATRICES_H
#define ARGHMM_MATRICES_H


// c++ includes
#include <list>
#include <vector>
#include <string.h>

// arghmm includes
#include "common.h"
#include "emit.h"
#include "local_tree.h"
#include "logging.h"
#include "model.h"
#include "seq.h"
#include "states.h"
#include "trans.h"



namespace arghmm {

using namespace std;
using namespace spidir;


class ArgHmmMatrices
{
public:
    ArgHmmMatrices() :
        nstates1(0),
        nstates2(0),
        blocklen(0),
        transmat_compress(NULL),
        transmat_switch_compress(NULL),
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
                   TransMatrixSwitchCompress *transmat_switch_compress,
                   double **transmat, double **transmat_switch, double **emit):
        nstates1(nstates1),
        nstates2(nstates2),
        blocklen(blocklen),
        transmat_compress(transmat_compress),
        transmat_switch_compress(transmat_switch_compress),
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
        if (transmat_switch_compress) {
            delete transmat_switch_compress;
            transmat_switch_compress = NULL;
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
    TransMatrixSwitchCompress* transmat_switch_compress;
    double **transmat;
    double **transmat_switch;
    double **emit;
};


// iterates through matricies for the ArgHmm
class ArgHmmMatrixIter
{
public:
    ArgHmmMatrixIter(ArgModel *model, Sequences *seqs, LocalTrees *trees, 
                     int _new_chrom=-1, bool calc_full=true) :
        model(model),
        seqs(seqs),
        trees(trees),
        new_chrom(_new_chrom),
        calc_full(calc_full),
        pos(trees->start_coord),
        lineages(model->ntimes)        
    {
        if (new_chrom == -1)
            new_chrom = trees->get_num_leaves();
    }

    virtual ~ArgHmmMatrixIter() {}

    // initializes iterator
    virtual void begin()
    {
        begin(trees->begin(), trees->start_coord);
    }
    
    virtual void begin(LocalTrees::iterator start, int start_coord)
    {
        tree_iter = start;
        pos = start_coord;
        
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
        rbegin(it, trees->end_coord - it->blocklen);
    }
    
    virtual void rbegin(LocalTrees::iterator start, int start_coord)
    {
        begin(start, start_coord);
        last_states = &states2;
    }


    virtual bool next()
    {
        // update pointers
        last_tree = tree_iter->tree;
        last_states = states;
        states = ((states == &states1) ? &states2 : &states1);
        pos += tree_iter->blocklen;

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
            pos -= it->blocklen;
        } else {
            last_tree = NULL;
            pos = trees->start_coord;
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

    int get_position()
    {
        return pos;
    }

    
    // calculate transition and emission matrices for current block
    void calc_matrices(ArgHmmMatrices *matrices)
    {
        int blocklen = tree_iter->blocklen;
        LocalTree *tree = tree_iter->tree;
        get_coal_states(tree, model->ntimes, *states);
        int nstates = states->size();
        matrices->blocklen = blocklen;
        int nleaves = trees->get_num_leaves();
        
        // calculate emissions
        if (seqs) {
            char *subseqs[seqs->nseqs];
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
            matrices->nstates1 = matrices->nstates2 = nstates;
            
        } else {
            matrices->nstates1 = last_states->size();
            matrices->nstates2 = nstates;
            lineages.count(last_tree);
                
            // calculate transmat_switch
            matrices->transmat_switch_compress = new TransMatrixSwitchCompress(
                matrices->nstates1, matrices->nstates2);

            calc_transition_probs_switch_compress(tree, last_tree, 
                tree_iter->spr, tree_iter->mapping,
                *last_states, *states, model, &lineages, 
                matrices->transmat_switch_compress);

            if (calc_full) {
                matrices->transmat_switch = new_matrix<double>(
                    matrices->nstates1, matrices->nstates2);
                calc_transition_switch_probs(matrices->transmat_switch_compress,
                                             matrices->transmat_switch);
            } else {
                matrices->transmat_switch = NULL;
            }
            //calc_transition_probs_switch(
            //    tree, last_tree, tree_iter->spr, tree_iter->mapping,
            //    *last_states, *states,
            //    model, &lineages, matrices->transmat_switch);
        }

        // update lineages to current tree
        lineages.count(tree);
        
        // calculate transmat and use it for rest of block
        matrices->transmat_compress = new TransMatrixCompress(
            model->ntimes, nstates);
        calc_transition_probs_compress(tree, model, *states, &lineages, 
                                       matrices->transmat_compress);
        if (calc_full) {
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
    bool calc_full;
    
    int pos;
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
                     int new_chrom=-1, bool calc_full=true) :
        ArgHmmMatrixIter(model, seqs, trees, 
                         new_chrom, calc_full)
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
        matrix_index = 0;
    }
    
    // initializes iterator
    virtual void begin(LocalTrees::iterator start, int start_coord)
    {
        ArgHmmMatrixIter::begin(start, start_coord);
        matrix_index = 0;
    }

    virtual void rbegin()
    {
        ArgHmmMatrixIter::rbegin();
        matrix_index = matrices.size() - 1;
    }
    
    // initializes iterator
    virtual void rbegin(LocalTrees::iterator start, int start_coord)
    {
        ArgHmmMatrixIter::rbegin(start, start_coord);
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


} // namespace arghmm


#endif // ARGHMM_MATRICES_H

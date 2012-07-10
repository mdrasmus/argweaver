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
#include "sequences.h"
#include "states.h"
#include "trans.h"



namespace arghmm {

using namespace std;

class ArgHmmMatrices
{
public:
    ArgHmmMatrices() :
        nstates1(0),
        nstates2(0),
        blocklen(0),
        transmat(NULL),
        transmat_switch(NULL),
        transprobs(NULL),
        transprobs_switch(NULL),
        emit(NULL)
    {}

    ArgHmmMatrices(int nstates1, int nstates2, int blocklen,
                   double **transprobs, double **transprobs_switch, 
                   double **emit):
        nstates1(nstates1),
        nstates2(nstates2),
        blocklen(blocklen),
        transmat(NULL),
        transmat_switch(NULL),
        transprobs(transprobs),
        transprobs_switch(transprobs_switch),
        emit(emit)
    {}

    ArgHmmMatrices(int nstates1, int nstates2, int blocklen,
                   TransMatrix *transmat,
                   TransMatrixSwitch *transmat_switch,
                   double **transprobs, double **transprobs_switch, 
                   double **emit):
        nstates1(nstates1),
        nstates2(nstates2),
        blocklen(blocklen),
        transmat(transmat),
        transmat_switch(transmat_switch),
        transprobs(transprobs),
        transprobs_switch(transprobs_switch),
        emit(emit)
    {}

    ~ArgHmmMatrices() 
    {}
    
    // delete all matrices
    void clear()
    {
        if (transmat) {
            delete transmat;
            transmat = NULL;
        }
        if (transmat_switch) {
            delete transmat_switch;
            transmat_switch = NULL;
        }
        if (transprobs) {
            delete_matrix<double>(transprobs, nstates2);
            transprobs = NULL;
        }
        if (transprobs_switch) {
            delete_matrix<double>(transprobs_switch, nstates1);
            transprobs_switch = NULL;
        }
        if (emit) {
            delete_matrix<double>(emit, blocklen);
            emit = NULL;
        }
    }


    // release ownership of underlying data
    void detach() 
    {
        transmat = NULL;
        transmat_switch = NULL;
        transprobs = NULL;
        transprobs_switch = NULL;
        emit = NULL;
    }

    
    int nstates1;
    int nstates2;
    int blocklen;
    TransMatrix* transmat;
    TransMatrixSwitch* transmat_switch;
    double **transprobs;
    double **transprobs_switch;
    double **emit;
};


// iterates through matricies for the ArgHmm
class ArgHmmMatrixIter
{
public:
    ArgHmmMatrixIter(ArgModel *model, Sequences *seqs, LocalTrees *trees, 
                     int _new_chrom=-1, bool calc_full=false) :
        model(model),
        seqs(seqs),
        trees(trees),
        new_chrom(_new_chrom),
        calc_full(calc_full),
        internal(false),
        pos(trees->start_coord),
        lineages(model->ntimes)
    {
        if (new_chrom == -1)
            new_chrom = trees->get_num_leaves();
    }

    virtual ~ArgHmmMatrixIter() 
    {
        mat.clear();
    }

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
        if (start_coord == trees->start_coord) {
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
            if (internal)
                get_coal_states_internal(
                    last_tree, model->ntimes, *last_states);
            else
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
            last_states = &states2;
            if (internal)
                get_coal_states_internal(
                    last_tree, model->ntimes, *last_states);
            else
                get_coal_states(last_tree, model->ntimes, *last_states);
            pos -= it->blocklen;
        } else {
            last_tree = NULL;
            last_states = NULL;
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
        mat.clear();
        calc_matrices(&mat);
        *matrices = mat;
    }

    LocalTrees::iterator get_tree_iter()
    {
        return tree_iter;
    }

    int get_position()
    {
        return pos;
    }

    void set_internal(bool _internal)
    {
        internal = _internal;
    }

    
    // calculate transition and emission matrices for current block
    void calc_matrices(ArgHmmMatrices *matrices)
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
            char *subseqs[seqs->get_nseqs()];
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
            // TODO: need to make a internal branch version
            matrices->transprobs = new_matrix<double>(nstates, nstates);
            get_transition_probs(tree, model, *states, &lineages,
                matrices->transmat, matrices->transprobs);
        } else {
            matrices->transprobs = NULL;
        }
    }



    // calculate transition and emission matrices for current block
    void calc_matrices_internal(ArgHmmMatrices *matrices)
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
            matrices->emit = new_matrix<double>(blocklen, nstates);
            calc_emissions_internal(*states, tree, subseqs, nleaves, 
                                    blocklen, model, matrices->emit);
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

            calc_transition_probs_switch_internal(tree, last_tree, 
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
        calc_transition_probs_internal(tree, model, *states, &lineages, 
                                       matrices->transmat);
        
        if (calc_full) {
            // TODO: need to make a internal branch version
            matrices->transprobs = new_matrix<double>(nstates, nstates);
            get_transition_probs(tree, model, *states, &lineages,
                matrices->transmat, matrices->transprobs);
        } else {
            matrices->transprobs = NULL;
        }
    }
    
    

protected:
    ArgModel *model;
    Sequences *seqs;
    LocalTrees *trees;
    int new_chrom;
    bool calc_full;
    bool internal;

    ArgHmmMatrices mat;
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
                     int new_chrom=-1, bool calc_full=false) :
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

    
protected:
    int matrix_index;
    vector<ArgHmmMatrices> matrices;
};


} // namespace arghmm


#endif // ARGHMM_MATRICES_H

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


// Transition and emission matrices for one non-recombining block in the ArgHmm
class ArgHmmMatrices
{
public:
    ArgHmmMatrices() :
        nstates1(0),
        nstates2(0),
        blocklen(0),
        transmat(NULL),
        transmat_switch(NULL),
        emit(NULL)
    {}

    ArgHmmMatrices(int nstates1, int nstates2, int blocklen,
                   TransMatrix *transmat,
                   TransMatrixSwitch *transmat_switch,
                   double **emit):
        nstates1(nstates1),
        nstates2(nstates2),
        blocklen(blocklen),
        transmat(transmat),
        transmat_switch(transmat_switch),
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
        emit = NULL;
    }

    
    int nstates1; // number of states in previous block
    int nstates2; // number of states in this block
    int blocklen; // block length
    TransMatrix* transmat; // transition matrix within this block
    TransMatrixSwitch* transmat_switch; // transition matrix from previous block
    double **emit; // emission matrix
};


void calc_matrices(
    const ArgModel *model, const Sequences *seqs, const LocalTrees *trees,
    const LocalTreeSpr *last_tree_spr, const LocalTreeSpr *tree_spr,
    const int start, const int end, const int new_chrom,
    const bool internal, ArgHmmMatrices *matrices);



// A block of the ARG and model
class ArgModelBlock
{
public:
    ArgModelBlock() {}

    ArgModelBlock(int start, int end, const LocalTreeSpr *tree_spr, 
                  int model_index) :
        start(start),
        end(end),
        tree_spr(tree_spr),
        model_index(model_index)
    {}

    int length() const {
        return end - start;
    }

    int start;
    int end;
    const LocalTreeSpr* tree_spr;
    int model_index;
};


// The common blocks of the ARG and model
class ArgModelBlocks
{
public:
    ArgModelBlocks(const ArgModel *model, const LocalTrees *trees) :
        model(model),
        trees(trees)
    {
        assert(model->has_recombmap() && model->has_mutmap());
    }
        
    void setup() {
        int start = trees->start_coord;
        int end = start;
        int tree_start = start;
        int model_index = model->mutmap.index(start);
        LocalTrees::const_iterator tree_iter = trees->begin();
        int tree_end = tree_start + tree_iter->blocklen;
        int model_end = model->mutmap[model_index].end; 

        // record all common blocks
        while (start < trees->end_coord) {
            start = end;
            
            if (tree_end <= model_end) {
                // tree block end is next
                tree_start += tree_iter->blocklen;
                ++tree_iter;
            } 

            if (tree_end >= model_end) {
                // model block end is next
                model_index++;
            }
            
            // check if no more blocks
            if (tree_iter == trees->end())
                break;
        
            // compute new block ends
            tree_end = tree_start + tree_iter->blocklen;
            model_end = model->mutmap[model_index].end; 
            
            // determine next block end
            end = min(tree_end, model_end);
            
            // record block
            blocks.push_back(ArgModelBlock(start, end, 
                                           &(*tree_iter), model_index));
        }
    }

    ArgModelBlock &at(int i) {
        return blocks[i];
    }

    const ArgModelBlock &at(int i) const {
        return blocks[i];
    }

    int size() const {
        return blocks.size();
    }

protected:
    const ArgModel *model;
    const LocalTrees *trees;
    
    vector<ArgModelBlock> blocks;
};



// iterates through matricies for the ArgHmm
class ArgHmmMatrixIter2
{
public:
    ArgHmmMatrixIter2(const ArgModel *model, const Sequences *seqs, 
                      const LocalTrees *trees, 
                      int _new_chrom=-1) :
        model(model),
        seqs(seqs),
        trees(trees),
        new_chrom(_new_chrom),
        internal(false),
        blocks(model, trees)
    {
        if (new_chrom == -1)
            new_chrom = trees->get_num_leaves();
    }

    virtual ~ArgHmmMatrixIter2() 
    {
        mat.clear();
    }

    void setup() {
        // determine all blocks
        blocks.setup();
    }

    // initializes iterator
    virtual void begin()
    {
        block_index = 0;
    }
    
    virtual void rbegin()
    {
        block_index = blocks.size() - 1;
    }

    virtual bool next()
    {
        block_index++;
        return block_index < int(blocks.size());
    }


    // moves iterator to previous block
    virtual bool prev()
    {
        block_index--;
        return block_index >= 0;
    }
    

    // returns true if there are more blocks
    virtual bool more() const {
        return block_index >= 0 && block_index < blocks.size();
    }


    virtual void get_matrices(ArgHmmMatrices *matrices)
    {
        mat.clear();
        calc_matrices(&mat);
        *matrices = mat;
    }

    void calc_matrices(ArgHmmMatrices *matrices)
    {
        ArgModel local_model;
        ArgModelBlock &block = blocks.at(block_index);
        
        model->get_local_model(block.start, local_model);
        LocalTreeSpr const* last_tree_spr = NULL;
        if (block_index > 0)
            last_tree_spr = blocks.at(block_index-1).tree_spr;
        
        arghmm::calc_matrices(
            &local_model, seqs, trees, last_tree_spr, block.tree_spr,
            block.start, block.end, new_chrom, internal, matrices);
    }


    const LocalTreeSpr *get_tree_spr() const
    {
        return blocks.at(block_index).tree_spr;
    }

    int get_block_start() const {
        return blocks.at(block_index).start;
    }

    int get_block_end() const {
        return blocks.at(block_index).end;
    }

    int get_blocklen() const {
        return blocks.at(block_index).length();
    }

    int get_model_index() const {
        return blocks.at(block_index).model_index;
    }

    void get_local_model(ArgModel &local_model) {
        model->get_local_model(blocks.at(block_index).start, local_model);
    }

    void set_internal(bool _internal)
    {
        internal = _internal;
    }
    
protected:
    // references to model, arg, sequences
    const ArgModel *model;
    const Sequences *seqs;
    const LocalTrees *trees;
    int new_chrom;
    bool internal;
    
    ArgHmmMatrices mat;

    // record of common blocks
    ArgModelBlocks blocks;
    int block_index;
};



// iterates through matricies for the ArgHmm
class ArgHmmMatrixIter
{
public:
    ArgHmmMatrixIter(const ArgModel *model, const Sequences *seqs, 
                     const LocalTrees *trees, 
                     int _new_chrom=-1) :
        model(model),
        seqs(seqs),
        trees(trees),
        new_chrom(_new_chrom),
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
    
    virtual void begin(LocalTrees::const_iterator start, int start_coord)
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
            LocalTrees::const_iterator tree_iter2 = tree_iter;
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
        LocalTrees::const_iterator it = --trees->end();
        rbegin(it, trees->end_coord - it->blocklen);
    }
    
    virtual void rbegin(LocalTrees::const_iterator start, int start_coord)
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
        LocalTrees::const_iterator it = tree_iter;
        --it;

        if (it != trees->end()) {
            last_tree = it->tree;
            last_states = &states2;
            if (internal)
                get_coal_states_internal(
                    last_tree, model->ntimes, *last_states);
            else
                get_coal_states(last_tree, model->ntimes, *last_states);
            pos -= tree_iter->blocklen;
        } else {
            last_tree = NULL;
            last_states = NULL;
            pos = trees->start_coord;
        }
        
        return tree_iter != trees->end();
    }
    

    // returns true if there are more blocks
    virtual bool more() const
    {
        return tree_iter != trees->end();
    }


    virtual void get_matrices(ArgHmmMatrices *matrices)
    {
        mat.clear();
        calc_matrices(&mat);
        *matrices = mat;
    }

    virtual ArgHmmMatrices &ref_matrices()
    {
        mat.clear();
        calc_matrices(&mat);
        return mat;
    }


    LocalTrees::const_iterator get_tree_iter() const
    {
        return tree_iter;
    }

    int get_position() const
    {
        return pos;
    }

    void set_internal(bool _internal)
    {
        internal = _internal;
    }

    
    // calculate transition and emission matrices for current block
    void calc_matrices(ArgHmmMatrices *matrices);

    // calculate transition and emission matrices for current block
    void calc_matrices_internal(ArgHmmMatrices *matrices);

protected:
    // references to model, arg, sequences
    const ArgModel *model;
    ArgModel local_model;
    const Sequences *seqs;
    const LocalTrees *trees;
    int new_chrom;
    bool internal;
    
    // data for current block
    ArgHmmMatrices mat;
    int pos;
    LocalTrees::const_iterator tree_iter;
    States states1;
    States states2;
    States *states;
    States *last_states;
    LocalTree *last_tree;
    LineageCounts lineages;
};



// compute all emission and transition matrices and store them in a list
class ArgHmmMatrixList : public ArgHmmMatrixIter
{
public:
    ArgHmmMatrixList(const ArgModel *model, const Sequences *seqs, 
                     const LocalTrees *trees, 
                     int new_chrom=-1) :
        ArgHmmMatrixIter(model, seqs, trees, 
                         new_chrom)
    {}
    virtual ~ArgHmmMatrixList() 
    {
        clear();
    }

    // precompute all matrices
    void setup() 
    {
        setup(trees->begin(), trees->end());
    }

    // precompute all matrices
    void setup(LocalTrees::const_iterator start, LocalTrees::const_iterator end)
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

    virtual ArgHmmMatrices &ref_matrices()
    {
        return matrices[matrix_index];
    }

    
protected:
    int matrix_index;
    vector<ArgHmmMatrices> matrices;
};


} // namespace arghmm


#endif // ARGHMM_MATRICES_H

//=============================================================================
// transition and emission matrices

#ifndef ARGWEAVER_MATRICES_H
#define ARGWEAVER_MATRICES_H


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



namespace argweaver {

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

    void set_states(int ntimes, bool internal, int minage=0)
    {
        states_model.set(ntimes, internal, minage);
    }


    int nstates1; // number of states in previous block
    int nstates2; // number of states in this block
    int blocklen; // block length
    StatesModel states_model;
    TransMatrix* transmat; // transition matrix within this block
    TransMatrixSwitch* transmat_switch; // transition matrix from previous block
    double **emit; // emission matrix
};


void calc_arghmm_matrices(
    const ArgModel *model, const Sequences *seqs, const LocalTrees *trees,
    const LocalTreeSpr *last_tree_spr, const LocalTreeSpr *tree_spr,
    const int start, const int end, const int new_chrom,
    const StatesModel &states_model, ArgHmmMatrices *matrices);



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
    {}

    void setup() {
        int start = trees->start_coord;
        int tree_start = start;
        LocalTrees::const_iterator tree_iter = trees->begin();
        int tree_end = tree_start + tree_iter->blocklen;

        int model_index = -1;
        int model_end = trees->end_coord;
        if (model->has_mutmap()) {
            model_index = model->mutmap.index(start);
            model_end = model->mutmap[model_index].end;
        }

        // record all common blocks
        while (start < trees->end_coord) {
            // determine next block end
            int end = min(tree_end, model_end);

            // record block
            blocks.push_back(ArgModelBlock(start, end,
                                           &(*tree_iter), model_index));

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
            if (model->has_mutmap()) {
                assert(model_index < int(model->mutmap.size()));
                model_end = model->mutmap[model_index].end;
            }

            start = end;
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
class ArgHmmMatrixIter
{
public:
    ArgHmmMatrixIter(const ArgModel *model, const Sequences *seqs,
                      const LocalTrees *trees,
                      int _new_chrom=-1) :
        states_model(model->ntimes),
        model(model),
        seqs(seqs),
        trees(trees),
        new_chrom(_new_chrom),
        blocks(model, trees)
    {
        if (new_chrom == -1)
            new_chrom = trees->get_num_leaves();
    }

    virtual ~ArgHmmMatrixIter()
    {
        mat.clear();
    }

    virtual void setup() {
        // determine all blocks
        blocks.setup();
    }

    virtual void clear() {}

    // calculate matrix for internal branch resampling
    void set_internal(bool internal, int minage=0)
    {
        states_model.set(model->ntimes, internal, minage);
    }

    //==================================================
    // iteration methods

    // initializes iterator
    virtual void begin()
    {
        if (blocks.size() == 0)
            setup();
        block_index = 0;
    }

    virtual void rbegin()
    {
        if (blocks.size() == 0)
            setup();
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

    //==================================================
    // accessors

    virtual ArgHmmMatrices &ref_matrices()
    {
        mat.clear();
        calc_matrices(&mat);
        return mat;
    }

    const LocalTreeSpr *get_tree_spr() const {
        return blocks.at(block_index).tree_spr;
    }

    const LocalTreeSpr *get_last_tree_spr() const {
        if (block_index > 0)
            return blocks.at(block_index-1).tree_spr;
        else
            return NULL;
    }

    // Returns true if block starts with switch transition matrix
    bool has_switch() const {
        const LocalTreeSpr * last_tree_spr = get_last_tree_spr();
        return last_tree_spr && last_tree_spr != get_tree_spr();
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

    void get_coal_states(const LocalTree *tree, States &states) const {
        states_model.get_coal_states(tree, states);
    }

    void get_coal_states(States &states) const {
        states_model.get_coal_states(get_tree_spr()->tree, states);
    }


    StatesModel states_model;

protected:

    void calc_matrices(ArgHmmMatrices *matrices)
    {
        ArgModel local_model;
        ArgModelBlock &block = blocks.at(block_index);

        model->get_local_model_index(block.model_index, local_model);
        const LocalTreeSpr * last_tree_spr = get_last_tree_spr();

        argweaver::calc_arghmm_matrices(
            &local_model, seqs, trees, last_tree_spr, block.tree_spr,
            block.start, block.end, new_chrom, states_model, matrices);
    }


    // references to model, arg, sequences
    const ArgModel *model;
    const Sequences *seqs;
    const LocalTrees *trees;
    int new_chrom;

    ArgHmmMatrices mat;

    // record of common blocks
    ArgModelBlocks blocks;
    int block_index;
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
    virtual void setup()
    {
        ArgHmmMatrixIter::setup();

        for (begin(); more(); next()) {
            matrices.push_back(ArgHmmMatrices());
            calc_matrices(&matrices.back());
        }
    }

    // free all computed matrices
    virtual void clear()
    {
        for (unsigned int i=0; i<matrices.size(); i++)
            matrices[i].clear();
        matrices.clear();
    }

    //==================================================
    // iteration methods

    virtual void begin()
    {
        ArgHmmMatrixIter::begin();
        matrix_index = 0;
    }

    virtual void rbegin()
    {
        ArgHmmMatrixIter::rbegin();
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

    //==================================================
    // accessors

    virtual ArgHmmMatrices &ref_matrices()
    {
        return matrices[matrix_index];
    }


protected:
    int matrix_index;
    vector<ArgHmmMatrices> matrices;
};


} // namespace argweaver


#endif // ARGWEAVER_MATRICES_H

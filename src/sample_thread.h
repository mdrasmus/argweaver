#ifndef ARGWEAVER_SAMPLE_THREAD_H
#define ARGWEAVER_SAMPLE_THREAD_H

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
#include "sequences.h"
#include "sequences.h"
#include "states.h"
#include "thread.h"
#include "trans.h"



namespace argweaver {

using namespace std;


//=============================================================================
// Forward tables


class ArgHmmForwardTable
{
public:
    ArgHmmForwardTable(int start_coord, int seqlen) :
        start_coord(start_coord),
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
        nstates = max(nstates, 1);
        int blocklen = end - start;
        double *block = new double [blocklen * nstates];
        blocks.push_back(block);

        // link block to fw table
        for (int i=start; i<end; i++) {
            assert(i-start_coord >= 0 && i-start_coord < seqlen);
            fw[i-start_coord] = &block[(i-start)*nstates];
        }
    }

    // delete all blocks
    virtual void delete_blocks()
    {
        for (unsigned int i=0; i<blocks.size(); i++)
            delete [] blocks[i];
        blocks.clear();
    }

    virtual double **get_table()
    {
        return &fw[-start_coord];
    }

    virtual double **detach_table()
    {
        double **ptr = fw;
        fw = NULL;
        return ptr;
    }

    int start_coord;
    int seqlen;

protected:
    double **fw;
    vector<double*> blocks;
};


// older style allocation for testing with python
class ArgHmmForwardTableOld : public ArgHmmForwardTable
{
public:
    ArgHmmForwardTableOld(int start_coord, int seqlen) :
        ArgHmmForwardTable(start_coord, seqlen)
    {}

    virtual ~ArgHmmForwardTableOld() {
        if (fw) {
            for (int i=0; i<seqlen; i++)
                delete [] fw[i];
            delete [] fw;
            fw = NULL;
        }
    }

    // allocate another block of the forward table
    virtual void new_block(int start, int end, int nstates)
    {
        // allocate block
        for (int i=start; i<end; i++)
            fw[i-start_coord] = new double [max(nstates, 1)];
    }
};


//=============================================================================
// Forward algorithm for thread path

void arghmm_forward_alg(const LocalTrees *trees, const ArgModel *model,
    const Sequences *sequences, ArgHmmMatrixIter *matrix_iter,
    ArgHmmForwardTable *forward, PhaseProbs *phase_pr=NULL,
    bool prior_given=false, bool internal=false, bool slow=false);

double stochastic_traceback(
    const LocalTrees *trees, const ArgModel *model,
    ArgHmmMatrixIter *matrix_iter,
    double **fw, int *path, bool last_state_given=false, bool internal=false);

//=============================================================================
// ARG thread sampling

void sample_arg_thread(
    const ArgModel *model, Sequences *sequences, LocalTrees *trees,
    int new_chrom);

void sample_arg_thread_internal(
   const ArgModel *model, const Sequences *sequences, LocalTrees *trees,
   int minage=0, PhaseProbs *phase_pr=NULL);

void resample_arg_thread(const ArgModel *model, Sequences *sequences,
    LocalTrees *trees, int chrom);

void cond_sample_arg_thread(const ArgModel *model, const Sequences *sequences,
                            LocalTrees *trees, int new_chrom,
                            State start_state, State end_state);

void cond_sample_arg_thread_internal(
    const ArgModel *model, const Sequences *sequences,
    LocalTrees *trees, State start_state, State end_state);


} // namespace argweaver

#endif // ARGWEAVER_SAMPLE_THREAD_H

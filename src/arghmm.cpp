
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
#include "model.h"
#include "ptree.h"
#include "seq.h"



using namespace std;

using namespace spidir;
using namespace dlcoal;


namespace arghmm {




//=============================================================================
// transitions


void calc_transition_probs(LocalTree *tree, ArgModel *model,
                           const States &states, LineageCounts *lineages,
                           double **transprob)
{
    // get model parameters
    const int ntimes = model->ntimes;
    const double *times = model->times;
    const double *time_steps = model->time_steps;
    const double *popsizes = model->popsizes;
    const double rho = model->rho;

    const int nstates = states.size();

    // get tree information
    LocalNode *nodes = tree->nodes;
    const double treelen = get_treelen(tree, times, ntimes);
    const int *nbranches = lineages->nbranches;
    const int *nrecombs = lineages->nrecombs;
    const int *ncoals = lineages->ncoals;


    // find root node
    int root = tree->root;
    const int root_age_index = nodes[root].age;
    const double root_age = times[root_age_index];
    
    // C_j = C_{j-1} + s'_{j-1} k_{j-1} / (2N)

    // B_{c,a} =& \sum_{k=0}^{c} \exp(- A_{k,a})
    //         =& B_{c-1,a} + \exp(- A_{c,a}).

    // S_{a,b} &= exp(C_b) * B_{min(a,b),b}
    double C[ntimes];
    double B[ntimes];
    double eC[ntimes];
    double D[ntimes];
    
    C[0] = 0.0;
    eC[0] = 1.0;
    B[0] = nbranches[0] * time_steps[0] / (nrecombs[0] + 1.0);
    D[0] = (1.0 - exp(-time_steps[0] * nbranches[0]
                          / (2.0 * popsizes[0]))) / ncoals[0];

    for (int b=1; b<ntimes; b++) {
        const int l = b - 1;
        C[b] = C[l] + time_steps[l] * nbranches[l] / (2.0 * popsizes[l]);
        eC[b] = exp(C[b]);
        B[b] = B[b-1] + nbranches[b] * time_steps[b] / 
            (nrecombs[b] + 1.0) * eC[b];
        D[b] = (1.0 - exp(-time_steps[b] * nbranches[b]
                          / (2.0 * popsizes[b]))) / ncoals[b];
    }

    
    for (int i=0; i<nstates; i++) {
        const int node1 = states[i].node;
        const int a = states[i].time;
        const int c = nodes[node1].age;
        
        double treelen2 = treelen + times[a];
        if (a > root_age_index) {
            treelen2 += times[a] - root_age;
            treelen2 += time_steps[a]; // add basal branch
        } else {
            treelen2 += time_steps[root_age_index]; // add basal branch
        }

        double const F = (1.0 - exp(-rho * treelen2)) / treelen2;
        
        for (int j=0; j<nstates; j++) {
            const int node2 = states[j].node;
            const int b = states[j].time;
            
            const double f = F * D[b];
            
            if (node1 != node2)
                transprob[i][j] = f * B[min(a,b)] / eC[b];
            else {
                transprob[i][j] = f * (2 * B[min(a,b)] - B[min(c,b)]) / eC[b];
                if (a == b)
                    transprob[i][j] += exp(-rho * treelen2);
            }
        }

        // normalize and convert to log scale
        double sum = 0.0;
        for (int j=0; j<nstates; j++)
            sum += transprob[i][j];
        for (int j=0; j<nstates; j++)
            transprob[i][j] = log(transprob[i][j] / sum);
    }
}


void calc_transition_probs_compressed(
    LocalTree *tree, ArgModel *model, LineageCounts *lineages,
    double **transprob)
{
    // get model parameters
    const int ntimes = model->ntimes;
    const double *times = model->times;
    const double *time_steps = model->time_steps;
    const double *popsizes = model->popsizes;
    const double rho = model->rho;
    
    // get tree information
    LocalNode *nodes = tree->nodes;
    const double treelen = get_treelen(tree, times, ntimes);
    const int *nbranches = lineages->nbranches;
    const int *nrecombs = lineages->nrecombs;
    const int *ncoals = lineages->ncoals;


    // find root node
    int root = tree->root;
    const int root_age_index = nodes[root].age;
    const double root_age = times[root_age_index];
    
    // C_j = C_{j-1} + s'_{j-1} k_{j-1} / (2N)

    // B_{c,a} =& \sum_{k=0}^{c} \exp(- A_{k,a})
    //         =& B_{c-1,a} + \exp(- A_{c,a}).

    // S_{a,b} &= exp(C_b) * B_{min(a,b),b}
    double C[ntimes];
    double B[ntimes];
    double eC[ntimes];
    double D[ntimes];
    
    C[0] = 0.0;
    eC[0] = 1.0;
    B[0] = nbranches[0] * time_steps[0] / (nrecombs[0] + 1.0);
    D[0] = (1.0 - exp(-time_steps[0] * nbranches[0]
                          / (2.0 * popsizes[0]))) / ncoals[0];

    for (int b=1; b<ntimes; b++) {
        const int l = b - 1;
        C[b] = C[l] + time_steps[l] * nbranches[l] / (2.0 * popsizes[l]);
        eC[b] = exp(C[b]);
        B[b] = B[b-1] + nbranches[b] * time_steps[b] / 
            (nrecombs[b] + 1.0) * eC[b];
        D[b] = (1.0 - exp(-time_steps[b] * nbranches[b]
                          / (2.0 * popsizes[b]))) / ncoals[b];
    }

    
    for (int a=0; a<ntimes; a++) {
        double treelen2 = treelen + times[a];
        if (a > root_age) {
            treelen2 += times[a] - root_age;
            treelen2 += time_steps[a]; // add basal branch
        } else {
            // TODO: is this double counted?
            treelen2 += time_steps[root_age_index]; // add basal branch
        }
        
        double const F = (1.0 - exp(-rho * treelen2)) /
            (exp(-rho * treelen) * treelen2);
        
        for (int b=0; b<ntimes; b++)
            transprob[a][b] = log(F * D[b] * B[min(a,b)] / eC[b]);
    }
}



void get_deterministic_transitions(
    LocalTree *tree, LocalTree *last_tree, const Spr &spr, int *mapping,
    const States &states1, const States &states2,
    int ntimes, int *next_states)
{
    // recomb_node in tree and last_tree
    // coal_node in last_tree
    
    const int nnodes = tree->nnodes;
    const LocalNode *nodes = tree->nodes;
    const LocalNode *last_nodes = last_tree->nodes;
    const int nstates1 = states1.size();

    // make state lookup
    NodeStateLookup state2_lookup(states2, nnodes);
    

    for (int i=0; i<nstates1; i++) {
        const int node1 = states1[i].node;
        const int time1 = states1[i].time;

        if (node1 == spr.coal_node && time1 == spr.coal_time) {
            // not a deterministic case
            next_states[i] = -1;
        
        } else if (node1 != spr.recomb_node) {
            // SPR only removes a subset of descendents, if any
            // trace up from remaining leaf to find correct new state

            int node2;
            const LocalNode *node = &last_nodes[node1];
            bool disrupt = false;
            
            if (node->child[0] == -1) {
                // SPR can't disrupt leaf branch
                node2 = node1;

            } else {
                int child1 = node->child[0];
                int child2 = node->child[1];
                
                if (spr.recomb_node == child1) {
                    // right child is not disrupted
                    node2 = mapping[child2];
                    disrupt = true;
                } else if (spr.recomb_node == child2) {
                    // left child is not disrupted
                    node2 = mapping[child1];
                    disrupt = true;
                } else {
                    // node is not disrupted
                    node2 = mapping[node1];
                }
            }

            // optionally walk up
            if ((spr.coal_node == node1 && spr.coal_time < time1) || 
                (mapping[spr.coal_node] == node2 && spr.coal_time < time1) ||
                (disrupt && mapping[spr.coal_node] == node2 && 
                 spr.coal_time <= time1))
            {
                // XXX: need to walk up for one more case
                // need to change < time1 to <== time for disrupt case
                // coal occurs under us
                node2 = nodes[node2].parent;
            }

            // DEBUG
            assert(nodes[node2].age <= time1);
            int p = nodes[node2].parent;
            if (p != -1)
                assert(time1 <= nodes[p].age);
            
            // set next state
            next_states[i] = state2_lookup.lookup(node2, time1);

        } else {
            // SPR is on same branch as new chromosome
            if (spr.recomb_time >= time1) {
                // we move with SPR subtree
                // TODO: we could probabilistically have subtree move
                // out from underneath.
                next_states[i] = state2_lookup.lookup(
                    mapping[spr.recomb_node], time1);

            } else {
                // SPR should not be able to coal back onto same branch
                // this would be a self cycle
                assert(spr.coal_node != node1);
                
                // SPR subtree moves out from underneath us
                // therefore therefore the new chromosome coalesces with
                // the branch above the subtree

                // search up for parent
                int parent = last_nodes[spr.recomb_node].parent;
                int time2 = last_nodes[parent].age;
                int node2;

                // find other child
                const int *c = last_nodes[parent].child;
                int other = (c[1] == spr.recomb_node ? c[0] : c[1]);

                // find new state in tree
                node2 = (other == spr.coal_node ? 
                         nodes[mapping[other]].parent : mapping[other]);
                next_states[i] = state2_lookup.lookup(node2, time2);
            }
        }
    }
    
}



void calc_transition_probs_switch(
    LocalTree *tree, LocalTree *last_tree, const Spr &spr, int *mapping,
    const States &states1, const States &states2,
    ArgModel *model, LineageCounts *lineages, double **transprob)
{
    // get tree information
    const LocalNode *nodes = tree->nodes;
    const LocalNode *last_nodes = last_tree->nodes;
    const int *nbranches = lineages->nbranches;
    const int *ncoals = lineages->ncoals;

    const int nstates1 = states1.size();
    const int nstates2 = states2.size();
    
    // get model parameters
    const int ntimes = model->ntimes;
    const double *time_steps = model->time_steps;
    const double *popsizes = model->popsizes;
    

    // get deterministic transitions
    int determ[nstates1];
    get_deterministic_transitions(tree, last_tree, spr, mapping,
                                  states1, states2, ntimes, determ);
    

    for (int i=0; i<nstates1; i++) {
        const int node1 = states1[i].node;
        const int time1 = states1[i].time;
        
        if (node1 != spr.coal_node || time1 != spr.coal_time) {
            // deterministic transition case
            assert (determ[i] != -1);
            for (int j=0; j<nstates2; j++)
                transprob[i][j] = -INFINITY;
            transprob[i][determ[i]] = 0.0;
            
        } else {
            // probabilistic transition case
            
            // determine if node1 is still here or not
            int node3;
            int last_parent = last_nodes[spr.recomb_node].parent;
            if (last_parent == node1) {
                // recomb breaks node1 branch, we need to use the other child
                const int *c = last_nodes[last_parent].child;
                node3 = mapping[c[1] == spr.recomb_node ? c[0] : c[1]];
            } else {
                node3 = mapping[node1];
            }

            // find parent of recomb_branch and node1
            int last_parent_age = last_nodes[last_parent].age;

            
            int parent = nodes[mapping[spr.recomb_node]].parent;
            assert(parent == nodes[node3].parent);

            for (int j=0; j<nstates2; j++) {
                const int node2 = states2[j].node;
                const int time2 = states2[j].time;
                
                transprob[i][j] = 0.0;
                if (!((node2 == mapping[spr.recomb_node] 
                       && time2 >= spr.recomb_time) ||
                      (node2 == node3 && time2 == time1) ||
                      (node2 == parent && time2 == time1)))
                    // not a probabilistic transition
                    continue;

                // get lineage counts
                // remove recombination branch and add new branch
                int kbn = nbranches[time2];
                int kcn = ncoals[time2] + 1;
                if (time2 < nodes[parent].age) {
                    kbn -= 1;
                    kcn -= 1;
                }
                if (time2 < time1)
                    kbn += 1;

                double sum = 0.0;
                for (int m=spr.recomb_time; m<time2; m++) {
                    const int nlineages = nbranches[m] + 1
                        - (m < last_parent_age ? 1 : 0);
                    sum += time_steps[m] * nlineages / (2.0 * popsizes[m]);
                }
                transprob[i][j] =
                    (1.0 - exp(- time_steps[time2] * kbn /  
                               (2.0 * popsizes[time2]))) / kcn * exp(-sum);
            }

            // normalize row to ensure they add up to one
            double sum = 0.0;
            for (int j=0; j<nstates2; j++)
                sum += transprob[i][j];
            for (int j=0; j<nstates2; j++) {
                double x = transprob[i][j];
                if (sum > 0.0 and x > 0.0)
                    transprob[i][j] = log(x / sum);
                else
                    transprob[i][j] = -INFINITY;
            }
        }
    }
}



void calc_state_priors(const States &states, LineageCounts *lineages, 
                       ArgModel *model, double *priors)
{
    const int nstates = states.size();
    const double *time_steps = model->time_steps;
    const double *popsizes = model->popsizes;
    const int *nbranches = lineages->nbranches;
    const int *ncoals = lineages->ncoals;
    
    for (int i=0; i<nstates; i++) {
        int b = states[i].time;

        double sum = 0.0;
        for (int m=0; m<b; m++)
            sum += time_steps[m] * nbranches[m] / (2.0 * popsizes[m]);
        
        priors[i] = log((1.0 - exp(- time_steps[b] * nbranches[b] /
                          (2.0 * popsizes[b]))) / ncoals[b] * exp(-sum)); 
    }
}


//=============================================================================
// C interface
extern "C" {

double **new_transition_probs(int nnodes, int *ptree, 
                              int *ages, double treelen,
                              intstate *istates, int nstates,
                              int ntimes, double *times, double *time_steps,
                              int *nbranches, int *nrecombs, int *ncoals, 
                              double *popsizes, double rho)
{

    // setup model, local tree, states
    ArgModel model(ntimes, times, popsizes, rho, 0.0);
    LocalTree tree(ptree, nnodes, ages);
    LineageCounts lineages(ntimes);
    lineages.count(&tree);
    States states;
    make_states(istates, nstates, states);

    double **transprob = new_matrix<double>(nstates, nstates);
    calc_transition_probs(&tree, &model, states, &lineages, transprob);
    return transprob;
}


double **new_transition_probs_switch(
    int *ptree, int *last_ptree, int nnodes, 
    int recomb_node, int recomb_time, int coal_node, int coal_time,
    int *ages_index, int *last_ages_index,
    double treelen, double last_treelen,
    intstate *istates1, int nstates1,
    intstate *istates2, int nstates2,
    
    int ntimes, double *times, double *time_steps,
    int *nbranches, int *nrecombs, int *ncoals, 
    double *popsizes, double rho)
{
    // setup model
    ArgModel model(ntimes, times, popsizes, rho, 0.0);
    
    // setup local trees
    LocalTree tree(ptree, nnodes, ages_index);
    LocalTree last_tree(last_ptree, nnodes, last_ages_index);
    Spr spr(recomb_node, recomb_time, coal_node, coal_time);
    int mapping[nnodes];
    make_node_mapping(mapping, nnodes, last_ptree, recomb_node);
    LineageCounts lineages(ntimes);
    lineages.count(&last_tree);
    
    // setup states
    States states1, states2;
    make_states(istates1, nstates1, states1);
    make_states(istates2, nstates2, states2);
    
    double **transprob = new_matrix<double>(nstates1, nstates2);
    calc_transition_probs_switch(&tree, &last_tree, spr, mapping,
        states1, states2, &model, &lineages, transprob);
    return transprob;
}


void delete_transition_probs(double **transmat, int nstates)
{
    delete_matrix<double>(transmat, nstates);
}

} // extern C




//=============================================================================

class ArgHmmMatrices
{
public:
    ArgHmmMatrices(int nstates1, int nstates2, int blocklen,
                   double **transmat, double **transmat_switch, double **emit):
        nstates1(nstates1),
        nstates2(nstates2),
        blocklen(blocklen),
        transmat(transmat),
        transmat_switch(transmat_switch),
        emit(emit)
    {}
    ~ArgHmmMatrices() 
    {}

    void clear()
    {
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
    double **transmat;
    double **transmat_switch;
    double **emit;
};


class ArgHmmMatrixList
{
public:
    ArgHmmMatrixList() {}
    ~ArgHmmMatrixList() {}
    
    void setup(ArgModel *model, Sequences *seqs, LocalTrees *trees)
    {
        setup(model, seqs, trees->begin(), trees->end());
    }

    void setup(ArgModel *_model, Sequences *_seqs, LocalTrees::iterator start, 
               LocalTrees::iterator end)
    {
        model = _model;
        char **seqs = NULL;
        if (_seqs)
            seqs = _seqs->seqs;
        int nseqs = start->tree->get_num_leaves() + 1;
    
        // allocate lineage counts
        LineageCounts lineages(model->ntimes);
    
        LocalTree *last_tree = NULL;
    
        // state spaces
        States states1;
        States states2;
        States *states = &states1;
        States *last_states = NULL;
        int nstates1, nstates2;
    
        // temporary subsequence of the sequences
        char *subseqs[nseqs];
        
        // matrices
        double **transmat_switch;
        double **transmat;
        
        // iterate over local trees
        for (LocalTrees::iterator it=start; it != end; it++) {
            int pos = it->block.start;
            int blocklen = it->block.end - it->block.start;
            LocalTree *tree = it->tree;
            get_coal_states(tree, model->ntimes, *states);
            int nstates = states->size();
                
            // calculate emissions
            double **emit = NULL;
            if (seqs) {
                for (int i=0; i<nseqs; i++)
                    subseqs[i] = &seqs[i][pos];
                emit = new_matrix<double>(blocklen, nstates);
                calc_emissions(*states, tree, subseqs, nseqs, blocklen, 
                               model, emit);
            }
        
            
            // use switch matrix for first column of forward table
            // if we have a previous state space (i.e. not first block)
            if (!last_states) {
                transmat_switch = NULL;
                nstates1 = nstates2 = nstates;
                lineages.count(tree);
            } else {
                nstates1 = last_states->size();
                nstates2 = states->size();
            
                // calculate transmat_switch
                lineages.count(last_tree);
                transmat_switch = new_matrix<double>(nstates1, nstates2);
            
                calc_transition_probs_switch(
                    tree, last_tree, it->spr, it->mapping,
                    *last_states, *states,
                    model, &lineages, transmat_switch);

                // update lineages to current tree
                lineages.count(tree);
            }
        
            // calculate transmat and use it for rest of block
            transmat = new_matrix<double>(nstates, nstates);
            calc_transition_probs(tree, model, *states, &lineages, transmat);

            // store matrices
            matrices.push_back(ArgHmmMatrices(nstates1, nstates2, blocklen,
                transmat, transmat_switch, emit));


            // update pointers
            last_tree = tree;
            last_states = states;
            states = ((states == &states1) ? &states2 : &states1);
        }

        
    }

    void clear()
    {
        for (unsigned int i=0; i<matrices.size(); i++)
            matrices[i].clear();
        matrices.clear();
    }

    ArgModel *model;
    vector<ArgHmmMatrices> matrices;
};


//=============================================================================
// HMM algorithms




void arghmm_forward_alg_block(
    int nleaves, int n, const States &states, int ntimes,
    double **trans, double **emit, double **fw)
{
    const int nstates = states.size();
    double vec[nstates];
    
    for (int i=1; i<n; i++) {
        double *col1 = fw[i-1];
        double *col2 = fw[i];
        double *emit2 = emit[i];

        for (int k=0; k<nstates; k++) {
            for (int j=0; j<nstates; j++)
                vec[j] = col1[j] + trans[j][k];
            col2[k] = logsum(vec, nstates) + emit2[k];
        }
    }
}


void arghmm_forward_alg_block(
    int n, int nstates, 
    double **trans, double **emit, double **fw)
{
    double vec[nstates];
    
    for (int i=1; i<n; i++) {
        double *col1 = fw[i-1];
        double *col2 = fw[i];
        double *emit2 = emit[i];

        for (int k=0; k<nstates; k++) {
            for (int j=0; j<nstates; j++)
                vec[j] = col1[j] + trans[j][k];
            col2[k] = logsum(vec, nstates) + emit2[k];
        }
    }
}



void arghmm_forward_alg_block3(
    int nleaves, int n, const States &states, int ntimes,
    double **trans, double **emit, double **fw)
{
    const int nstates = states.size();

    if (nleaves == 1)
        return forward_alg(n, nstates, trans, emit, fw);

    double vec[nstates];
    double fgroups[ntimes];
    int time2state[ntimes];
    //double rest[nstates];
    //int restlen = 0;

    for (int i=1; i<n; i++) {
        double *col1 = fw[i-1];
        double *col2 = fw[i];
        double *emit2 = emit[i];

        // precompute the fgroup sums
        memset(fgroups, sizeof(double) * ntimes, 0);
        for (int j=0; j<nstates; j++) {
            int a = states[j].time;
            fgroups[a] += col1[j];
            time2state[a] = j;
        }

        for (int k=0; k<nstates; k++) {
            // iterate only over the fgroups
            for (int a=0; a<ntimes; a++)
                vec[a] = fgroups[a] + trans[time2state[a]][k];

            // collect terms that don't change branch

            col2[k] = logsum(vec, ntimes) + emit2[k];
        }
    }
}



// run forward algorithm with low memory for matrices
double **arghmm_forward_alg(LocalTrees *trees, ArgModel *model,
                            Sequences *sequences, double **fw=NULL)
{
    LocalTree *last_tree = NULL;
    const int nleaves = (trees->nnodes + 1) / 2;
    const int ntimes = model->ntimes;
    

    // allocate lineage counts
    LineageCounts lineages(ntimes);
   
    // state spaces
    States states1;
    States states2;
    States *states = &states1;
    States *last_states = NULL;
    
    // temporary subsequence of the sequence
    char *subseqs[sequences->nseqs];
    
    // allocate the forward table if necessary
    if (fw == NULL) {
        fw = new double* [sequences->seqlen];
        for (int i=0; i<sequences->seqlen; i++)
            fw[i] = NULL;
    }
    

    // iterate over local trees to find largest blocklen and nstates
    int max_nstates = 0;
    int max_blocklen = 0;
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); it++) {
        int blocklen = it->block.end - it->block.start;
        get_coal_states(it->tree, ntimes, *states);
        int nstates = states->size();
        max_nstates = max(max_nstates, nstates);
        max_blocklen = max(max_blocklen, blocklen);
    }


    // allocate HMM matrices
    double **emit = new_matrix<double>(max_blocklen, max_nstates);
    double **transmat = new_matrix<double>(max_nstates, max_nstates);
    double **transmat_switch = new_matrix<double>(max_nstates, max_nstates);
    
    
    // iterate over local trees
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); it++) {
        int pos = it->block.start;
        int blocklen = it->block.end - it->block.start;
        LocalTree *tree = it->tree;
        get_coal_states(tree, ntimes, *states);
        int nstates = states->size();
        
        // allocate the forward table column if necessary
        if (fw[pos] == NULL) {
            for (int i=pos; i<pos+blocklen; i++) {
                fw[i] = new double [nstates];
            }
        }
        
        // calculate emissions
        for (int i=0; i<sequences->nseqs; i++)
            subseqs[i] = &sequences->seqs[i][pos];
        calc_emissions(*states, tree, subseqs, sequences->nseqs, 
                       blocklen, model, emit);
        
        
        // use switch matrix for first column of forward table
        // if we have a previous state space (i.e. not first block)
        if (!last_states) {
            // calculate prior of first state
            lineages.count(tree);
            calc_state_priors(*states, &lineages, model, fw[0]);  

        } else {
            // perform one column of forward algorithm with transmat_switch
            lineages.count(last_tree);
            calc_transition_probs_switch(tree, last_tree, it->spr, it->mapping,
                                         *last_states, *states,
                                         model, &lineages, transmat_switch);
            forward_step(fw[pos-1], fw[pos], last_states->size(), 
                         states->size(), transmat_switch, emit[0]);

            // update lineages to current tree
            lineages.count(tree);
        }
        
        // calculate transmat and use it for rest of block
        calc_transition_probs(tree, model, *states, &lineages, transmat);
        arghmm_forward_alg_block(nleaves, blocklen, *states, ntimes, 
                                 transmat, emit, &fw[pos]);
        
        // update pointers
        last_tree = tree;
        last_states = states;
        states = ((states == &states1) ? &states2 : &states1);
    }

    // cleanup
    delete_matrix<double>(emit, max_blocklen);
    delete_matrix<double>(transmat, max_nstates);
    delete_matrix<double>(transmat_switch, max_nstates);


    return fw;
}


// run forward algorithm with matrices precomputed
double **arghmm_forward_alg(LocalTrees *trees, ArgModel *model,
    Sequences *sequences, ArgHmmMatrixList *matrix_list, double **fw=NULL)
{
    LineageCounts lineages(model->ntimes);
    States states;
    
    // forward algorithm over local trees
    int blocki = 0;
    for (LocalTrees::iterator it=trees->begin(); 
         it != trees->end(); it++, blocki++) {
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
        arghmm_forward_alg_block(matrices.blocklen, matrices.nstates2, 
                                 matrices.transmat, matrices.emit, &fw[pos]);
    }
    
    return fw;
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
                             mat.emit, &fw[pos], &path[pos]);
        
        // use switch matrix for last col of next block
        if (pos > 0) {
            int i = pos - 1;
            double A[mat.nstates1];
            int k = path[i+1];
            for (int j=0; j<mat.nstates1; j++)
                A[j] = fw[i][j] + mat.transmat_switch[j][k];
            double total = logsum(A, mat.nstates1);
            for (int j=0; j<mat.nstates1; j++)
                A[j] = exp(A[j] - total);
            path[i] = sample(A, mat.nstates1);
            
            printf("trace %d, %d %d, %e\n", pos, path[i], k,
                   mat.transmat_switch[path[i]][k]);
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
    
    // loop through local blocks
    int blocki = 0;
    for (LocalTrees::iterator it=trees->begin(); 
         it != trees->end(); it++, blocki++) {

        // get local block information
        int start = it->block.start + 1;  // don't allow new recomb at start
        int end = it->block.end;
        ArgHmmMatrices matrices = matrix_list->matrices[blocki];
        LocalTree *tree = it->tree;
        double treelen = get_treelen(tree, model->times, model->ntimes);
        double treelen2;
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
                    treelen2 = get_treelen_branch(
                        tree, model->times, model->ntimes,
                        states[last_state].node,
                        states[last_state].time, treelen);
                    double self_trans = matrices.transmat[last_state][last_state];
                    double rate = max(
                        1.0 - exp(-model->rho * (treelen2 - treelen)
                              - self_trans), model->rho);

                    next_recomb = int(min(float(end), i + expovariate(rate)));

                    //printf("next %d %e %d\n", i, 
                    //       expovariate(rate), next_recomb);
                }

                if (i < next_recomb)
                    continue;
            }

            // sample recombination
            next_recomb = -1;
            statei = thread_path[i];
            State state = states[statei];
            State last_state = states[thread_path[i-1]];
            treelen2 = get_treelen_branch(tree, model->times, model->ntimes,
                                          last_state.node,
                                          last_state.time, treelen);
            

            // there must be a recombination
            // either because state changed or we choose to recombine
            // find candidates
            vector <NodePoint> candidates;
            vector <double> probs;
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

            int j = state.time;
            for (vector<NodePoint>::iterator it=candidates.begin(); 
                 it != candidates.end(); ++it) {
                int k = it->time;
                double p = 
                    (lineages.nbranches[k] + 1) * model->time_steps[k] /
                    (lineages.ncoals[j] * (lineages.nrecombs[k] + 1.0) * 
                     treelen2) * 
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


void apply_spr(LocalTree *tree, Spr *spr)
{
    // before SPR:
    //       bp          cp
    //      / \           \       .
    //     rc              c
    //    / \                     .
    //   r   rs

    // after SPR:
    //    bp         cp
    //   /  \         \           .
    //  rs             rc
    //                /  \        .
    //               r    c

    // key:
    // r = recomb branch
    // rs = sibling of recomb branch
    // rc = recoal node (broken node)
    // bp = parent of broken node
    // c = coal branch
    // cp = parent of coal branch


    LocalNode *nodes = tree->nodes;

    // recoal is also the node we are breaking
    int recoal = nodes[spr->recomb_node].parent;

    // find recomb node sibling and broke node parent
    int *c = nodes[recoal].child;
    int other = (c[0] == spr->recomb_node ? 1 : 0);
    int recomb_sib = c[other];
    int broke_parent =  nodes[recoal].parent;


    // fix recomb sib pointer
    nodes[recomb_sib].parent = broke_parent;

    // fix parent of broken node
    int x = 0;
    if (broke_parent != -1) {
        c = nodes[broke_parent].child;
        x = (c[0] == recoal ? 0 : 1);
        nodes[broke_parent].child[x] = recomb_sib;
    }

    // reuse node as recoal
    if (spr->coal_node == recoal) {
        // we just broke coal_node, so use recomb_sib
        nodes[recoal].child[other] = recomb_sib;
        nodes[recoal].parent = nodes[recomb_sib].parent;
        nodes[recomb_sib].parent = recoal;
        if (broke_parent != -1)
            nodes[broke_parent].child[x] = recoal;
    } else {
        nodes[recoal].child[other] = spr->coal_node;
        nodes[recoal].parent = nodes[spr->coal_node].parent;
        nodes[spr->coal_node].parent = recoal;
        
        // fix coal_node parent
        int parent = nodes[recoal].parent;
        if (parent != -1) {
            c = nodes[parent].child;
            if (c[0] == spr->coal_node) 
                c[0] = recoal;
            else
                c[1] = recoal;
        }
    }
    nodes[recoal].age = spr->coal_time;   
    
    // set tree data
    tree->set_root();
}


// add a thread to an ARG
bool assert_trees(LocalTrees *trees)
{
    LocalTree *last_tree = NULL;

    // loop through blocks
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it) {
        printf("check %d,%d\n", it->block.start, it->block.end);
        LocalTree *tree = it->tree;
        Spr *spr = &it->spr;
        int *mapping = it->mapping;
        if (last_tree)
            assert(assert_spr(last_tree, tree, spr, mapping));
        last_tree = tree;
    }
    
    return true;
}


// add a thread to an ARG
bool assert_trees_thread(LocalTrees *trees, int *thread_path, int ntimes)
{
    LocalTree *last_tree = NULL;
    States states1, states2;
    States *states = &states1;
    States *last_states = &states2;

    // loop through blocks
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it) {
        int start = it->block.start;
        get_coal_states(it->tree, ntimes, *states);

        // check spr
        if (last_tree) {
            assert(assert_spr(last_tree, it->tree, &it->spr, it->mapping));

            int determ[last_states->size()];
            get_deterministic_transitions(
                 it->tree, last_tree, it->spr, it->mapping,
                 *last_states, *states, ntimes, determ);

            int a = thread_path[start-1];
            int b = determ[thread_path[start-1]];
            int c = thread_path[start];

            printf(">> %d (%d,%d) --> (%d,%d) == (%d,%d)\n", start, 
                   (*last_states)[a].node, (*last_states)[a].time, 
                   (*states)[b].node, (*states)[b].time,
                   (*states)[c].node, (*states)[c].time);
        }
        
        // set last tree and state pointers
        last_tree = it->tree;
        last_states = states;
        if (states == &states1)
            states = &states2;
        else
            states = &states1;
    }

    return true;
}



// add a thread to an ARG
void add_arg_thread(LocalTrees *trees, int ntimes, int *thread_path, 
                    vector<int> &recomb_pos, vector<NodePoint> &recombs)
{
    unsigned int irecomb = 0;
    int nleaves = trees->get_num_leaves();
    int nnodes = trees->nnodes;
    int nnodes2 = nnodes + 2;
    
    int newleaf = nleaves;
    int displaced = nnodes;
    int newcoal = nnodes + 1;

    States states;
    LocalTree *last_tree = NULL;
    int last_node = -1;


    // loop through blocks
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it) {
        LocalTree *tree = it->tree;
        Spr *spr = &(it->spr);
        const int start = it->block.start;
        const int end = it->block.end;
        get_coal_states(tree, ntimes, states);

        // DBEUG
        printf("add %d\n", start);        
        
        // add new branch to local tree
        it->ensure_capacity(nnodes2);
        LocalNode *nodes = tree->nodes;
        State state = states[thread_path[start]];
        
        // determine node displacement
        int node2 = (state.node != newleaf ? state.node : displaced);
        int parent = nodes[state.node].parent;
        int parent2 = (parent != newleaf ? parent : displaced);

        // displace node
        if (newleaf < displaced) {
            nodes[displaced] = nodes[newleaf]; // copy displaced node
            if (nodes[displaced].parent != -1) {
                int *c = nodes[nodes[displaced].parent].child;
                if (c[0] == newleaf)
                    c[0] = displaced;
                else
                    c[1] = displaced;
            }
            int *c = nodes[displaced].child;
            nodes[c[0]].parent = displaced;
            nodes[c[1]].parent = displaced;
        }
   

        // add new leaf
        nodes[newleaf].parent = newcoal;
        nodes[newleaf].child[0] = -1;
        nodes[newleaf].child[1] = -1;
        nodes[newleaf].age = 0;
        
        // add new coal node
        nodes[newcoal].parent = parent2;
        nodes[newcoal].child[0] = newleaf;
        nodes[newcoal].child[1] = node2;
        nodes[newcoal].age = state.time;

        // fix pointers
        nodes[node2].parent = newcoal;
        if (parent2 != -1) {
            int *c = nodes[parent2].child;
            if (c[0] == node2)
                c[0] = newcoal;
            else
                c[1] = newcoal;
        }

        // fix up tree data
        tree->nnodes = nnodes2;
        tree->set_root();
        assert(assert_tree(tree));


        // update mapping and spr
        int *mapping = it->mapping;
        if (mapping) {
            // update mapping due to displacement
            mapping[displaced] = mapping[newleaf];
            mapping[newleaf] = newleaf;
            
            // set default new node mapping 
            mapping[newcoal] = newcoal;

            for (int i=newleaf+1; i<tree->nnodes; i++) {
                if (mapping[i] == newleaf)
                    mapping[i] = displaced;
            }


            // update spr due to displacement
            if (spr->recomb_node == newleaf)
                spr->recomb_node = displaced;
            if (spr->coal_node == newleaf)
                spr->coal_node = displaced;

            LocalNode *last_nodes = last_tree->nodes;
        
            // parent of recomb node should be the recoal point
            // however, if it equals newcoal, then recomb branch is 
            // either renamed or we have mediation
            int recoal = nodes[mapping[spr->recomb_node]].parent;
            if (recoal == newcoal) {
                if (mapping[last_node] == node2) {
                    // recomb is above coal state, we rename spr recomb node
                    spr->recomb_node = newcoal;
                } else {
                    // this is a mediated coal, rename coal node and time
                    spr->coal_node = newleaf;
                    spr->coal_time = state.time;
                }
            } else {
                // the other possibility is that newcoal is under recoal point
                // if newcoal is child of recoal, then coal is renamed
                int *c = nodes[recoal].child;
                if (c[0] == newcoal || c[1] == newcoal) {
                    // we either coal above the newcoal or our existing
                    // node just broke and newcoal was underneath.

                    // if newcoal was previously above spr->coal_node
                    // then we rename the spr coal node
                    if (last_nodes[spr->coal_node].parent == newcoal)
                        spr->coal_node = newcoal;
                }
            }
            
            // determine if mapping of new node needs to be changed
            // newcoal was parent of recomb, it is broken
            if (last_nodes[spr->recomb_node].parent == newcoal) {
                mapping[newcoal] = -1;
                int p = last_nodes[newcoal].parent;
                if (p != -1)
                    mapping[p] = newcoal;
            } else {
                // newcoal was not broken
                // find child without recomb or coal on it
                int x = newcoal;
                while (true) {
                    int y = last_nodes[x].child[0];
                    printf("x=%d, y=%d\n", x, y);
                    if (y == spr->coal_node || y == spr->recomb_node)
                        y = last_nodes[x].child[1];
                    x = y;
                    printf("  x=%d, mapping[x] = %d\n", x, mapping[x]);
                    if (mapping[x] != -1)
                        break;
                }
                printf(":: %d -> %d, x=%d\n", 
                       newcoal, nodes[mapping[x]].parent, x);
                mapping[newcoal] = nodes[mapping[x]].parent;
            }

            // assert SPR
            if (!assert_spr(last_tree, tree, spr, mapping)) {
                printf("!!! %d spr fail\n", start);
            }
        }

        // assert new branch is where it should be
        assert(tree->nodes[newcoal].age == states[thread_path[start]].time);


        // break this block for each new recomb within this block
        for (;irecomb < recombs.size() && 
              recomb_pos[irecomb] < end; irecomb++) {
            int pos = recomb_pos[irecomb];
            printf("start %d pos %d\n", start, pos);

            assert(tree->nodes[newcoal].age == states[thread_path[pos-1]].time);

            // determine real name of recomb node
            // it may be different due to displacement
            Spr spr2;
            spr2.recomb_node = recombs[irecomb].node;
            spr2.recomb_time = recombs[irecomb].time;
            if (spr2.recomb_node == newleaf)
                spr2.recomb_node = displaced;            
            assert(spr2.recomb_time <= tree->nodes[newcoal].age);

            // determine coal node and time
            int istate = thread_path[pos];
            if (spr2.recomb_node == -1) {
                // recomb on new branch, coal given thread
                spr2.recomb_node = newleaf;
                spr2.coal_node = states[istate].node;

                // fix coal node due to displacement
                if (spr2.coal_node == newleaf)
                    spr2.coal_node = displaced;

                // rename due to newcoal
                if (states[istate].node == states[thread_path[pos-1]].node &&
                    states[istate].time > states[thread_path[pos-1]].time)
                    spr2.coal_node = newcoal;

            } else {
                // recomb in ARG, coal on new branch
                // fix recomb node due to displacement
                if (spr2.recomb_node == newleaf)
                    spr2.recomb_node = displaced;

                if (states[istate].time > states[thread_path[pos-1]].time)
                    spr2.coal_node = nodes[newleaf].parent;
                else
                    spr2.coal_node = newleaf;
            }
            spr2.coal_time = states[istate].time;

            // determine mapping
            int *mapping2 = new int [tree->capacity];
            for (int j=0; j<nnodes2; j++)
                mapping2[j] = j;
            mapping2[nodes[spr2.recomb_node].parent] = -1;

            // make new local tree
            LocalTree *new_tree = new LocalTree(nnodes2, tree->capacity);
            LocalNode *new_nodes = new_tree->nodes;
            for (int j=0; j<nnodes2; j++) {
                new_nodes[j].parent = nodes[j].parent;
                new_nodes[j].child[0] = nodes[j].child[0];
                new_nodes[j].child[1] = nodes[j].child[1];
                new_nodes[j].age = nodes[j].age;
            }

            // apply spr to new tree
            apply_spr(new_tree, &spr2);

            int block_end;
            if (irecomb < recombs.size() - 1)
                block_end = min(recomb_pos[irecomb+1], end);
            else
                block_end = end;
            
            // assert tree and SPR
            assert(assert_tree(new_tree));
            assert(new_tree->nodes[newcoal].age == 
                   states[thread_path[pos]].time);
            if (!assert_spr(tree, new_tree, &spr2, mapping2)) {
                printf("!!! %d spr fail\n", pos);
            }


            // insert new tree into local trees list
            it->block.end = pos;
            ++it;
            it = trees->trees.insert(it, 
                LocalTreeSpr(pos, block_end, new_tree, spr2, mapping2));
            
            tree = new_tree;
            nodes = tree->nodes;
        }

        last_tree = tree;
        last_node = states[thread_path[end-1]].node;
        if (last_node == newleaf)
            last_node = displaced;
    }

    
    // update number of nodes
    trees->nnodes = nnodes2;

    assert_trees(trees);
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

    return arghmm_forward_alg(&trees, &model, &sequences, fw);
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
    ArgHmmMatrixList matrix_list;
    matrix_list.setup(&model, &sequences, &trees);
    
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
    matrix_list.clear();

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
    ArgHmmMatrixList matrix_list;
    matrix_list.setup(&model, &sequences, trees);
    
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
    assert_trees_thread(trees, thread_path, ntimes);
    add_arg_thread(trees, model.ntimes, thread_path, recomb_pos, recombs);

    // clean up
    delete_matrix<double>(fw, seqlen);
    delete [] thread_path;
    matrix_list.clear();

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
    
    LocalTrees *trees = new LocalTrees();
    const int capacity = 2 * nseqs - 1;
    trees->make_trunk(0, seqlen, capacity); 

    // allocate temp variables
    double **fw = new double* [seqlen];
    int *thread_path = new int [seqlen];

    // add more chromosomes one by one
    for (int nchroms=2; nchroms<=nseqs; nchroms++) {
        // use first nchroms sequences
        Sequences sequences(seqs, nchroms, seqlen);
    
        // build matrices
        ArgHmmMatrixList matrix_list;
        matrix_list.setup(&model, &sequences, trees);
    
        // compute forward table
        arghmm_forward_alg(trees, &model, &sequences, &matrix_list, fw);

        // traceback
        stochastic_traceback(&matrix_list, fw, thread_path, seqlen);

        // sample recombination points
        vector<int> recomb_pos;
        vector<NodePoint> recombs;
        sample_recombinations(trees, &model, &matrix_list,
                              thread_path, recomb_pos, recombs);

        // add thread to ARG
        add_arg_thread(trees, model.ntimes, thread_path, recomb_pos, recombs);

        // cleanup
        matrix_list.clear();
        for (int i=0; i<seqlen; i++)
            delete [] fw[i];
    }

    // clean up    
    delete [] thread_path;

    return trees;
}



int get_local_trees_ntrees(LocalTrees *trees)
{
    return trees->trees.size();
}


int get_local_trees_nnodes(LocalTrees *trees)
{
    return trees->nnodes;
}


void get_local_trees_ptrees(LocalTrees *trees, int **ptrees, int **ages,
                            int **sprs, int *blocklens)
{
    int i = 0;
    for (LocalTrees::iterator it=trees->begin(); it!=trees->end(); ++it, i++) {
        LocalTree *tree = it->tree;

        for (int j=0; j<tree->nnodes; j++) {
            ptrees[i][j] = tree->nodes[j].parent;
            ages[i][j] = tree->nodes[j].age;
            blocklens[i] = it->block.length();
            sprs[i][0] = it->spr.recomb_node;
            sprs[i][1] = it->spr.recomb_time;
            sprs[i][2] = it->spr.coal_node;
            sprs[i][3] = it->spr.coal_time;
        }
    }
}


void delete_local_trees(LocalTrees *trees)
{
    delete trees;
}


void delete_path(int *path)
{
    delete [] path;
}


void delete_double_matrix(double **mat, int nrows)
{
    delete_matrix<double>(mat, nrows);
}


//=============================================================================
// State spaces


// Returns state-spaces, useful for calling from python
intstate **get_state_spaces(int **ptrees, int **ages, int **sprs, 
                            int *blocklens, int ntrees, int nnodes, int ntimes)
{
    LocalTrees trees(ptrees, ages, sprs, blocklens, ntrees, nnodes);
    States states;
    
    // allocate state space
    intstate **all_states = new intstate* [ntrees];

    // iterate over local trees
    int i = 0;
    for (LocalTrees::iterator it=trees.begin(); it != trees.end(); it++) {
        LocalTree *tree = it->tree;
        get_coal_states(tree, ntimes, states);
        int nstates = states.size();
        all_states[i] = new intstate [nstates];
        
        for (int j=0; j<nstates; j++) {
            all_states[i][j][0] = states[j].node;
            all_states[i][j][1] = states[j].time;
        }
        i++;
    }

    return all_states;
}


// Deallocate state space memory
void delete_state_spaces(intstate **all_states, int ntrees)
{
    for (int i=0; i<ntrees; i++)
        delete [] all_states[i];
    delete [] all_states;
}



} // extern C

}

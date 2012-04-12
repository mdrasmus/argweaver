
// c++ includes
#include <list>
#include <vector>
#include <string.h>

// arghmm includes
#include "common.h"
#include "hmm.h"
#include "itree.h"
#include "ptree.h"
#include "seq.h"
#include "local_tree.h"

using namespace std;

using namespace spidir;
using namespace dlcoal;


namespace arghmm {


//=============================================================================
// ArgHmm model


// The model parameters and time discretization scheme
class ArgModel 
{
public:
    ArgModel(int ntimes, double *times, double *popsizes, 
             double rho, double mu) :
        ntimes(ntimes),
        times(times),
        popsizes(popsizes),
        rho(rho),
        mu(mu)
    {
        time_steps = new double [ntimes];
        for (int i=0; i<ntimes-1; i++)
            time_steps[i] = times[i+1] - times[i];
        time_steps[ntimes-1] = INFINITY;
    }
    ~ArgModel()
    {
        delete [] time_steps;
    }

    // time points
    int ntimes;
    double *times;
    double *time_steps;

    // parameters
    double *popsizes;
    double rho;
    double mu;
};




//=============================================================================
// state methods


// A state in the ArgHmm
//
// Each state represents a node and time where coalescing is allowed
class State
{
public:
    State(int node=0, int time=0) :
        node(node), time(time) {}

    int node;
    int time;
};

// A state space for a local block
typedef vector<State> States;



// This data structure provides a mapping from (node, time) tuples to
// the corresponding state index.
class NodeStateLookup
{
public:
    NodeStateLookup(const States &states, int nnodes) :
        nstates(states.size()),
        nnodes(nnodes)
    {
        const int nstates = states.size();

        // allocate lookup arrays
        node_offset = new int[nnodes];
        state_lookup = new int[nstates];


        // count number of states per node and mintime per node
        int nstates_per_node[nnodes];
        int node_mintimes[nnodes];

        // initialize arrays
        for (int i=0; i<nnodes; i++) {
            nstates_per_node[i] = 0;
            node_mintimes[i] = nstates;
        }

        for (int i=0; i<nstates; i++) {
            nstates_per_node[states[i].node]++;
            node_mintimes[states[i].node] = min(node_mintimes[states[i].node], 
                                                states[i].time);
        }

        // setup node_offsets
        int offset = 0;
        for (int i=0; i<nnodes; i++) {
            node_offset[i] = offset - node_mintimes[i];
            offset += nstates_per_node[i];
        }

        // set states
        for (int i=0; i<nstates; i++)
            state_lookup[node_offset[states[i].node] + states[i].time] = i;
    }

    ~NodeStateLookup()
    {
        // clean up lookup arrays
        delete [] node_offset;
        delete [] state_lookup;
    }

    // Returns the state index for state (node, time)
    inline int lookup(int node, int time) {
        return state_lookup[node_offset[node] + time];
    }

    int nstates;
    int nnodes;
    int *node_offset;
    int *state_lookup;
};


// A simple representation of a state, useful for passing from python
typedef int intstate[2];


// Converts integer-based states to State class
void make_states(intstate *istates, int nstates, States &states) {
    states.clear();
    for (int i=0; i<nstates; i++)
        states.push_back(State(istates[i][0], istates[i][1]));
}


// Converts state class represent to integer-based
void make_intstates(States states, intstate *istates)
{
    const int nstates = states.size();
    for (int i=0; i<nstates; i++) {
        istates[i][0] = states[i].node;
        istates[i][1] = states[i].time;
    }
}



// Retruns the possible coalescing states for a tree
//
// NOTE: Do not allow coalescing at top time
void get_coal_states(LocalTree *tree, int ntimes, States &states)
{
    states.clear();
    LocalNode *nodes = tree->nodes;
    
    // iterate over the branches of the tree
    for (int i=0; i<tree->nnodes; i++) {
        int time = nodes[i].age;
        const int parent = nodes[i].parent;
        
        if (parent == -1) {
            // no parent, allow coalescing up basal branch until ntimes-2
            for (; time<ntimes-1; time++)
                states.push_back(State(i, time));
        } else {
            // allow coalescing up branch until parent
            const int parent_age = nodes[parent].age;
            for (; time<=parent_age; time++)
                states.push_back(State(i, time));
        }
    }
}



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
    LocalTree *tree, LocalTree *last_tree, const Spr &spr,
    const States &states1, const States &states2,
    int ntimes, const double *times, int *next_states)
{
    // recomb_node in tree and last_tree
    // coal_node in last_tree
    
    const int nnodes = tree->nnodes;
    const LocalNode *nodes = tree->nodes;
    const LocalNode *last_nodes = last_tree->nodes;
    const int nstates1 = states1.size();

    // make state lookup
    NodeStateLookup state2_lookup(states2, nnodes);

    // find old node
    int old_node = last_nodes[spr.recomb_node].parent;


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

            if (node->child[0] == -1) {
                // SPR can't disrupt leaf branch
                node2 = node1;

            } else {
                int child1 = node->child[0];
                int child2 = node->child[1];
                
                if (spr.recomb_node == child1) {
                    // right child is not disrupted
                    node2 = child2;
                } else if (spr.recomb_node == child2) {
                    // left child is not disrupted
                    node2 = child1;
                } else {
                    // node is not disrupted
                    node2 = node1;
                }
            }

            // optionally walk up
            if ((spr.coal_node == node1 || 
                 (spr.coal_node == node2 && spr.coal_node != old_node)) && 
                spr.coal_time < time1)
            {
                // coal occurs under us
                // TODO: make this probabilistic (is this true?)
                node2 = nodes[node2].parent;
            }
            
            // set next state
            next_states[i] = state2_lookup.lookup(node2, time1);

        } else {
            // SPR is on same branch as new chromosome
            if (spr.recomb_time >= time1) {
                // we move with SPR subtree
                // TODO: we could probabilistically have subtree move
                // out from underneath.
                next_states[i] = state2_lookup.lookup(spr.recomb_node, time1);

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
                node2 = (other == spr.coal_node ? nodes[other].parent : other);
                next_states[i] = state2_lookup.lookup(node2, time2);
            }
        }
    }
    
}



void calc_transition_probs_switch(
    LocalTree *tree, LocalTree *last_tree, const Spr &spr,
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
    const double *times = model->times;
    const double *time_steps = model->time_steps;
    const double *popsizes = model->popsizes;
    

    // get deterministic transitions
    int determ[nstates1];
    get_deterministic_transitions(tree, last_tree, spr,
                                  states1, states2, 
                                  ntimes, times, determ);
    

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
                node3 = (c[1] == spr.recomb_node ? c[0] : c[1]);
            } else {
                node3 = node1;
            }

            // find parent of recomb_branch and node1
            int last_parent_age = last_nodes[last_parent].age;
            int parent = nodes[spr.recomb_node].parent;
            assert(parent == nodes[node3].parent);

            for (int j=0; j<nstates2; j++) {
                const int node2 = states2[j].node;
                const int time2 = states2[j].time;

                transprob[i][j] = 0.0;
                if (!((node2 == spr.recomb_node && time2 >= spr.recomb_time) ||
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
    LineageCounts lineages(ntimes);
    lineages.count(&last_tree);
    
    // setup states
    States states1, states2;
    make_states(istates1, nstates1, states1);
    make_states(istates2, nstates2, states2);
    
    double **transprob = new_matrix<double>(nstates1, nstates2);
    calc_transition_probs_switch(&tree, &last_tree, spr,
        states1, states2, &model, &lineages, transprob);
    return transprob;
}


void delete_transition_probs(double **transmat, int nstates)
{
    delete_matrix<double>(transmat, nstates);
}

} // extern C



//============================================================================
// emissions


void parsimony_ancestral_seq(LocalTree *tree, char **seqs, 
                             int nseqs, int seqlen, int pos, char *ancestral) 
{
    const int nnodes = tree->nnodes;
    LocalNode *nodes = tree->nodes;
    const int nleaves = tree->get_num_leaves();
    char sets[nnodes];
    int pchar;
    
    // clear sets
    for (int node=0; node<nnodes; node++)
        sets[node] = 0;

    // do unweighted parsimony by postorder traversal
    int postorder[nnodes];
    tree->get_postorder(postorder);
    for (int i=0; i<nnodes; i++) {
        int node = postorder[i];
        if (node < nleaves) {
            sets[node] = 1 << dna2int[(int)seqs[node][pos]];
        } else {
            char lset = sets[nodes[node].child[0]];
            char rset = sets[nodes[node].child[1]];
            char intersect = lset & rset;
            if (intersect > 0)
                sets[node] = intersect;
            else
                sets[node] = lset | rset;
        }
    }

    // traceback
    // arbitrary choose root base from set
    int root = postorder[nnodes-1];
    char rootset = sets[root];
    ancestral[root] = (rootset & 1) ? int2dna[0] :
        (rootset & 2) ? int2dna[1] :
        (rootset & 4) ? int2dna[2] : int2dna[3];
        
    // traceback with preorder traversal
    for (int i=nnodes-2; i>=0; i--) {
        int node = postorder[i];
        char s = sets[node];
            
        switch (s) {
        case 1: // just A
            ancestral[node] = int2dna[0];
            break;
        case 2: // just C
            ancestral[node] = int2dna[1];
            break;
        case 4: // just G
            ancestral[node] = int2dna[2];
            break;
        case 8: // just T
            ancestral[node] = int2dna[3];
            break;
        default:
            pchar = ancestral[nodes[node].parent];
            if (dna2int[pchar] & s) {
                // use parent char if possible
                ancestral[node] = pchar;
            } else {
                // use arbitrary char otherwise
                ancestral[node] = (s & 1) ? int2dna[0] :
                    (s & 2) ? int2dna[1] :
                    (s & 4) ? int2dna[2] : int2dna[3];
            }
        }
    }
}



void calc_emissions(const States &states, LocalTree *tree,
                    char **seqs, int nseqs, int seqlen, 
                    ArgModel *model, double **emit)
{
    const double *times = model->times;
    const double mintime = times[1];
    const double maxtime = times[model->ntimes - 1];
    const double mu = model->mu;
    const int nnodes = tree->nnodes;
    LocalNode *nodes = tree->nodes;
    
    double t1, t2, t2a, t2b, t3;
    double parent_age;
    int parent;
    int newnode = nseqs - 1;
    
    // compute ages
    double ages[nnodes];
    for (int i=0; i<nnodes; i++)
        ages[i] = times[nodes[i].age];

    // parsimony ancestral sequences
    char ancestral[nnodes];


    // base variables
    // v = new chromosome
    // x = current branch
    // p = parent of current branch
    char v, x, p;

    // iterate through positions
    for (int i=0; i<seqlen; i++) {
        v = seqs[newnode][i];

        parsimony_ancestral_seq(tree, seqs, nseqs, seqlen, i, ancestral);
        
        // iterate through states
        for (unsigned int j=0; j<states.size(); j++) {
            int node = states[j].node;
            int timei = states[j].time;
            double time = times[timei];
            double node_age = ages[node];

            x = ancestral[node];

            if (nodes[node].parent != -1) {
                parent = nodes[node].parent;
                parent_age = ages[parent];

                if (nodes[parent].parent == -1) {
                    // unwrap top branch
                    int *c = nodes[parent].child;
                    int sib = (node == c[0] ? c[1] : c[0]);
                    p = ancestral[sib];

                    // modify (x,p) length to (x,p) + (sib,p)
                    parent_age = 2 * parent_age - ages[sib];

                } else {
                    p = ancestral[parent];
                }
            } else {
                // adjust time by unwrapping branch e(v)
                parent = -1;
                parent_age = -1;
                time = 2 * time - node_age;
                p = x;
            }
            //printf(" %d %d %c %c %c\n", i, j, v, x, p);


            // ensure mintime
            if (time < mintime) 
                time = mintime;

            if (v == x && x == p) {
                // no mutation
                emit[i][j] = - mu * time;

            } else if (v != p && p == x) {
                // mutation on v
                emit[i][j] = log(.33 - .33 * exp(-mu * time));

            } else if (v == p && p != x) {
                // mutation on x
                t1 = max(parent_age - node_age, mintime);
                t2 = max(time - node_age, mintime);

                emit[i][j] = log((1 - exp(-mu *t2)) / (1 - exp(-mu * t1))
                                 * exp(-mu * (time + t2 - t1)));

            } else if (v == x && x != p) {
                // mutation on (y,p)
                t1 = max(parent_age - node_age, mintime);
                t2 = max(parent_age - time, mintime);

                emit[i][j] = log((1 - exp(-mu * t2)) / (1 - exp(-mu * t1))
                                 * exp(-mu * (time + t2 - t1)));

            } else {
                // two mutations (v,x)
                // mutation on x
                if (parent != -1) {
                    t1 = max(parent_age - node_age, mintime);
                    t2a = max(parent_age - time, mintime);
                } else {
                    t1 = max(maxtime - node_age, mintime);
                    t2a = max(maxtime - time, mintime);
                }
                t2b = max(time - node_age, mintime);
                t2 = max(t2a, t2b);
                t3 = time;

                emit[i][j] = log((1 - exp(-mu *t2)) * (1 - exp(-mu *t3))
                                 / (1 - exp(-mu * t1))
                                 * exp(-mu * (time + t2 + t3 - t1)));
            }
        }
    }
}


// C interface
extern "C" {

double **new_emissions(intstate *istates, int nstates, 
                       int *ptree, int nnodes, int *ages_index,
                       char **seqs, int nseqs, int seqlen, 
                       double *times, int ntimes,
                       double mu)
{
    States states;
    make_states(istates, nstates, states);
    LocalTree tree(ptree, nnodes, ages_index);
    ArgModel model(ntimes, times, NULL, 0.0, mu);
    
    double **emit = new_matrix<double>(seqlen, nstates);
    calc_emissions(states, &tree, seqs, nseqs, seqlen, &model, emit);

    return emit;
}


void delete_emissions(double **emit, int seqlen)
{
    delete_matrix<double>(emit, seqlen);
}

} // extern "C"

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
        if (transmat)
            delete_matrix<double>(transmat, nstates2);
        if (transmat_switch)
            delete_matrix<double>(transmat_switch, nstates1);
        if (emit)
            delete_matrix<double>(emit, blocklen);
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
        char **seqs = _seqs->seqs;
        int nseqs = _seqs->nseqs;
    
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
            for (int i=0; i<nseqs; i++)
                subseqs[i] = &seqs[i][pos];
            double **emit = new_matrix<double>(blocklen, nstates);
            calc_emissions(*states, tree, subseqs, nseqs, blocklen, 
                           model, emit);
        
            
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
                    tree, last_tree, it->spr, *last_states, *states,
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
            calc_transition_probs_switch(tree, last_tree, it->spr,
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
        }
    }
}


/*

def sample_recombinations_thread(model, thread, use_times=True):
    """Samples new recombination for a thread"""
    
    r = 0
    
    # assumes that recomb_pos starts with -1 and ends with arg.end
    arg_recomb = model.recomb_pos
    time_lookup = util.list2lookup(model.times)
    minlen = model.time_steps[0]
    
    tree = model.arg.get_marginal_tree(-.5)
    treelen = get_treelen(tree, model.times)
    new_node = model.new_name
    transmat = None
    nstates = 0
    selftrans = None

    next_recomb = -1
    
    for pos, state in enumerate(thread):
        node, node_time = state
        timei = time_lookup[node_time]
        
        # update local tree if needed
        while r < len(arg_recomb) and arg_recomb[r] < pos:
            r += 1
            tree = model.arg.get_marginal_tree(pos-.5)
            treelen = get_treelen(tree, model.times)
            nlineages = get_nlineages_recomb_coal(tree, model.times)
            nbranches, nrecombs, ncoals = nlineages

            if transmat is not None:
                delete_transition_probs(transmat, nstates)
            transmat = calc_transition_probs_c(
                tree, model.states[pos], nlineages,
                model.times, model.time_steps, model.popsizes, model.rho)
            nstates = len(model.states[pos])
            statei = model.states[pos].index((node, timei))
            selftrans = transmat[statei][statei]
            

        if pos == 0 or arg_recomb[r-1] == pos - 1:
            # previous arg recomb is right behind us, sample no recomb
            next_recomb = -1
            continue

        # get information about pos-1
        # since their no recomb in G_{n-1}, last_tree == tree
        last_state = thread[pos-1]
        last_node, last_time = last_state
        last_timei = time_lookup[last_time]
        last_tree = tree
        last_treelen = treelen
        
        if state == last_state:
            if pos > next_recomb:
                # sample the next recomb pos
                last_treelen2 = get_treelen_branch(
                    last_tree, model.times, last_node, last_time)
                rate = max(1.0 - exp(-model.rho * (last_treelen2 - last_treelen)
                                     - selftrans), model.rho)
                next_recomb = pos + int(random.expovariate(rate))
                
            if pos < next_recomb:
                continue


        next_recomb = -1
        last_treelen2 = get_treelen_branch(
            last_tree, model.times, last_node, last_time)
        statei = model.states[pos].index((node, timei))
        selftrans = transmat[statei][statei]

        # there must be a recombination
        # either because state changed or we choose to recombine
        if node == last_node:
            if timei == last_timei:
                # y = v, k in [0, min(timei, last_timei)]
                # y = node, k in Sr(node)
                # if node.parent.age == model.times[timei],
                #   y = sis(last_tree, node.name), k in Sr(y)
                node_timei = time_lookup[tree[node].age]
                recombs = [(new_node, k) for k in
                           range(0, min(timei, last_timei)+1)] + \
                          [(node, k) for k in
                           range(node_timei, min(timei, last_timei)+1)]
            else:
                # y = v, k in [0, min(timei, last_timei)]
                # y = node, k in Sr(node)
                node_timei = time_lookup[tree[node].age]
                recombs = [(new_node, k) for k in
                           range(0, min(timei, last_timei)+1)] + \
                          [(node, k) for k in
                           range(node_timei, min(timei, last_timei)+1)]
        else:
            # y = v, k in [0, min(timei, last_timei)]
            recombs = [(new_node, k)
                       for k in range(0, min(timei, last_timei)+1)]

        if len(recombs) == 0:
            print ((last_node, last_timei), (node, timei))
            raise Exception("recomb not sampled!")

        C = calc_C(model.time_steps, nbranches, model.popsizes)
        j = timei
        probs = []
        for recomb in recombs:
            k = recomb[1]
            probs.append((nbranches[k] + 1) * model.time_steps[k] /
                         (ncoals[j] * (nrecombs[k] + 1.0) * last_treelen2) *
                         (1.0 - exp(-model.time_steps[j] * nbranches[j] /
                                    (2.0 * model.popsizes[j-1]))) *
                         (1.0 - exp(-model.rho * last_treelen2)) *
                         exp(-C[k] + C[k]))
        recomb_node, recomb_time = recombs[stats.sample(probs)]
        
        if use_times:
            recomb_time = model.times[recomb_time]
        yield (pos, recomb_node, recomb_time)

*/



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
    for (LocalTrees::iterator it=trees.begin(); it != trees.end(); it++) {
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

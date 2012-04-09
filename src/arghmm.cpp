
// c++ includes
#include <list>
#include <vector>

// arghmm includes
#include "common.h"
#include "hmm.h"
#include "itree.h"
#include "ptree.h"
#include "seq.h"

using namespace std;

using namespace spidir;
using namespace dlcoal;


namespace arghmm {

extern "C" {

typedef int intstate[2];



//=============================================================================
// Local trees


// A block within a sequence alignment
class Block 
{
public:
    Block(int start, int end) : 
        start(start), end(end) {}
    int start;
    int end;
};


// A Subtree Pruning and Regrafting operation
class Spr
{
public:
    Spr() {}
    Spr(int recomb_node, int recomb_time, 
        int coal_node, int coal_time) :
        recomb_node(recomb_node), recomb_time(recomb_time),
        coal_node(coal_node), coal_time(coal_time) {}

    int recomb_node;
    int recomb_time;
    int coal_node;
    int coal_time;
};



// A node in a local tree
class LocalNode 
{
public:
    LocalNode() {}
    LocalNode(int parent, int left_child, int right_child, int age=-1) :
        parent(parent), age(age)
    {
        child[0] = left_child;
        child[1] = right_child;
    }

    inline bool is_leaf()
    {
        return child[0] == -1;
    }

    int parent;
    int child[2];
    int age;
};


// A local tree in a set of local trees
class LocalTree
{
public:
    LocalTree() :
        nnodes(0),
        root(-1),
        nodes(NULL) 
    {}

    LocalTree(int *ptree, int nnodes, int *ages=NULL) :
        nodes(NULL) 
    {
        set_ptree(ptree, nnodes, ages);
    }

    ~LocalTree() {
        if (nodes)
            delete [] nodes;
    }

    // initialize a local tree by on a parent array
    void set_ptree(int *ptree, int _nnodes, int *ages=NULL) 
    {
        nnodes = _nnodes;
        if (nodes)
            delete [] nodes;
        nodes = new LocalNode [nnodes];

        // initialize nodes
        for (int i=0; i<nnodes; i++) {
            nodes[i].parent = ptree[i];
            nodes[i].child[0] = -1;
            nodes[i].child[1] = -1;
        }

        if (ages)
            for (int i=0; i<nnodes; i++)
                nodes[i].age = ages[i];

    
        // populate children
        for (int i=0; i<nnodes; i++) {
            const int parent = ptree[i];
        
            if (parent != -1) {
                int *child = nodes[parent].child;
                if (child[0] == -1)
                    child[0] = i;
                else
                    child[1] = i;
            } else {
                root = i;
            }
        }
    }


    void get_postorder(int *order)
    {
        char visit[nnodes];
        int i;
        for (i=0; i<nnodes; i++)
            visit[i] = 0;

        // list leaves first
        for (i=0; i<nnodes; i++) {
            if (!nodes[i].is_leaf())
                break;
            order[i] = i;
        }
        
        // add the remaining nodes
        int end = i;
        for (i=0; i<nnodes; i++) {
            int parent = nodes[order[i]].parent;
            if (parent != -1) {
                visit[parent]++;
                
                // add parent to queue if both children have been seen
                if (visit[parent] == 2)
                    order[end++] = parent;
            }
        }
    }
    

    int nnodes;
    int root;
    LocalNode *nodes;
};



bool assert_tree_postorder(LocalTree *tree, int *order)
{
    if (tree->root != order[tree->nnodes-1])
        return false;

    char seen[tree->nnodes];
    for (int i=0; i<tree->nnodes; i++)
        seen[i] = 0;

    for (int i=0; i<tree->nnodes; i++) {
        int node = order[i];
        seen[node] = 1;
        if (!tree->nodes[node].is_leaf()) {
            if (! seen[tree->nodes[node].child[0]] ||
                ! seen[tree->nodes[node].child[1]])
                return false;
        }
    }
    
    return true;
}




class LocalTreeSpr
{
public:
    LocalTreeSpr(int start, int end, LocalTree *tree, int *ispr) :
        block(start, end),
        tree(tree),
        spr(ispr[0], ispr[1], ispr[2], ispr[3])
    {}

    LocalTreeSpr(int start, int end, LocalTree *tree, Spr spr) :
        block(start, end),
        tree(tree),
        spr(spr)
    {}

    Block block;
    LocalTree *tree;
    Spr spr;
};


class LocalTrees 
{
public:
    LocalTrees(int **ptrees, int**ages, int **isprs, int *blocklens,
               int ntrees, int nnodes, int start=0) : 
        start_coord(start),
        nnodes(nnodes)
    {
        // copy data
        int pos = start;
        for (int i=0; i<ntrees; i++) {
            end_coord = pos + blocklens[i];
            trees.push_back(
                LocalTreeSpr(pos, end_coord,
                             new LocalTree(ptrees[i], nnodes, ages[i]),
                             isprs[i]));
            pos = end_coord;
        }

    }
    ~LocalTrees() {}

    typedef list<LocalTreeSpr>::iterator iterator;

    iterator begin()
    {
        return trees.begin();
    }

    iterator end()
    {
        return trees.end();
    }


    int start_coord;
    int end_coord;
    int nnodes;
    list<LocalTreeSpr> trees;
};



void count_lineages(LocalTree *tree, int ntimes,
                    int *nbranches, int *nrecombs, int *ncoals)
{
    const LocalNode *nodes = tree->nodes;

    // initialize counts
    for (int i=0; i<ntimes; i++) {
        nbranches[i] = 0;
        nrecombs[i] = 0;
        ncoals[i] = 0;
    }

    for (int i=0; i<tree->nnodes; i++) {
        const int parent = nodes[i].parent;
        const int parent_age = ((parent == -1) ? ntimes - 2 : 
                                nodes[parent].age);
        
        // add counts for every segment along branch
        for (int j=nodes[i].age; j<parent_age; j++) {
            nbranches[j]++;
            nrecombs[j]++;
            ncoals[j]++;
        }

        // recomb and coal is also allowed at the top of a branch
        nrecombs[parent_age]++;
        ncoals[parent_age]++;
        if (parent == -1)
            nbranches[parent_age]++;
    }
    
    nbranches[ntimes - 1] = 1;
}



class LineageCounts
{
public:
    LineageCounts(int ntimes) :
        ntimes(ntimes)
    {
        nbranches = new int [ntimes];
        nrecombs = new int [ntimes];
        ncoals = new int [ntimes];
    }

    ~LineageCounts()
    {
        delete [] nbranches;
        delete [] nrecombs;
        delete [] ncoals;
    }

    inline void count(LocalTree *tree) {
        count_lineages(tree, ntimes, nbranches, nrecombs, ncoals);
    }

    int ntimes;
    int *nbranches;
    int *nrecombs;
    int *ncoals;
};


// Calculate tree length according to ArgHmm rules
double get_treelen(const LocalTree *tree, const double *times, int ntimes)
{
    double treelen = 0.0;
    const LocalNode *nodes = tree->nodes;
    
    for (int i=0; i<tree->nnodes; i++) {
        int parent = nodes[i].parent;
        int age = nodes[i].age;
        if (parent == -1) {
            // add basal stub
            treelen += times[age+1] - times[age];
        } else {
            treelen += times[nodes[parent].age] - times[age];
        }
    }
    
    return treelen;
}


//=============================================================================
// ArgHmm model


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


class Sequences
{
public:
    Sequences(char **seqs, int nseqs, int seqlen) :
        seqs(seqs), nseqs(nseqs), seqlen(seqlen)
    {}
    
    char **seqs;
    int nseqs;
    int seqlen;
};



//=============================================================================
// state methods


class State
{
public:
    State(int node=0, int time=0) :
        node(node), time(time) {}

    int node;
    int time;
};


typedef vector<State> States;


//
//  This data structure provides a mapping from (node, time) tuples
//  to the corresponding state index.
//
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

    inline int lookup(int node, int time) {
        return state_lookup[node_offset[node] + time];
    }

    int nstates;
    int nnodes;
    int *node_offset;
    int *state_lookup;
};


// convert integer-based states to State class
void make_states(intstate *istates, int nstates, States &states) {
    states.clear();
    for (int i=0; i<nstates; i++)
        states.push_back(State(istates[i][0], istates[i][1]));
}


void make_intstates(States states, intstate *istates)
{
    const int nstates = states.size();
    for (int i=0; i<nstates; i++) {
        istates[i][0] = states[i].node;
        istates[i][1] = states[i].time;
    }
}


//
// get the possible coalescing states for a tree
// NOTE: Do not allow coalescing at top time
//
void get_coal_states(LocalTree *tree, int ntimes, States &states)
{
    states.clear();
    LocalNode *nodes = tree->nodes;
    
    for (int i=0; i<tree->nnodes; i++) {
        int time = nodes[i].age;
        const int parent = nodes[i].parent;
        
        if (parent == -1) {
            for (; time<ntimes-1; time++)
                states.push_back(State(i, time));
        } else {
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
    double D[ntimes];
    double S[ntimes * ntimes];
    
    C[0] = 0.0;
    B[0] = nbranches[0] * time_steps[0] / (nrecombs[0] + 1.0);
    D[0] = (1.0 - exp(-time_steps[0] * nbranches[0]
                          / (2.0 * popsizes[0]))) / ncoals[0];

    for (int b=1; b<ntimes; b++) {
        const int l = b - 1;
        C[b] = C[l] + time_steps[l] * nbranches[l] / (2.0 * popsizes[l]);
        B[b] = B[b-1] + nbranches[b] * time_steps[b] / 
            (nrecombs[b] + 1.0) * exp(C[b]);
        D[b] = (1.0 - exp(-time_steps[b] * nbranches[b]
                          / (2.0 * popsizes[b]))) / ncoals[b];
    }


    for (int b=0; b<ntimes; b++) {
        int a;
        const double ec = exp(-C[b]);
        for (a=0; a<b; a++)
            S[matind(ntimes, a, b)] = ec * B[a];
        const double sab = ec * B[b];
        for (; a<ntimes; a++)
            S[matind(ntimes, a, b)] = sab;
    }
    
    // f =\frac{[1 - \exp(- \rho (|T^{n-1}_{i-1}| + s_a))] 
    //       [1 - \exp(- s'_b k_b / (2N))]}
    //      {\exp(-\rho |T^{n-1}_{i-1}|) (|T^{n-1}_{i-1}| + s_a) k^C_b}
    for (int i=0; i<nstates; i++) {
        const int node1 = states[i].node;
        const int a = states[i].time;
        const int c = nodes[node1].age;
        assert(a < ntimes);

        double treelen2 = treelen + times[a];
        if (node1 == root) {
            treelen2 += times[a] - root_age;
            treelen2 += time_steps[a]; // add basal branch
        } else {
            treelen2 += time_steps[root_age_index]; // add basal branch
        }

        double const F = (1.0 - exp(-rho * treelen2)) /
            (exp(-rho * treelen) * treelen2);
        
        for (int j=0; j<nstates; j++) {
            const int node2 = states[j].node;
            const int b = states[j].time;
            assert(b < ntimes);
            
            const double f = F * D[b];
            
            if (node1 != node2)
                transprob[i][j] = f * S[matind(ntimes, a, b)];
            else {
                transprob[i][j] = f * (2*S[matind(ntimes, a, b)] - 
                                       S[matind(ntimes, c, b)]);
                if (a == b)
                    transprob[i][j] += exp(-rho * (treelen2 - treelen));
            }
        }

        // compute diagonal
        double sum = 0.0;
        for (int j=0; j<nstates; j++)
            sum += transprob[i][j];
        //transprob[i][i] = 1.0 - sum;

        // convert to log scale
        for (int j=0; j<nstates; j++)
            transprob[i][j] = log(transprob[i][j] / sum);
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




double **new_transition_probs(int nnodes, int *ptree, 
                              int *ages, double treelen,
                              intstate *istates, int nstates,
                              int ntimes, double *times, double *time_steps,
                              int *nbranches, int *nrecombs, int *ncoals, 
                              double *popsizes, double rho)
{

    // setup model
    ArgModel model(ntimes, times, popsizes, rho, 0.0);
    
    // setup local tree
    LocalTree tree(ptree, nnodes, ages);
    LineageCounts lineages(ntimes);
    lineages.count(&tree);
    
    /*
    for (int i=0; i<ntimes; i++) {
        //printf("b[%d] %d %d\n", i, nbranches[i], lineages.nbranches[i]);
        //printf("r[%d] %d %d\n", i, nrecombs[i], lineages.nrecombs[i]);
        //printf("c[%d] %d %d\n", i, ncoals[i], lineages.ncoals[i]);
        assert(nbranches[i] == lineages.nbranches[i]);
        assert(nrecombs[i] == lineages.nrecombs[i]);
        assert(ncoals[i] == lineages.ncoals[i]);
    }*/

    // setup states
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
    
    /*
    for (int i=0; i<ntimes; i++) {
        //printf("b[%d] %d %d\n", i, nbranches[i], lineages.nbranches[i]);
        //printf("r[%d] %d %d\n", i, nrecombs[i], lineages.nrecombs[i]);
        //printf("c[%d] %d %d\n", i, ncoals[i], lineages.ncoals[i]);
        assert(nbranches[i] == lineages.nbranches[i]);
        assert(nrecombs[i] == lineages.nrecombs[i]);
        assert(ncoals[i] == lineages.ncoals[i]);
    }
    */

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
                                   (2.0 * popsizes[b]))) / ncoals[b] * 
                        exp(-sum)); 
    }
}



//============================================================================
// emissions


void parsimony_ancestral_seq(LocalTree *tree, char **seqs, 
                             int nseqs, int seqlen, int pos, char *ancestral) 
{
    const int nnodes = tree->nnodes;
    LocalNode *nodes = tree->nodes;
    const int nleaves = (nnodes + 1) / 2;
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
        
        //printf("%d ", i);
        //for (int j=0; j<nnodes; j++)
        //    putc(ancestral[j], stdout);
        //printf("\n");

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
    ~ArgHmmMatrices() {}

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

    ArgModel *model;
    vector<ArgHmmMatrices> matrices;
};


//=============================================================================
// HMM algorithms



double **arghmm_forward_alg(
    int **ptrees, int **ages, int **sprs, int *blocklens,
    int ntrees, int nnodes, double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, double **fw=NULL)
{
    
    // allocate lineage counts
    LineageCounts lineages(ntimes);
    
    // setup model
    ArgModel model(ntimes, times, popsizes, rho, mu);
    
    // setup local trees
    LocalTrees trees(ptrees, ages, sprs, blocklens, ntrees, nnodes);
    LocalTree *last_tree = NULL;
    
    // state spaces
    States states1;
    States states2;
    States *states = &states1;
    States *last_states = NULL;
    
    // temporary subsequence of the sequence
    char *subseqs[nseqs];
    
    // allocate the forward table if necessary
    if (fw == NULL) {
        fw = new double* [seqlen];
        for (int i=0; i<seqlen; i++)
            fw[i] = NULL;
    }
    

    // iterate over local trees to find largest blocklen and nstates
    int max_nstates = 0;
    int max_blocklen = 0;
    for (LocalTrees::iterator it=trees.begin(); it != trees.end(); it++) {
        int blocklen = it->block.end - it->block.start;
        get_coal_states(it->tree, ntimes, *states);
        int nstates = states->size();
        max_nstates = max(max_nstates, nstates);
        max_blocklen = max(max_blocklen, blocklen);
    }


    // HMM matrices
    double **emit = new_matrix<double>(max_blocklen, max_nstates);
    double **transmat = new_matrix<double>(max_nstates, max_nstates);
    double **transmat_switch = new_matrix<double>(max_nstates, max_nstates);
    
    
    // iterate over local trees
    for (LocalTrees::iterator it=trees.begin(); it != trees.end(); it++) {
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
        for (int i=0; i<nseqs; i++)
            subseqs[i] = &seqs[i][pos];
        calc_emissions(*states, tree, subseqs, nseqs, blocklen, &model, emit);
        
        
        // use switch matrix for first column of forward table
        // if we have a previous state space (i.e. not first block)
        if (!last_states) {
            // calculate prior of first state
            lineages.count(tree);
            calc_state_priors(*states, &lineages, &model, fw[0]);  

        } else {
            double *col1 = fw[pos-1];
            double *col2 = fw[pos];
            int nstates1 = last_states->size();
            int nstates2 = states->size();
            
            // calculate transmat_switch
            lineages.count(last_tree);
            calc_transition_probs_switch(tree, last_tree, it->spr,
                                         *last_states, *states,
                                         &model, &lineages, transmat_switch);

            // perform one column of forward algorithm
            double tmp[nstates1];
            for (int k=0; k<nstates2; k++) {
                double e = emit[0][k];
                for (int j=0; j<nstates1; j++)
                    tmp[j] = col1[j] + transmat_switch[j][k] + e;
                col2[k] = logsum(tmp, nstates1);
            }
            
            // update lineages to current tree
            lineages.count(tree);
        }
        
        // calculate transmat and use it for rest of block
        calc_transition_probs(tree, &model, *states, &lineages, transmat);
        
        double **subfw = &fw[pos];
        forward_alg(blocklen, nstates, transmat, emit, subfw);
        
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



double **arghmm_posterior_sample(
    int **ptrees, int **ages, int **sprs, int *blocklens,
    int ntrees, int nnodes, double *times, int ntimes,
    double *popsizes, double rho, double mu,
    char **seqs, int nseqs, int seqlen, double **fw=NULL)
{
    
    // allocate lineage counts
    LineageCounts lineages(ntimes);
    
    // setup model
    ArgModel model(ntimes, times, popsizes, rho, mu);
    
    // setup local trees
    LocalTrees trees(ptrees, ages, sprs, blocklens, ntrees, nnodes);
    LocalTree *last_tree = NULL;
    
    // state spaces
    States states1;
    States states2;
    States *states = &states1;
    States *last_states = NULL;
    
    // temporary subsequence of the sequence
    char *subseqs[nseqs];
    
    // allocate the forward table if necessary
    if (fw == NULL) {
        fw = new double* [seqlen];
        for (int i=0; i<seqlen; i++)
            fw[i] = NULL;
    }
    
    
    // iterate over local trees
    for (LocalTrees::iterator it=trees.begin(); it != trees.end(); it++) {
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
        for (int i=0; i<nseqs; i++)
            subseqs[i] = &seqs[i][pos];
        double **emit = new_matrix<double>(blocklen, nstates);
        calc_emissions(*states, tree, subseqs, nseqs, blocklen, &model, emit);
        
        
        // use switch matrix for first column of forward table
        // if we have a previous state space (i.e. not first block)
        if (!last_states) {
            // calculate prior of first state
            lineages.count(tree);
            calc_state_priors(*states, &lineages, &model, fw[0]);  

        } else {
            double *col1 = fw[pos-1];
            double *col2 = fw[pos];
            int nstates1 = last_states->size();
            int nstates2 = states->size();
            
            // calculate transmat_switch
            lineages.count(last_tree);
            double **transmat_switch = new_matrix<double>(nstates1, nstates2);
            
            calc_transition_probs_switch(tree, last_tree, it->spr,
                                         *last_states, *states,
                                         &model, &lineages, transmat_switch);

            // perform one column of forward algorithm
            double tmp[nstates1];
            for (int k=0; k<nstates2; k++) {
                double e = emit[0][k];
                for (int j=0; j<nstates1; j++)
                    tmp[j] = col1[j] + transmat_switch[j][k] + e;
                col2[k] = logsum(tmp, nstates1);
            }

            delete_matrix<double>(transmat_switch, nstates1);

            // update lineages to current tree
            lineages.count(tree);
        }
        
        // calculate transmat and use it for rest of block
        double **transmat = new_matrix<double>(nstates, nstates);
        calc_transition_probs(tree, &model, *states, &lineages, transmat);
        
        double **subfw = &fw[pos];
        forward_alg(blocklen, nstates, transmat, emit, subfw);
        

        // clean up
        delete_matrix<double>(transmat, nstates);
        delete_matrix<double>(emit, blocklen);

        // update pointers
        last_tree = tree;
        last_states = states;
        states = ((states == &states1) ? &states2 : &states1);
    }

    return fw;
}




void delete_double_matrix(double **mat, int nrows)
{
    delete_matrix<double>(mat, nrows);
}



intstate **get_state_spaces(int **ptrees, int **ages, int **sprs, 
                            int *blocklens, int ntrees, int nnodes, int ntimes)
{
    LocalTrees trees(ptrees, ages, sprs, blocklens, ntrees, nnodes);
    States states;
    
    // allocate state space
    intstate **all_states = new intstate* [ntrees];

    // iterate over local trees
    int i = 0;
    for (LocalTrees::iterator it=trees.begin(); it != trees.end(); it++)
    {
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


void delete_state_spaces(intstate **all_states, int ntrees)
{
    for (int i=0; i<ntrees; i++) {
        delete [] all_states[i];
    }
}



} // extern C

}

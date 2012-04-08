
// c++ includes
#include <list>
#include <vector>

// arghmm includes
#include "common.h"
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
            }
        }
    }
    

    int nnodes;
    LocalNode *nodes;
};


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
        const int parent = tree->nodes[i].parent;
        const int parent_age = ((parent == -1) ? ntimes - 2 : 
                                nodes[parent].age);
        
        for (int j=nodes[i].age; j<parent_age; j++) {
            nbranches[j]++;
            nrecombs[j]++;
            ncoals[j]++;
        }

        if (parent_age > nodes[i].age) {
            nrecombs[parent_age]++;
            ncoals[parent_age]++;
        }            
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
    NodeStateLookup(intstate *states, int nstates, int nnodes) :
        nstates(nstates),
        nnodes(nnodes)
    {
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
            nstates_per_node[states[i][0]]++;
            node_mintimes[states[i][0]] = min(node_mintimes[states[i][0]], 
                                              states[i][1]);
        }

        // setup node_offsets
        int offset = 0;
        for (int i=0; i<nnodes; i++) {
            node_offset[i] = offset - node_mintimes[i];
            offset += nstates_per_node[i];
        }

        // set states
        for (int i=0; i<nstates; i++)
            state_lookup[node_offset[states[i][0]] + states[i][1]] = i;
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



void calc_transition_probs(int nnodes, 
                           int *ptree, int *ages_index, double treelen,
                           intstate *states, int nstates,
                           int ntimes, double *times, double *time_steps,
                           int *nbranches, int *nrecombs, int *ncoals, 
                           double *popsizes, double rho,
                           double **transprob)
{
    // find root node
    int root = 0;
    for (int i=0; i<nnodes; i++) {
        if (ptree[i] == -1) {
            root = i;
            break;
        }
    }
    const int root_age_index = ages_index[root];
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
        const int node1 = states[i][0];
        const int a = states[i][1];
        const int c = ages_index[node1];
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
            const int node2 = states[j][0];
            const int b = states[j][1];
            assert(b < ntimes);
            
            const double f = F * D[b];
            if (node1 != node2)
                transprob[i][j] = f * S[matind(ntimes, a, b)];
            else if (a != b) 
                transprob[i][j] = f * (2*S[matind(ntimes, a, b)] - 
                                       S[matind(ntimes, c, b)]);
            else
                transprob[i][j] = 0.0;
        }

        // compute diagonal
        double sum = 0.0;
        for (int j=0; j<nstates; j++)
            sum += transprob[i][j];
        transprob[i][i] = 1.0 - sum;

        // convert to log scale
        for (int j=0; j<nstates; j++)
            transprob[i][j] = log(transprob[i][j]);
    }
}



void get_deterministic_transitions(
    intnode *itree, intnode *last_itree, int nnodes, 
    int *last_ages_index,
    int recomb_name, int recomb_time, int coal_name, int coal_time,
    intstate *states1, int nstates1, intstate *states2, int nstates2,
    int ntimes, double *times, int *next_states)
{

    // recomb_node in tree and last_tree
    // coal_node in last_tree


    // make state lookup
    NodeStateLookup state2_lookup(states2, nstates2, nnodes);

    // find old node
    int old_node = last_itree[recomb_name].parent;


    for (int i=0; i<nstates1; i++) {
        const int node1 = states1[i][0];
        const int time1 = states1[i][1];

        if (node1 == coal_name && time1 == coal_time) {
            // not a deterministic case
            next_states[i] = -1;
        
        } else if (node1 != recomb_name) {
            // SPR only removes a subset of descendents, if any
            // trace up from remaining leaf to find correct new state

            int node2;
            intnode *node = &last_itree[node1];

            if (node->child[0] == -1) {
                // SPR can't disrupt leaf branch
                node2 = node1;

            } else {
                int child1 = node->child[0];
                int child2 = node->child[1];
                
                if (recomb_name == child1) {
                    // right child is not disrupted
                    node2 = child2;
                } else if (recomb_name == child2) {
                    // left child is not disrupted
                    node2 = child1;
                } else {
                    // node is not disrupted
                    node2 = node1;
                }
            }

            // optionally walk up
            if ((coal_name == node1 || 
                 (coal_name == node2 && coal_name != old_node)) && 
                coal_time < time1)
            {
                // coal occurs under us
                // TODO: make this probabilistic (is this true?)
                node2 = itree[node2].parent;
            }
            
            // set next state
            next_states[i] = state2_lookup.lookup(node2, time1);

        } else {
            // SPR is on same branch as new chromosome
            if (recomb_time >= time1) {
                // we move with SPR subtree
                // TODO: we could probabilistically have subtree move
                // out from underneath.
                next_states[i] = state2_lookup.lookup(recomb_name, time1);

            } else {
                // SPR should not be able to coal back onto same branch
                // this would be a self cycle
                assert(coal_name != node1);
                
                // SPR subtree moves out from underneath us
                // therefore therefore the new chromosome coalesces with
                // the branch above the subtree

                // search up for parent
                int parent = last_itree[recomb_name].parent;
                int time2 = last_ages_index[parent];
                int node2;

                // find other child
                int *c = last_itree[parent].child;
                int other = (c[1] == recomb_name ? c[0] : c[1]);

                // find new state in tree
                node2 = (other == coal_name ? itree[other].parent : other);
                next_states[i] = state2_lookup.lookup(node2, time2);
            }
        }
    }
    
}



void calc_transition_probs_switch(
    int *ptree, int *last_ptree, int nnodes, 
    int recomb_name, int recomb_time, int coal_name, int coal_time,
    int *ages_index, int *last_ages_index,
    double treelen, double last_treelen,
    intstate *states1, int nstates1,
    intstate *states2, int nstates2,
    
    int ntimes, double *times, double *time_steps,
    int *nbranches, int *nrecombs, int *ncoals, 
    double *popsizes, double rho,
    double **transprob)
{
    
    // create inttree data structures
    intnode *last_itree = make_itree(nnodes, last_ptree);
    intnode *itree = make_itree(nnodes, ptree);

    // find root node
    int last_root;
    for (int i=0; i<nnodes; i++) {
        if (itree[i].parent == -1) {
            last_root = i;
            break;
        }
    }
    

    // get deterministic transitions
    int determ[nstates1];
    get_deterministic_transitions(itree, last_itree, nnodes,
                                  last_ages_index,
                                  recomb_name, recomb_time,
                                  coal_name, coal_time,
                                  states1, nstates1, states2, nstates2,
                                  ntimes, times, determ);
    

    for (int i=0; i<nstates1; i++) {
        const int node1 = states1[i][0];
        const int time1 = states1[i][1];

        if (node1 != coal_name || time1 != coal_time) {
            // deterministic transition case
            assert (determ[i] != -1);
            for (int j=0; j<nstates2; j++)
                transprob[i][j] = -INFINITY;
            transprob[i][determ[i]] = 0.0;

        } else {
            // probabilistic transition case

            // determine if node1 is still here or not
            int node3;
            int last_parent = last_itree[recomb_name].parent;
            if (last_parent == node1) {
                // recomb breaks node1 branch, we need to use the other child
                int *c = last_itree[last_parent].child;
                node3 = (c[1] == recomb_name ? c[0] : c[1]);
            } else {
                node3 = node1;
            }

            // find parent of recomb_branch and node1
            int last_parent_age = last_ages_index[last_parent];
            int parent = itree[recomb_name].parent;
            assert(parent == itree[node3].parent);

            // treelen of T^n_{i-1}
            double blen = times[time1];
            double treelen2 = treelen + blen;
            if (node1 == last_root) {
                treelen2 += blen - times[last_ages_index[last_root]];
                treelen2 += time_steps[time1];
            } else {
                treelen2 += time_steps[last_ages_index[last_root]];
            }


            for (int j=0; j<nstates2; j++) {
                int node2 = states2[j][0];
                int time2 = states2[j][1];

                transprob[i][j] = 0.0;
                if (! ((node2 == recomb_name && time2 >= recomb_time) ||
                       (node2 == node3 && time2 == time1) ||
                       (node2 == parent && time2 == time1)))
                    // not a probabilistic transition
                    continue;

                // get lineage counts
                // remove recombination branch and add new branch
                int kbn = nbranches[time2];
                int kcn = ncoals[time2] + 1;
                if (time2 < ages_index[parent]) {
                    kbn -= 1;
                    kcn -= 1;
                }
                if (time2 < time1)
                    kbn += 1;

                double sum = 0.0;
                for (int m=recomb_time; m<time2; m++) {
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

    free_itree(last_itree);
    free_itree(itree);
}




double **new_transition_probs(int nnodes, int *ptree, 
                              int *ages_index, double treelen,
                              intstate *states, int nstates,
                              int ntimes, double *times, double *time_steps,
                              int *nbranches, int *nrecombs, int *ncoals, 
                              double *popsizes, double rho)
{
    double **transprob = new_matrix<double>(nstates, nstates);
    calc_transition_probs(nnodes, ptree, ages_index, treelen, 
                          states, nstates,
                          ntimes, times, time_steps,
                          nbranches, nrecombs, ncoals, 
                          popsizes, rho, transprob);
    return transprob;
}


double **new_transition_probs_switch(
    int *ptree, int *last_ptree, int nnodes, 
    int recomb_name, int recomb_time, int coal_name, int coal_time,
    int *ages_index, int *last_ages_index,
    double treelen, double last_treelen,
    intstate *states1, int nstates1,
    intstate *states2, int nstates2,
    
    int ntimes, double *times, double *time_steps,
    int *nbranches, int *nrecombs, int *ncoals, 
    double *popsizes, double rho)
{
    double **transprob = new_matrix<double>(nstates1, nstates2);
    calc_transition_probs_switch(
        ptree, last_ptree, nnodes, 
        recomb_name, recomb_time, coal_name, coal_time,
        ages_index, last_ages_index,
        treelen, last_treelen,
        states1, nstates1, states2, nstates2,
        ntimes, times, time_steps,
        nbranches, nrecombs, ncoals, 
        popsizes, rho, transprob);
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
        for (int m=0; m<b; b++)
            sum += time_steps[m] * nbranches[m] / (2.0 * popsizes[m]);
        
        priors[i] = log((1.0 - exp(- time_steps[b] * nbranches[b] /
                                   (2.0 * popsizes[b]))) / ncoals[b] * 
                        exp(-sum)); 
    }
}



//============================================================================
// emissions

void parsimony_ancestral_seq_all(intnode *itree, int nnodes, char **seqs, 
                                 int nseqs, int seqlen, char **ancestral) 
{
    char *sets = new char [nnodes];
    const int nleaves = (nnodes + 1) / 2;
    int pchar;
    
    for (int i=0; i<seqlen; i++) {
        // clear sets
        for (int node=0; node<nnodes; node++)
            sets[node] = 0;

        // do unweighted parsimony by postorder traversal
        for (int node=0; node<nnodes; node++) {
            if (node < nleaves) {
                sets[node] = 1 << dna2int[(int)seqs[node][i]];
            } else {
                char lset = sets[itree[node].child[0]];
                char rset = sets[itree[node].child[1]];
                char intersect = lset & rset;

                assert(lset != 0 && rset != 0);

                if (intersect > 0)
                    sets[node] = intersect;
                else
                    sets[node] = lset | rset;
            }
        }

        // traceback
        // arbitrary choose root base from set
        char rootset = sets[nnodes-1];
        ancestral[nnodes-1][i] = (rootset & 1) ? int2dna[0] :
            (rootset & 2) ? int2dna[1] :
            (rootset & 4) ? int2dna[2] : int2dna[3];
        
        // traceback with preorder traversal
        for (int node=nnodes-2; node>-1; node--) {
            char s = sets[node];
            
            switch (s) {
            case 1: // just A
                ancestral[node][i] = int2dna[0];
                break;
            case 2: // just C
                ancestral[node][i] = int2dna[1];
                break;
            case 4: // just G
                ancestral[node][i] = int2dna[2];
                break;
            case 8: // just T
                ancestral[node][i] = int2dna[3];
                break;
            default:
                pchar = ancestral[itree[node].parent][i];
                if (dna2int[pchar] & s) {
                    // use parent char if possible
                    ancestral[node][i] = pchar;
                } else {
                    // use arbitrary char otherwise
                    ancestral[node][i] = (s & 1) ? int2dna[0] :
                        (s & 2) ? int2dna[1] :
                        (s & 4) ? int2dna[2] : int2dna[3];
                }
            }
        }
    }

    delete [] sets;
}


void get_posterior_order(intnode *itree, int nnodes)
{
    
}

void parsimony_ancestral_seq(intnode *itree, int nnodes, char **seqs, 
                             int nseqs, int seqlen, int pos, char *ancestral) 
{
    char sets[nnodes];
    const int nleaves = (nnodes + 1) / 2;
    int pchar;
    
        // clear sets
        for (int node=0; node<nnodes; node++)
            sets[node] = 0;

        // do unweighted parsimony by postorder traversal
        for (int node=0; node<nnodes; node++) {
            if (node < nleaves) {
                sets[node] = 1 << dna2int[(int)seqs[node][pos]];
            } else {
                char lset = sets[itree[node].child[0]];
                char rset = sets[itree[node].child[1]];
                char intersect = lset & rset;
                if (intersect > 0)
                    sets[node] = intersect;
                else
                    sets[node] = lset | rset;
            }
        }

        // traceback
        // arbitrary choose root base from set
        char rootset = sets[nnodes-1];
        ancestral[nnodes-1] = (rootset & 1) ? int2dna[0] :
            (rootset & 2) ? int2dna[1] :
            (rootset & 4) ? int2dna[2] : int2dna[3];
        
        // traceback with preorder traversal
        for (int node=nnodes-2; node>-1; node--) {
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
                pchar = ancestral[itree[node].parent];
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


void calc_emissions(intstate *states, int nstates, 
                    int *ptree, int nnodes, double *ages, 
                    char **seqs, int nseqs, int seqlen, 
                    double *times, int ntimes,
                    double mu, double **emit)
{
    const double mintime = times[1];
    const double maxtime = times[ntimes - 1];
    double t1, t2, t2a, t2b, t3;
    double parent_age;
    int parent;
    int newnode = nseqs - 1;

    // create inttree data structure
    intnode *itree = make_itree(nnodes, ptree);

    // infer parsimony ancestral sequences
    char **ancestral = new_matrix<char>(nnodes, seqlen);
    parsimony_ancestral_seq_all(itree, nnodes, seqs, nseqs, seqlen, ancestral);


    // base variables
    // v = new chromosome
    // x = current branch
    // p = parent of current branch
    char v, x, p;

    // iterate through positions
    for (int i=0; i<seqlen; i++) {
        v = seqs[newnode][i];
        
        // iterate through states
        for (int j=0; j<nstates; j++) {
            int node = states[j][0];
            int timei = states[j][1];
            double time = times[timei];
            double node_age = ages[node];

            x = ancestral[node][i];

            if (itree[node].parent != -1) {
                parent = itree[node].parent;
                parent_age = ages[parent];

                if (itree[parent].parent == -1) {
                    // unwrap top branch
                    int *c = itree[parent].child;
                    int sib = (node == c[0] ? c[1] : c[0]);
                    p = ancestral[sib][i];

                    // modify (x,p) length to (x,p) + (sib,p)
                    parent_age = 2 * parent_age - ages[sib];

                } else {
                    p = ancestral[parent][i];
                }
            } else {
                // adjust time by unwrapping branch e(v)
                parent = -1;
                parent_age = -1;
                time = 2 * time - node_age;
                p = x;
            }

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

    free_itree(itree);
    delete_matrix(ancestral, nnodes);
}


void calc_emissions2(const States &states,
                    int *ptree, int nnodes, int *ages_index, 
                    char **seqs, int nseqs, int seqlen, 
                    double *times, int ntimes,
                    double mu, double **emit)
{
    const double mintime = times[1];
    const double maxtime = times[ntimes - 1];
    double t1, t2, t2a, t2b, t3;
    double parent_age;
    int parent;
    int newnode = nseqs - 1;
    

    double ages[nnodes];
    for (int i=0; i<nnodes; i++) {
        ages[i] = times[ages_index[i]];
        //rintf("%d %d %f\n", i, ages_index[i], ages[i]);
    }

    // create inttree data structure
    intnode *itree = make_itree(nnodes, ptree);

    // infer parsimony ancestral sequences
    //char **ancestral = new_matrix<char>(nnodes, seqlen);
    //parsimony_ancestral_seq(itree, nnodes, seqs, nseqs, seqlen, ancestral);
    char ancestral[nnodes];


    // base variables
    // v = new chromosome
    // x = current branch
    // p = parent of current branch
    char v, x, p;

    // iterate through positions
    for (int i=0; i<seqlen; i++) {
        v = seqs[newnode][i];

        parsimony_ancestral_seq(itree, nnodes, seqs, nseqs, seqlen, 
                                 i, ancestral);
        
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

            if (itree[node].parent != -1) {
                parent = itree[node].parent;
                parent_age = ages[parent];

                if (itree[parent].parent == -1) {
                    // unwrap top branch
                    int *c = itree[parent].child;
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

    free_itree(itree);
}


double **new_emissions(intstate *istates, int nstates, 
                       int *ptree, int nnodes, int *ages_index,
                       char **seqs, int nseqs, int seqlen, 
                       double *times, int ntimes,
                       double mu)
{
    double **emit = new_matrix<double>(seqlen, nstates);

    States states;
    make_states(istates, nstates, states);
    //LocalTree tree(ptree, nnodes, ages);

    double ages[nnodes];
    for (int i=0; i<nnodes; i++) {
        ages[i] = times[ages_index[i]];
        //printf("%d %d %f\n", i, ages_index[i], ages[i]);
    }

    /*
    calc_emissions(istates, nstates,
                   ptree, nnodes, ages,
                   seqs, nseqs, seqlen,
                   times, ntimes,
                   mu, emit);
    */

    calc_emissions2(states, 
                   ptree, nnodes, ages_index,
                   seqs, nseqs, seqlen,
                   times, ntimes,
                   mu, emit);
    

    return emit;
}


void delete_emissions(double **emit, int seqlen)
{
    delete_matrix<double>(emit, seqlen);
}



//=============================================================================



void arghmm_forward_alg(int **ptrees, int **ages, int **sprs, int *blocklens,
                        int ntrees, int nnodes, 
                        double *times, int ntimes,
                        double *popsizes, double rho, double mu,
                        double **fw)
{
    
    // allocate lineage counts
    LineageCounts lineages(ntimes);
    
    // setup model
    ArgModel model(ntimes, times, popsizes, rho, mu);
    
    // setup local trees
    LocalTrees local_trees(ptrees, ages, sprs, blocklens, ntrees, nnodes);
    
    // state spaces
    States states1;
    States states2;
    States *states = &states1;
    States *last_states = NULL;
    

    double **emit;
    
    // iterate over local trees
    for (LocalTrees::iterator it=local_trees.begin();
         it != local_trees.end(); it++)
    {
        int pos = it->block.start;
        LocalTree *tree = it->tree;
        get_coal_states(tree, ntimes, *states);
        lineages.count(tree);
        
        // TODO: calculate emissions


        // use switch matrix for first column of forward table
        // if we have a previous state space (i.e. not first block)
        if (!last_states) {
            // calculate prior of first state
            calc_state_priors(*states, &lineages, &model, fw[0]);  

        } else {
            double *col1 = fw[pos-1];
            double *col2 = fw[pos];

            int nstates1 = last_states->size();
            int nstates2 = states->size();
            double tmp[nstates1];

            // TODO: calculate transmat_switch
            double **transmat_switch;

            for (int k=0; k<nstates2; k++) {
                double e = emit[0][k];
                for (int j=0; j<nstates1; j++)
                    tmp[j] = col1[j] + transmat_switch[j][k] + e;
                col2[k] = logsum(tmp, nstates1);
            }
            
        }

        // update state space pointers
        last_states = states;
        states = ((states == &states1) ? &states2 : &states1);
    }

    /*
        # iterate over blocks
    for block, nstates, transmat, transmat_switch, emit in matrices:
        if verbose:
            util.logger(" pos %d" % block[0])

        blocklen = block[1] - block[0]

        # use switch matrix for first col
        if block[0] > 0:
            nstates1 = len(transmat_switch)
            nstates2 = len(transmat_switch[0])
            
            col1 = probs[-1]
            col2 = []
            for k in xrange(nstates2):
                e = emit[0][k]
                col2.append(stats.logsum([col1[j] + transmat_switch[j][k] + e
                                          for j in xrange(nstates1)]))
            probs.append(col2)

        # use transmat for rest of block
        # make forward table for block
        fw = [probs[-1]]
        for pos in xrange(block[0]+1, block[1]):
            fw.append([0.0 for k in xrange(nstates)])
        
        forward_alg(blocklen, nstates, transmat, emit, fw)

        if not matrices_given:
            delete_emissions(emit, blocklen)
            delete_transition_probs(transmat, nstates)
        
        for col in fw[1:]:
            probs.append(col[:nstates])
*/


}



} // extern C

}

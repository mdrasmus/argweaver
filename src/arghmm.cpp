
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

typedef int State[2];



//=============================================================================
// ArgHMM model

class Block 
{
public:
    Block(int start, int end) : 
        start(start), end(end) {}
    int start;
    int end;
};


class LocalTrees 
{
public:
    LocalTrees() {}
    ~LocalTrees() {}

    int nnodes;
    list<intnode *> trees;
    list<Block> blocks;
    
};


class ArgHmm 
{
public:
    ArgHmm() {}
    ~ArgHmm() {}

    // time points
    int ntimes;
    double *times;
    double *time_steps;

    // parameters
    double *popsizes;
    double rho;
    double mut;

};


class NodeStateLookup
{
public:
    NodeStateLookup(State *states, int nstates, int nnodes) :
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
            node_mintimes[i] = min(node_mintimes[i], states[i][1]);
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



//=============================================================================
// transitions



void calc_transition_probs(int nnodes, int *ages_index, double treelen,
                           State *states, int nstates,
                           int ntimes, double *times, double *time_steps,
                           int *nbranches, int *nrecombs, int *ncoals, 
                           double *popsizes, double rho,
                           double **transprob)
{
    const int root = nnodes - 1;
    const double root_age = times[ages_index[root]];
    

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
            treelen2 += time_steps[ages_index[root]]; // add basal branch
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
    State *states1, int nstates1, State *states2, int nstates2,
    int ntimes, double *times, int *next_states)
{

    // recomb_node in tree and last_tree
    // coal_node in last_tree


    // make state lookup
    NodeStateLookup state2_lookup(states2, nstates2, nnodes);

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
            int ignore = -1;
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
                    ignore = node1;
                } else if (recomb_name == child2) {
                    // left child is not disrupted
                    node2 = child1;
                    ignore = node1;
                } else {
                    // node is not disrupted
                    node2 = node1;
                }
            }

            // optionally walk up
            if ((coal_name == node1 or coal_name == node2) and
                coal_time < time1)
            {
                // coal occurs under us
                // TODO: make this probabilistic
                node2 = itree[node2].parent;
                if (node2 == ignore)
                    node2 = itree[node2].parent;
            }
            
            // set next state
            next_states[i] = state2_lookup.lookup(node2, time1);

        } else {
            // SPR should not be able to coal back onto same branch
            // this would be a self cycle
            assert (coal_name != node1);
                
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
            if (other == coal_name)
                node2 = itree[other].parent;
            else
                node2 = other;
               
            // set next state 
            next_states[i] = state2_lookup.lookup(node2, time2);
        }
    }
    
}



void calc_transition_probs_switch(
    int *ptree, int *last_ptree, int nnodes, 
    int recomb_name, int recomb_time, int coal_name, int coal_time,
    int *ages_index, int *last_ages_index,
    double treelen, double last_treelen,
    State *states1, int nstates1,
    State *states2, int nstates2,
    
    int ntimes, double *times, double *time_steps,
    int *nbranches, int *nrecombs, int *ncoals, 
    double *popsizes, double rho,
    double **transprob)
{

    // create inttree data structures
    intnode *last_itree = make_itree(nnodes, last_ptree);
    intnode *itree = make_itree(nnodes, ptree);


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
            // deterministic transition
            assert (determ[i] != -1);
            transprob[i][determ[i]] = 0.0;

        } else {
            // probabilistic transition case
            
    }

    /*
            # probabilistic transition case
            
            # \frac{k^{(n-1)}_{j-1}}{k^{(n-1)}_{j-1} + 1}
            # \frac{[1 - \exp(- s'_{j-1} k^{(n)}_{j-1} / (2N))]}
            #      {[1 - \exp(- s'_{j-1} k^{(n-1)}_{j-1} / (2N))]}
            # \frac{|T^{n-1}_{i-1}|
            #       [1 - \exp(- \rho (|T^{n-1}_{i-1}| + t_{i-1}))]}
            #      {[|T^{n-1}_{i-1}| + t_{i-1}]
            #       [1 - \exp(- \rho |T^{n-1}_{i-1}|)]} 
            # \exp(- \sum_{m=k}^{j-2} s'_k / (2N))
            
            if (node1, a) in states2:
                j = states2.index((node1, a))

                # TODO: add ncoals and nrecomb
                # - 1 lineage if recomb br in time segment b-1
                b = a                
                if (recomb_branch, b) in states2:
                    kbn1 = max(nbranches[b-1] - 1, 1.0)
                else:
                    kbn1 = nbranches[b-1]
                kbn  = kbn1 + 1

                transprob[i][j] = (
                    (kbn1/float(kbn)) *
                    ((1.0 - exp(-time_steps[b-1] * kbn /
                                (2.0*popsizes[b-1])))/
                     (1.0 - exp(-time_steps[b-1] * kbn1 /
                                (2.0*popsizes[b-1]))))*
                    (treelen / (treelen + times[a])) *
                    ((1.0 - exp(-rho * (treelen + times[a]))) /
                     (1.0 - exp(-rho * treelen))) *
                    exp(- sum(time_steps[m] / (2.0 * popsizes[m])
                              for m in xrange(k, b-1))))

            for j, (node2, b) in enumerate(states2):
                transprob[i][j] = 0.0
                if node2 != recomb_branch:
                    continue

                # require coal above recomb
                if b < k:
                    continue

                # TODO: fix
                # - 1 lineage if recomb br in time segment b-1
                if (recomb_branch, b) in states2:
                    kbn1 = max(nbranches[b-1] - 1, 1.0)
                else:
                    kbn1 = nbranches[b-1]
                kbn  = kbn1 + 1

                transprob[i][j] = (
                    (kbn1/float(kbn)) *
                    ((1.0 - exp(-time_steps[b-1] * kbn /
                                (2.0*popsizes[b-1])))/
                     (1.0 - exp(-time_steps[b-1] * kbn1 /
                                (2.0*popsizes[b-1]))))*
                    (treelen / (treelen + times[a])) *
                    ((1.0 - exp(-rho * (treelen + times[a]))) /
                     (1.0 - exp(-rho * treelen))) *
                    exp(- sum(time_steps[m] / (2.0 * popsizes[m])
                              for m in xrange(k, b-1))))            

            # HACK for now:  renormalize row to ensure they add up to one
            tot = sum(transprob[i])
            for j in xrange(len(states2)):
                x = transprob[i][j]
                if tot > 0.0 and x > 0.0:
                    transprob[i][j] = log(x / tot)
                else:
                    transprob[i][j] = -1e1000
    */
}




double **new_transition_probs(int nnodes, int *ages_index, double treelen,
                              State *states, int nstates,
                              int ntimes, double *times, double *time_steps,
                              int *nbranches, int *nrecombs, int *ncoals, 
                              double *popsizes, double rho)
{
    double **transprob = new_matrix<double>(nstates, nstates);
    calc_transition_probs(nnodes, ages_index, treelen, 
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
    State *states1, int nstates1,
    State *states2, int nstates2,
    
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



//============================================================================
// emissions

void parsimony_ancestral_seq(intnode *itree, int nnodes, char **seqs, 
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


void calc_emissions(State *states, int nstates, 
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
    parsimony_ancestral_seq(itree, nnodes, seqs, nseqs, seqlen, ancestral);


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
                    int sib = (node == c[0] ? c[1] : c[1]);
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


double **new_emissions(State *states, int nstates, 
                       int *ptree, int nnodes, double *ages,
                       char **seqs, int nseqs, int seqlen, 
                       double *times, int ntimes,
                       double mu)
{
    double **emit = new_matrix<double>(seqlen, nstates);
    calc_emissions(states, nstates, 
                   ptree, nnodes, ages,
                   seqs, nseqs, seqlen,
                   times, ntimes,
                   mu, emit);
    return emit;
}


void delete_emissions(double **emit, int seqlen)
{
    delete_matrix<double>(emit, seqlen);
}



} // extern C

}

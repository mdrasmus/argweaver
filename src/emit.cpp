
#include "common.h"
#include "emit.h"
#include "seq.h"

namespace arghmm {


//============================================================================
// emissions


void parsimony_ancestral_seq(const LocalTree *tree, const char * const *seqs, 
                             int nseqs, int pos, char *ancestral) 
{
    const int nnodes = tree->nnodes;
    const LocalNode *nodes = tree->nodes;
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




void calc_emissions(const States &states, const LocalTree *tree,
                    const char *const *seqs, int nseqs, int seqlen, 
                    const ArgModel *model, double **emit)
{
    const double *times = model->times;
    const double mintime = times[1];
    const double maxtime = times[model->ntimes - 1];
    const double mu = model->mu;
    const int nnodes = tree->nnodes;
    const LocalNode *nodes = tree->nodes;
    
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

        // check for no mutation case
        char c = seqs[0][i];
        bool mut = false;
        for (int j=1; j<nseqs; j++) {
            if (seqs[j][i] != c) {
                mut = true;
                break;
            }
        }

        // handle no mutation case
        if (!mut) {
            for (unsigned int j=0; j<states.size(); j++) {
                int node = states[j].node;
                int timei = states[j].time;
                double time = times[timei];
                double node_age = ages[node];

                if (nodes[node].parent == -1) {
                    // adjust time by unwrapping branch e(v)
                    time += time - node_age;
                }

                emit[i][j] = - mu * max(time, mintime);
            }
            continue;
        }

        parsimony_ancestral_seq(tree, seqs, nseqs, i, ancestral);
        
        // iterate through states
        for (unsigned int j=0; j<states.size(); j++) {
            int node = states[j].node;
            int timei = states[j].time;
            double time = times[timei];
            double node_age = ages[node];

            x = ancestral[node];

            // get bases and ages
            if (nodes[node].parent != -1) {
                parent = nodes[node].parent;
                parent_age = ages[parent];

                if (nodes[parent].parent == -1) {
                    // unwrap top branch
                    const int *c = nodes[parent].child;
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
            

            // ensure mintime
            if (time < mintime) 
                time = mintime;

            // handle cases
            if (v == x && x == p) {
                // no mutation
                emit[i][j] = - mu * time;

            } else if (v != p && p == x) {
                // mutation on v
                emit[i][j] = log(.333333 - .333333 * exp(-mu * time));

            } else if (v == p && p != x) {
                // mutation on x
                t1 = max(parent_age - node_age, mintime);
                t2 = max(time - node_age, mintime);

                emit[i][j] = log((1 - exp(-mu *t2)) / (1 - exp(-mu * t1)))
                    -mu * (time + t2 - t1);

            } else if (v == x && x != p) {
                // mutation on (y,p)
                t1 = max(parent_age - node_age, mintime);
                t2 = max(parent_age - time, mintime);

                emit[i][j] = log((1 - exp(-mu * t2)) / (1 - exp(-mu * t1)))
                    -mu * (time + t2 - t1);

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
                                 / (1 - exp(-mu * t1)))
                    -mu * (time + t2 + t3 - t1);
            }
        }
    }
}



void calc_emissions_internal(const States &states, const LocalTree *tree,
                             const char *const *seqs, int nseqs, int seqlen, 
                             const ArgModel *model, double **emit)
{
    const double *times = model->times;
    const double mintime = times[1];
    const double maxtime = times[model->ntimes - 1];
    const double mu = model->mu;
    const int nnodes = tree->nnodes;
    const LocalNode *nodes = tree->nodes;

    // ignore fully specified local tree
    if (states.size() == 0)
        return;
    
    double t1, t2, t2a, t2b, t3;
    double parent_age;
    int parent;
    int subtree_root = nodes[tree->root].child[0];
    double subtree_root_age = times[nodes[subtree_root].age];
    
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
        // check for no mutation case
        char c = seqs[0][i];
        bool mut = false;
        for (int j=1; j<nseqs; j++) {
            if (seqs[j][i] != c) {
                mut = true;
                break;
            }
        }

        // handle no mutation case
        if (!mut) {
            for (unsigned int j=0; j<states.size(); j++) {
                int node = states[j].node;
                int timei = states[j].time;
                double time = times[timei] - subtree_root_age;
                double node_age = ages[node];

                if (nodes[node].parent == -1) {
                    // adjust time by unwrapping branch e(v)
                    time += time - node_age;
                }
                time = max(time, mintime);

                emit[i][j] = - mu * time;
            }
            continue;
        }

        // TODO: don't do parsimony directly on subtree-maintree represenation
        parsimony_ancestral_seq(tree, seqs, nseqs, i, ancestral);
        v = ancestral[subtree_root];

        
        // iterate through states
        for (unsigned int j=0; j<states.size(); j++) {
            int node = states[j].node;
            int timei = states[j].time;
            double time = times[timei] - subtree_root_age;
            double node_age = ages[node];

            x = ancestral[node];

            // get bases and ages
            if (nodes[node].parent != -1) {
                parent = nodes[node].parent;
                parent_age = ages[parent];

                if (nodes[parent].parent == -1) {
                    // unwrap top branch
                    const int *c = nodes[parent].child;
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
                time += time - node_age;
                p = x;
            }
            

            // ensure mintime
            if (time < mintime) 
                time = mintime;

            // handle cases
            if (v == x && x == p) {
                // no mutation
                emit[i][j] = - mu * time;

            } else if (v != p && p == x) {
                // mutation on v
                emit[i][j] = log(.333333 - .333333 * exp(-mu * time));

            } else if (v == p && p != x) {
                // mutation on x
                t1 = max(parent_age - node_age, mintime);
                t2 = max(time - node_age, mintime);

                emit[i][j] = log((1 - exp(-mu *t2)) / (1 - exp(-mu * t1)))
                    -mu * (time + t2 - t1);

            } else if (v == x && x != p) {
                // mutation on (y,p)
                t1 = max(parent_age - node_age, mintime);
                t2 = max(parent_age - time, mintime);

                emit[i][j] = log((1 - exp(-mu * t2)) / (1 - exp(-mu * t1)))
                    -mu * (time + t2 - t1);

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
                                 / (1 - exp(-mu * t1)))
                    -mu * (time + t2 + t3 - t1);
            }
        }
    }
}




//=============================================================================
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



} // namespace arghmm

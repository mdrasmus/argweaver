
#include "common.h"
#include "emit.h"
#include "seq.h"
#include "thread.h"

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


// Juke-Cantor
inline double prob_branch(double t, double mu, bool mut)
{
    const double f = 4. / 3.;
    if (!mut)
        return .25 * (1.0 + 3. * exp(-f*mu*t));
    else
        return .25 * (1.0 - exp(-f*mu*t));
}


double likelihood_site(const LocalTree *tree, const char *const *seqs, 
                       int pos, const int *order, const double *muts,
                       const double *nomuts)
{
    const int nnodes = tree->nnodes;
    const LocalNode* nodes = tree->nodes;
    double table[nnodes][4];

    // iterate postorder through nodes
    for (int i=0; i<nnodes; i++) {
        const int j = order[i];
        
        if (nodes[j].is_leaf()) {
            // leaf case
            const char c = seqs[j][pos];
            table[j][0] = 0.0;
            table[j][1] = 0.0;
            table[j][2] = 0.0;
            table[j][3] = 0.0;
            table[j][dna2int[(int) c]] = 1.0;
        } else {
            // internal node case
            int c1 = nodes[j].child[0];
            int c2 = nodes[j].child[1];
            
            for (int a=0; a<4; a++) {
                double p1 = 0.0;
                double p2 = 0.0;
                
                for (int b=0; b<4; b++) {
                    if (a == b) {
                        p1 += table[c1][b] * nomuts[c1];
                        p2 += table[c2][b] * nomuts[c2];
                    } else {
                        p1 += table[c1][b] * muts[c1];
                        p2 += table[c2][b] * muts[c2];
                    }
                }

                table[j][a] = p1 * p2;
            }
        }
    }

    // sum over root node
    double p = 0.0;
    int root = tree->root;
    for (int a=0; a<4; a++)
        p += table[root][a]  * .25;
    
    return log(p);
}


void likelihood_sites(const LocalTree *tree, const ArgModel *model,
                      const char *const *seqs, 
                      const int seqlen, const int statei, 
                      const bool *invariant,
                      double **emit)
{
    const double *times = model->times;
    const LocalNode *nodes = tree->nodes;
    const double mintime = times[1];
    double invariant_lk = -1;

    // get postorder
    int order[tree->nnodes];
    tree->get_postorder(order);

    // get mutation probabilities
    double muts[tree->nnodes];
    double nomuts[tree->nnodes];
    for (int i=0; i<tree->nnodes; i++) {
        if (i != tree->root) {
            double t = max(times[nodes[nodes[i].parent].age] - 
                           times[nodes[i].age], mintime);
            muts[i] = prob_branch(t, model->mu, true);
            nomuts[i] = prob_branch(t, model->mu, false);
        }
    }

    // calculate emissions for tree at each site
    for (int i=0; i<seqlen; i++) {
        if (invariant && invariant[i] && invariant_lk > 0)
            emit[i][statei] = invariant_lk;
        else {
            emit[i][statei] = likelihood_site(tree, seqs, i, order, 
                                              muts, nomuts);
            // save invariant likelihood
            if (invariant && invariant[i]) {
                invariant_lk = emit[i][statei];
            }
        }
    }
}


void find_invariant_sites(const char *const *seqs, int nseqs, int seqlen, 
                          bool *invariant)
{
    // find invariant sites
    for (int i=0; i<seqlen; i++) {
        char c = seqs[0][i];
        bool mut = false;
        for (int j=1; j<nseqs; j++) {
            if (seqs[j][i] != c) {
                mut = true;
                break;
            }
        }
        invariant[i] = !mut;
    }
}


void calc_emissions2(const States &states, const LocalTree *tree,
                     const char *const *seqs, int nseqs, int seqlen, 
                     const ArgModel *model, double **emit)
{
    const int nstates = states.size();
    const int newleaf = tree->get_num_leaves();
    bool *invariant = new bool [seqlen];

    // create local tree we can edit
    LocalTree tree2(tree->nnodes, tree->nnodes + 2);
    tree2.copy(*tree);

    // find invariant sites
    find_invariant_sites(seqs, nseqs, seqlen, invariant);

    for (int j=0; j<nstates; j++) {
        State state = states[j];
        add_tree_branch(&tree2, state.node, state.time);
        likelihood_sites(&tree2, model, seqs, seqlen, j, invariant, emit);
        remove_tree_branch(&tree2, newleaf, NULL);
    }

    delete [] invariant;
}


void calc_emissions_internal2(const States &states, const LocalTree *tree,
                     const char *const *seqs, int nseqs, int seqlen, 
                     const ArgModel *model, double **emit)
{
    const int nstates = states.size();
    const int subtree_root = tree->nodes[tree->root].child[0];
    const int subtree_root_age = tree->nodes[subtree_root].age;
    const int maxtime = model->ntimes + 1;
    bool *invariant = new bool [seqlen];

    // ignore fully specified local tree
    if (nstates == 0)
        return;

    // create local tree we can edit
    LocalTree tree2(tree->nnodes, tree->nnodes + 2);
    tree2.copy(*tree);

    // find invariant sites
    find_invariant_sites(seqs, nseqs, seqlen, invariant);

    for (int j=0; j<nstates; j++) {
        State state = states[j];
        assert(subtree_root != tree2.root);
        Spr add_spr(subtree_root, subtree_root_age, state.node, state.time);
        apply_spr(&tree2, add_spr);

        likelihood_sites(&tree2, model, seqs, seqlen, j, invariant, emit);

        Spr remove_spr(subtree_root, subtree_root_age, 
                       tree2.root, maxtime);
        apply_spr(&tree2, remove_spr);
    }

    delete [] invariant;
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

                if (nodes[node].parent == tree->root) {
                    // adjust time by unwrapping branch e(v)
                    time += times[timei] - node_age;
                }
                time = max(time, mintime);

                emit[i][j] = - mu * time;
                assert(!isnan(emit[i][j]));
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
            double time = times[timei];
            double node_age = ages[node];
            double blen = time - subtree_root_age;

            x = ancestral[node];

            // get bases and ages
            if (nodes[node].parent != tree->root) {
                parent = nodes[node].parent;
                parent_age = ages[parent];

                if (nodes[parent].parent == tree->root) {
                    // unwrap top branch
                    const int *c = nodes[parent].child;
                    int sib = (node == c[0] ? c[1] : c[0]);
                    p = ancestral[sib];

                    // modify (x,p) length to (x,p) + (sib,p)
                    parent_age += parent_age - ages[sib];
                } else {
                    p = ancestral[parent];
                }
            } else {
                // adjust time by unwrapping branch e(v)
                parent = tree->root;
                parent_age = -1;
                blen += time - node_age;
                p = x;
            }
            

            // ensure mintime
            if (blen < mintime) 
                blen = mintime;

            // handle cases
            if (v == x && x == p) {
                // no mutation
                emit[i][j] = - mu * blen;

            } else if (v != p && p == x) {
                // mutation on v
                emit[i][j] = log(.333333 - .333333 * exp(-mu * blen));

            } else if (v == p && p != x) {
                // mutation on x
                t1 = max(parent_age - node_age, mintime);
                t2 = max(time - node_age, mintime);

                emit[i][j] = log((1 - exp(-mu *t2)) / (1 - exp(-mu * t1)))
                    -mu * (blen + t2 - t1);

            } else if (v == x && x != p) {
                // mutation on (y,p)
                t1 = max(parent_age - node_age, mintime);
                t2 = max(parent_age - time, mintime);

                emit[i][j] = log((1 - exp(-mu * t2)) / (1 - exp(-mu * t1)))
                    -mu * (blen + t2 - t1);

            } else {
                // two mutations (v,x)
                // mutation on x
                if (parent != tree->root) {
                    t1 = max(parent_age - node_age, mintime);
                    t2a = max(parent_age - time, mintime);
                } else {
                    t1 = max(maxtime - node_age, mintime);
                    t2a = max(maxtime - time, mintime);
                }
                t2b = max(time - node_age, mintime);
                t2 = max(t2a, t2b);
                t3 = blen;

                emit[i][j] = log((1 - exp(-mu *t2)) * (1 - exp(-mu *t3))
                                 / (1 - exp(-mu * t1)))
                    -mu * (blen + t2 + t3 - t1);
            }

            assert(!isnan(emit[i][j]));
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

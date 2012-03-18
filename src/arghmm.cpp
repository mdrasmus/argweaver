
#include "common.h"
#include "itree.h"
#include "ptree.h"
#include "seq.h"

using namespace spidir;
using namespace dlcoal;

namespace arghmm {

extern "C" {

typedef int State[2];



//=============================================================================
// ArgHMM model




//=============================================================================
// transitions

void calc_transition_probs(int nnodes, int *ages_index, double treelen,
                           State *states, int nstates,
                           int ntimes, double *times, double *time_steps,
                           int *nbranches, int *nrecombs, int *ncoals, 
                           double *popsizes, double rho,
                           double **transprob)
{
    
    // TODO: correctly use root node to calc treelen2
    
    //double time_steps[ntimes];
    //for (int i=1; i<=ntimes; i++)
    //    time_steps[i-1] = times[i] - times[i-1];
    const double mintime = time_steps[0];
    const int root = nnodes - 1;
    const double root_age = times[ages_index[root]];
    const double topbranch = time_steps[ages_index[root]];
    
    
    // A_{k,j} =& s'_{j-1} k_{j-1} / (2N) + \sum_{m=k}^{j-2} s'_m k_m / (2N)
    //         =& s'_{j-1} k_{j-1} / (2N) + A_{k,j-1}.
    double **A = new_matrix<double>(ntimes, ntimes);
    for (int k=0; k<ntimes; k++) {
        A[k][k] = A[k][min(k+1, ntimes-1)] = 0.0;
        for (int j=k+2; j<ntimes; j++) {
            const int l = j - 2;
            A[k][j] = A[k][j-1] + time_steps[l] * nbranches[l] / 
                (2.0 * popsizes[l]);
        }
    }

    // B_{c,a} =& \sum_{k=0}^{c} \exp(- A_{k,a})
    //         =& B_{c-1,a} + \exp(- A_{c,a}).
    double **B = new_matrix<double>(ntimes, ntimes);
    for (int b=0; b<ntimes; b++) {
        B[0][b] = nbranches[0] * time_steps[0] / nrecombs[0] * exp(-A[0][b]);
        for (int c=1; c<b; c++) {
            B[c][b] = B[c-1][b] + nbranches[c] * time_steps[c] / nrecombs[c]
                * exp(-A[c][b]);
        }
    }

    // S_{a,b} &= B_{min(a-1,b-1),b}
    double **S = new_matrix<double>(ntimes, ntimes);
    for (int a=1; a<ntimes; a++) {
        for (int b=1; b<ntimes; b++)
            S[a][b] = B[min(a-1, b-1)][b];
    }
    
    // f =\frac{[1 - \exp(- \rho (|T^{n-1}_{i-1}| + s_a))] 
    //       [1 - \exp(- s'_b k_b / (2N))]}
    //      {\exp(-\rho |T^{n-1}_{i-1}|) (|T^{n-1}_{i-1}| + s_a) k^C_b}
    // |T^{n-1}_{i-1}| = treelen
    for (int i=0; i<nstates; i++) {
        const int node1 = states[i][0];
        const int a = states[i][1];
        const int c = ages_index[node1];
        assert(a < ntimes);

        double treelen2 = treelen + times[a];
        if (node1 == nnodes-1)
            treelen2 += times[a] - root_age;
        if (treelen2 <= treelen)
            // because of discritization we need to enfore
            // nonzero branch change
           treelen2 += treelen + mintime;

        for (int j=0; j<nstates; j++) {
            const int node2 = states[j][0];
            const int b = states[j][1];
            assert(b < ntimes);
            
            //double treelen2 = treelen + max(times[a], mintime);
            double f = ((1.0 - exp(-rho * treelen2)) /
                              (exp(-rho * treelen) * 
                               (treelen2 + topbranch) * ncoals[b]));
            if (b > 0)
                f *= (1.0 - exp(-time_steps[b-1] * nbranches[b-1]
                                / (2.0 * popsizes[b-1])));

            if (node1 != node2)
                transprob[i][j] = f * S[a][b];
            else if (a != b) 
                transprob[i][j] = f * (2*S[a][b] - S[c][b]);
            else
                transprob[i][j] = 0.0;
        }

        double sum = 0.0;
        for (int j=0; j<nstates; j++)
            sum += transprob[i][j];
        transprob[i][i] = 1.0 - sum;
        for (int j=0; j<nstates; j++)
            transprob[i][j] = log(transprob[i][j]);
    }

    delete_matrix<double>(A, ntimes);
    delete_matrix<double>(B, ntimes);
    delete_matrix<double>(S, ntimes);
}


void calc_transition_probs2(int nnodes, int *ages_index, double treelen,
                           State *states, int nstates,
                           int ntimes, double *times, double *time_steps,
                           int *nbranches, int *nrecombs, int *ncoals, 
                           double *popsizes, double rho,
                           double **transprob)
{
    
    // TODO: correctly use root node to calc treelen2
    
    //double time_steps[ntimes];
    //for (int i=1; i<=ntimes; i++)
    //    time_steps[i-1] = times[i] - times[i-1];
    const double mintime = time_steps[0];
    const int root = nnodes - 1;
    const double root_age = times[ages_index[root]];
    const double topbranch = time_steps[ages_index[root]];
    
    
    // A_{k,j} =& s'_{j-1} k_{j-1} / (2N) + \sum_{m=k}^{j-2} s'_m k_m / (2N)
    //         =& s'_{j-1} k_{j-1} / (2N) + A_{k,j-1}.
    double **A = new_matrix<double>(ntimes, ntimes);
    for (int k=0; k<ntimes; k++) {
        A[k][k] = 0.0;
        for (int j=k+1; j<ntimes; j++) {
            const int l = j - 1;
            A[k][j] = A[k][j-1] + time_steps[l] * nbranches[l] / 
                (2.0 * popsizes[l]);
        }
    }

    // B_{c,a} =& \sum_{k=0}^{c} \exp(- A_{k,a})
    //         =& B_{c-1,a} + \exp(- A_{c,a}).
    double **B = new_matrix<double>(ntimes, ntimes);
    for (int b=0; b<ntimes; b++) {
        B[0][b] = nbranches[0] * time_steps[0] / nrecombs[0] * exp(-A[0][b]);
        for (int c=1; c<=min(b,ntimes-2); c++) {
            B[c][b] = B[c-1][b] + nbranches[c] * time_steps[c] / nrecombs[c]
                * exp(-A[c][b]);
        }
    }

    // S_{a,b} &= B_{min(a,b),b}
    double **S = new_matrix<double>(ntimes, ntimes);
    for (int a=0; a<ntimes; a++) {
        for (int b=0; b<ntimes; b++)
            S[a][b] = B[min(a, b)][b];
    }
    
    // f =\frac{[1 - \exp(- \rho (|T^{n-1}_{i-1}| + s_a))] 
    //       [1 - \exp(- s'_b k_b / (2N))]}
    //      {\exp(-\rho |T^{n-1}_{i-1}|) (|T^{n-1}_{i-1}| + s_a) k^C_b}
    // |T^{n-1}_{i-1}| = treelen
    for (int i=0; i<nstates; i++) {
        const int node1 = states[i][0];
        const int a = states[i][1];
        const int c = ages_index[node1];
        assert(a < ntimes);

        double treelen2 = treelen + times[a];
        if (node1 == nnodes-1)
            treelen2 += times[a] - root_age;
        if (treelen2 <= treelen)
            // because of discritization we need to enfore
            // nonzero branch change
            treelen2 += treelen + mintime;

        for (int j=0; j<nstates; j++) {
            const int node2 = states[j][0];
            const int b = states[j][1];
            assert(b < ntimes);
            
            //double treelen2 = treelen + max(times[a], mintime);
            const double f = ((1.0 - exp(-rho * treelen2)) *
                              (1.0 - exp(-time_steps[b] * nbranches[b]
                                         / (2.0 * popsizes[b]))) /
                              (exp(-rho * treelen) * (treelen2 + topbranch) * 
                               ncoals[b]));
            if (node1 != node2)
                transprob[i][j] = f * S[a][b];
            else if (a != b) 
                transprob[i][j] = f * (2*S[a][b] - S[c][b]);
            else
                transprob[i][j] = exp(-rho * (treelen2 - treelen));
        }

        double sum = 0.0;
        for (int j=0; j<nstates; j++)
            sum += transprob[i][j];
        for (int j=0; j<nstates; j++)
            transprob[i][j] = log(transprob[i][j] / sum);
    }

    delete_matrix<double>(A, ntimes);
    delete_matrix<double>(B, ntimes);
    delete_matrix<double>(S, ntimes);
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


#include "common.h"
#include "itree.h"
#include "ptree.h"
#include "seq.h"

using namespace spidir;
using namespace dlcoal;

namespace arghmm {

extern "C" {

typedef int State[2];


void forward_step(int i, double *col1, double* col2, 
                  int nstates1, int nstates2, double **trans, double *emit)
{
    for (int k=0; k<nstates2; k++) {
        double tot = -INFINITY;
        for (int j=0; j<nstates1; j++) {
            tot = logadd(tot, col1[j] + trans[j][k]);
        }
        col2[k] = tot + emit[k];
    }
}


void forward_alg(int n, int nstates1, int nstates2, double **fw, 
                 double **trans, double **emit)
{
    double *vec = new double [nstates1];

    for (int i=1; i<n; i++) {
        double *col1 = fw[i-1];
        double *col2 = fw[i];
        double *emit2 = emit[i];

        for (int k=0; k<nstates2; k++) {
            for (int j=0; j<nstates1; j++)
                vec[j] = col1[j] + trans[j][k];
            col2[k] = logsum(vec, nstates1) + emit2[k];
        }
    }

    delete [] vec;
}


void backward_alg(int n, int nstates1, int nstates2, double **bw, 
                  double **trans, double **emit)
{
    double *vec = new double [nstates2];

    for (int i=n-2; i>-1; i--) {
        double *col1 = bw[i];
        double *col2 = bw[i+1];
        double *emit2 = emit[i+1];

        for (int j=0; j<nstates1; j++) {
            for (int k=0; k<nstates2; k++)
                vec[k] = trans[j][k] + col2[k] + emit2[k];
            col1[j] = logsum(vec, nstates2);
        }
    }

    delete []  vec;
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




//=========================================================================

void forward_alg2(int n, int nstates1, int nstates2, double **fw, 
                 double **trans, double **emit)
{
    for (int i=1; i<n; i++) {
        double *col1 = fw[i-1];
        double *col2 = fw[i];
        double *emit2 = emit[i];

        for (int k=0; k<nstates2; k++) {
            double tot = -INFINITY;
            for (int j=0; j<nstates1; j++) {
                tot = logadd(tot, col1[j] + trans[j][k]);
            }
            col2[k] = tot + emit2[k];
        }
    }
}



void backward_alg2(int n, int nstates1, int nstates2, double **bw, 
                  double **trans, double **emit)
{
    for (int i=n-2; i>-1; i--) {
        double *col1 = bw[i];
        double *col2 = bw[i+1];
        double *emit2 = emit[i+1];

        for (int j=0; j<nstates1; j++) {
            double tot = -INFINITY;
            for (int k=0; k<nstates2; k++) {
                tot = logadd(tot, trans[j][k] + col2[k] + emit2[k]);
            }
            col1[j] = tot;
        }
    }
}



} // extern C

}


#include "common.h"

using namespace spidir;

namespace arghmm {

extern "C" {

//=============================================================================
// Hidden Markov Models (HMM)

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


void forward_alg(int n, int nstates, double **trans, double **emit, 
                 double **fw)
{
    double *vec = new double [nstates];

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

    delete [] vec;
}


void backward_alg(int n, int nstates, double **trans, double **emit, 
                  double **bw)
{
    double *vec = new double [nstates];

    for (int i=n-2; i>-1; i--) {
        double *col1 = bw[i];
        double *col2 = bw[i+1];
        double *emit2 = emit[i+1];

        for (int j=0; j<nstates; j++) {
            for (int k=0; k<nstates; k++)
                vec[k] = trans[j][k] + col2[k] + emit2[k];
            col1[j] = logsum(vec, nstates);
        }
    }

    delete []  vec;
}


void sample_hmm_posterior(int n, int nstates, double **trans, double **emit, 
                          double **fw, int *path)
{
    // NOTE: path[n-1] must already be sampled
    
    double *A = new double [nstates];
    double *C = new double [nstates];
    
    //printf("path[n-1] = %d, %d\n", path[n-1], nstates);

    // recurse
    for (int i=n-2; i>=0; i--) {
        double *e = emit[i+1];
        int k = path[i+1];
        for (int j=0; j<nstates; j++) {
            // C_{i,j} = trans(j, Y[i+1]) * emit(X[i+1], Y[i+1])
            // A_{j,i} = F_{i,j} C_{i,j}
            C[j] = trans[j][k] + e[k];
            A[j] = fw[i][j] + C[j];
        }
        double total = logsum(A, nstates);
        for (int j=0; j<nstates; j++)
            A[j] = exp(A[j] - total);
        path[i] = sample(A, nstates);
    }

    delete [] A;
    delete [] C;
}


//=========================================================================
// old code

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
} // namespace arghmm




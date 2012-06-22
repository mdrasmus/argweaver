
#include "common.h"

namespace arghmm {

extern "C" {

//=============================================================================
// Hidden Markov Models (HMM)


// run forward algorithm for one column of the table
void forward_step(double *col1, double* col2, 
                  int nstates1, int nstates2, double **trans, double *emit)
{
    double tmp[nstates1];
    
    for (int k=0; k<nstates2; k++) {
        for (int j=0; j<nstates1; j++)
            tmp[j] = col1[j] + trans[j][k];
        col2[k] = logsum(tmp, nstates1) + emit[k];
    }
}

// run forward algorithm 
// NOTE: first column of fw table must already be calculated
void forward_alg(int n, int nstates, double **trans, double **emit, 
                 double **fw)
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


// run backward algorithm
// NOTE: last column of bw table must already be calculated
void backward_alg(int n, int nstates, double **trans, double **emit, 
                  double **bw)
{
    double vec[nstates];

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
}


void sample_hmm_posterior(int n, int nstates, double **trans, 
                          double **fw, int *path)
{
    // NOTE: path[n-1] must already be sampled
    
    double A[nstates];

    // recurse
    for (int i=n-2; i>=0; i--) {
        int k = path[i+1];
        for (int j=0; j<nstates; j++)
            A[j] = fw[i][j] + trans[j][k];
        double total = logsum(A, nstates);
        for (int j=0; j<nstates; j++)
            A[j] = exp(A[j] - total);
        path[i] = sample(A, nstates);
    }
}


int sample_hmm_posterior_step(int nstates1, double **trans, double *col1,
                              int state2)
{
    double A[nstates1];
    
    for (int j=0; j<nstates1; j++)
        A[j] = col1[j] + trans[j][state2];
    double total = logsum(A, nstates1);
    for (int j=0; j<nstates1; j++)
        A[j] = exp(A[j] - total);
    return sample(A, nstates1);
}


} // extern C
} // namespace arghmm




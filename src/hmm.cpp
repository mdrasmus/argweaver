
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

} // extern C
} // namespace arghmm



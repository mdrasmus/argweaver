#ifndef ARGHMM_HMM_H
#define ARGHMM_HMM_H


namespace arghmm {

extern "C" {

//=============================================================================
// Hidden Markov Models (HMM)

void forward_step(int i, double *col1, double* col2, 
                  int nstates1, int nstates2, double **trans, double *emit);
void forward_alg(int n, int nstates, 
                 double **trans, double **emit, double **fw);
void backward_alg(int n, int nstates,
                  double **trans, double **emit, double **bw);
void sample_hmm_posterior(int n, int nstates, double **trans, double **emit, 
                          double **fw, int *path);


} // extern C
} // namespace arghmm



#endif // ARGHMM_HMM_H

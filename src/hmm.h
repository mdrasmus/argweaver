#ifndef ARGHMM_COMMON_H
#define ARGHMM_COMMON_H


namespace arghmm {

extern "C" {

//=============================================================================
// Hidden Markov Models (HMM)

void forward_step(int i, double *col1, double* col2, 
                  int nstates1, int nstates2, double **trans, double *emit);
void forward_alg(int n, int nstates1, int nstates2, double **fw, 
                 double **trans, double **emit);
void backward_alg(int n, int nstates1, int nstates2, double **bw, 
                  double **trans, double **emit);

} // extern C
} // namespace arghmm



#endif // ARGHMM_COMMON_H

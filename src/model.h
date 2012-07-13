//=============================================================================
// ArgHmm model


#ifndef ARGHMM_MODEL_H
#define ARGHMM_MODEL_H

#include "common.h"

namespace arghmm {


// The model parameters and time discretization scheme
class ArgModel 
{
public:
    ArgModel(int ntimes, double *_times, double *_popsizes, 
             double rho, double mu) :
        ntimes(ntimes),
        rho(rho),
        mu(mu)
    {
        times = new double [ntimes];
        copy(_times, _times + ntimes, times);

        if (_popsizes) {
            popsizes = new double [ntimes];
            copy(_popsizes, _popsizes + ntimes, popsizes);
        } else {
            popsizes = NULL;
        }

        // setup time steps
        time_steps = new double [ntimes];
        for (int i=0; i<ntimes-1; i++)
            time_steps[i] = times[i+1] - times[i];
        time_steps[ntimes-1] = INFINITY;
    }
    ~ArgModel()
    {
        delete [] times;
        delete [] time_steps;
        if (popsizes)
            delete [] popsizes;
    }

    // time points
    int ntimes;
    double *times;
    double *time_steps;

    // parameters
    double *popsizes;
    double rho;
    double mu;
};



} // namespace arghmm


#endif // ARGHMM_MODEL_H

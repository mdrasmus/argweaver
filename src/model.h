//=============================================================================
// ArgHmm model


#ifndef ARGHMM_MODEL_H
#define ARGHMM_MODEL_H

#include "common.h"

namespace arghmm {


// Returns a discretized time point
inline double get_time_point(int i, int ntimes, double maxtime, double delta=10)
{
    return (exp(i/double(ntimes) * log(1.0 + delta * maxtime)) - 1) / delta;
}

// Returns a list of discretized time points
inline void get_time_points(int ntimes, double maxtime, 
                            double *times, double delta=.01)
{
    for (int i=0; i<ntimes; i++)
        times[i] = get_time_point(i, ntimes-1, maxtime, delta);
}



// The model parameters and time discretization scheme
class ArgModel 
{
public:
    ArgModel(int ntimes, double maxtime, double popsize, 
             double rho, double mu) :
        ntimes(ntimes),
        rho(rho),
        mu(mu)
    {
        times = new double [ntimes+1];
        get_time_points(ntimes, maxtime, times);
        
        popsizes = new double [ntimes];
        fill(popsizes, popsizes + ntimes, popsize);

        // setup time steps
        setup_time_steps();
    }

    ArgModel(int ntimes, double maxtime, double *_popsizes, 
             double rho, double mu) :
        ntimes(ntimes),
        rho(rho),
        mu(mu)
    {
        times = new double [ntimes];
        get_time_points(ntimes, maxtime, times);

        if (_popsizes) {
            popsizes = new double [ntimes];
            copy(_popsizes, _popsizes + ntimes, popsizes);
        } else {
            popsizes = NULL;
        }

        // setup time steps
        setup_time_steps();
    }


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
        setup_time_steps();
    }
    ~ArgModel()
    {
        delete [] times;
        delete [] time_steps;
        if (popsizes)
            delete [] popsizes;
    }

protected:

    void setup_time_steps()
    {
        time_steps = new double [ntimes];
        for (int i=0; i<ntimes-1; i++)
            time_steps[i] = times[i+1] - times[i];
        time_steps[ntimes-1] = INFINITY;
    }

public:
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

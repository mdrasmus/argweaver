//=============================================================================
// ArgHmm model


#ifndef ARGHMM_MODEL_H
#define ARGHMM_MODEL_H

// c/c++ includes
#include <math.h>

// arghmm includes
#include "track.h"


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
    ArgModel(int ntimes, double rho, double mu) :
        ntimes(ntimes),
        times(NULL),
        time_steps(NULL),
        popsizes(NULL),
        rho(rho),
        mu(mu)
    {}

    ArgModel(int ntimes, double maxtime, double popsize, 
             double rho, double mu) :
        ntimes(ntimes),
        times(NULL),
        time_steps(NULL),
        popsizes(NULL),
        rho(rho),
        mu(mu)
    {
        set_log_times(maxtime, ntimes);
        set_popsizes(popsize, ntimes);
    }

    ArgModel(int ntimes, double maxtime, double *_popsizes, 
             double rho, double mu) :
        ntimes(ntimes),
        times(NULL),
        time_steps(NULL),
        popsizes(NULL),
        rho(rho),
        mu(mu)
    {
        set_log_times(maxtime, ntimes);
        if (_popsizes)
            set_popsizes(_popsizes, ntimes);
    }


    ArgModel(int ntimes, double *_times, double *_popsizes, 
             double rho, double mu) :
        ntimes(ntimes),
        times(NULL),
        time_steps(NULL),
        popsizes(NULL),
        rho(rho),
        mu(mu)
    {
        set_times(_times, ntimes);
        if (_popsizes)
            set_popsizes(_popsizes, ntimes);
    }
    ~ArgModel()
    {
        delete [] times;
        delete [] time_steps;
        if (popsizes)
            delete [] popsizes;
    }

    //=====================================================================
    // setting time points and popsizes

    void set_times(double *_times, int _ntimes) {
        ntimes = _ntimes;
        if (times)
            delete [] times;
        times = new double [ntimes];
        copy(_times, _times + ntimes, times);
        
        setup_time_steps();
    }

    void set_log_times(double maxtime, int _ntimes) {
        ntimes = _ntimes;
        if (times)
            delete [] times;
        times = new double [ntimes];
        get_time_points(ntimes, maxtime, times);
        setup_time_steps();
    }

    void set_linear_times(double time_step, int _ntimes) {
        ntimes = _ntimes;
        if (times)
            delete [] times;
        times = new double [ntimes];
        for (int i=0; i<ntimes; i++)
            times[i] = i * time_step;
        setup_time_steps();
    }

    void set_popsizes(double *_popsizes, int _ntimes) {
        ntimes = _ntimes;
        if (popsizes)
            delete [] popsizes;
        popsizes = new double [ntimes];
        copy(_popsizes, _popsizes + ntimes, popsizes);
    }

    void set_popsizes(double popsize, int _ntimes) {
        ntimes = _ntimes;
        if (popsizes)
            delete [] popsizes;
        popsizes = new double [ntimes];
        fill(popsizes, popsizes + ntimes, popsize);
    }

    //====================================================================
    // accessors

    bool has_mutmap()
    {
        return mutmap.size() > 0;
    }

    bool has_recombmap()
    {
        return recombmap.size() > 0;
    }

protected:

    void setup_time_steps()
    {
        if (time_steps)
            delete [] time_steps;
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
    Track<double> mutmap;
    Track<double> recombmap;
};




} // namespace arghmm


#endif // ARGHMM_MODEL_H

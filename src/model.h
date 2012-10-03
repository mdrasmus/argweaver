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
    ArgModel(int ntimes=0, double rho=0, double mu=0) :
        owned(true),
        ntimes(ntimes),
        times(NULL),
        time_steps(NULL),
        popsizes(NULL),
        rho(rho),
        mu(mu)
    {}

    // Model with constant population sizes and log-spaced time points
    ArgModel(int ntimes, double maxtime, double popsize, 
             double rho, double mu) :
        owned(true),
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

    // Model with variable population sizes and log-space time points
    ArgModel(int ntimes, double maxtime, double *_popsizes, 
             double rho, double mu) :
        owned(true),
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

        
    // Model with custom time points and variable population sizes
    ArgModel(int ntimes, double *_times, double *_popsizes, 
             double rho, double mu) :
        owned(true),
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

        
    // share data reference
    ArgModel(const ArgModel &other, double rho, double mu) :
        owned(false),
        ntimes(other.ntimes),
        times(other.times),
        time_steps(other.time_steps),
        popsizes(other.popsizes),
        rho(rho),
        mu(mu)
    {}
    

    // Copy constructor
    ArgModel(const ArgModel &other) :
        ntimes(other.ntimes),
        times(NULL),
        time_steps(NULL),
        popsizes(NULL),
        rho(other.rho),
        mu(other.mu)
    {
        copy(other);
    }


    ~ArgModel()
    {
        clear();
    }

    
    // deallocate all data
    void clear() {
        if (owned) {
            delete [] times;
            delete [] time_steps;
            if (popsizes)
                delete [] popsizes;
        }
    }

protected:    
    void clear_array(double **array) {
        if (owned && *array)
            delete [] *array;
        *array = NULL;
    }

    
public:
    // Copy parameters from another model
    void copy(const ArgModel &other)
    {
        owned = true;
        rho = other.rho;
        mu = other.mu;

        // copy popsizes and times
        set_times(other.times, ntimes);
        if (other.popsizes)
            set_popsizes(other.popsizes, ntimes);
        
        // copy maps
        mutmap.insert(mutmap.begin(), other.mutmap.begin(), other.mutmap.end());
        recombmap.insert(recombmap.begin(), 
                         other.recombmap.begin(), other.recombmap.end());
    }

    // Returns dummy time used for root of tree with internal branch removed
    int get_removed_root_time() const {
        return ntimes + 1;
    }

    //=====================================================================
    // setting time points and population sizes

    // Sets the model time points from an array
    void set_times(double *_times, int _ntimes) {
        ntimes = _ntimes;
        clear_array(&times);
        times = new double [ntimes];
        std::copy(_times, _times + ntimes, times);
        
        setup_time_steps();
    }

    // Sets the model time points linearily in log space
    void set_log_times(double maxtime, int _ntimes) {
        ntimes = _ntimes;
        clear_array(&times);
        times = new double [ntimes];
        get_time_points(ntimes, maxtime, times);
        setup_time_steps();
    }

    // Sets the model time points linearily
    void set_linear_times(double time_step, int _ntimes) {
        ntimes = _ntimes;
        clear_array(&times);
        times = new double [ntimes];
        for (int i=0; i<ntimes; i++)
            times[i] = i * time_step;
        setup_time_steps();
    }

    // Sets the model population sizes from an array
    void set_popsizes(double *_popsizes, int _ntimes) {
        ntimes = _ntimes;
        clear_array(&popsizes);
        popsizes = new double [ntimes];
        std::copy(_popsizes, _popsizes + ntimes, popsizes);
    }

    // Sets the model populations to be constant over all time points
    void set_popsizes(double popsize, int _ntimes) {
        ntimes = _ntimes;
        clear_array(&popsizes);
        popsizes = new double [ntimes];
        fill(popsizes, popsizes + ntimes, popsize);
    }

    //====================================================================
    // maps

    // Returns true if mutation map is present
    bool has_mutmap() const {
        return mutmap.size() > 0;
    }

    // Returns true if recombination map is present
    bool has_recombmap() const {
        return recombmap.size() > 0;
    }

    // Initializes mutation and recombination maps for use
    void setup_maps(string chrom, int start, int end);

    // set model parameters from map position
    void set_map_pos(int pos) {
        mu = mutmap.find(pos, mu);
        rho = recombmap.find(pos, rho);
    }

    // Returns a model customized for the local position
    void get_local_model(int pos, ArgModel &model) const {
        model.mu = mutmap.find(pos, mu);
        model.rho = recombmap.find(pos, rho);

        model.owned = false;
        model.times = times;
        model.ntimes = ntimes;
        model.time_steps = time_steps;
        model.popsizes = popsizes;
    }
    

protected:

    // Setup time steps between time points
    void setup_time_steps()
    {
        clear_array(&time_steps);
        time_steps = new double [ntimes];
        for (int i=0; i<ntimes-1; i++)
            time_steps[i] = times[i+1] - times[i];
        time_steps[ntimes-1] = INFINITY;
    }

public:
    bool owned; // if true, this object owns the array pointers

    // time points (presented in generations)
    int ntimes;
    double *times;
    double *time_steps;

    // parameters
    double *popsizes;        // population sizes
    double rho;              // recombination rate (recombs/generation/site)
    double mu;               // mutation rate (mutations/generation/site)
    Track<double> mutmap;    // mutation map
    Track<double> recombmap; // recombination map
};




} // namespace arghmm


#endif // ARGHMM_MODEL_H

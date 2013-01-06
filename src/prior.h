#ifndef ARGHMM_PRIOR_H
#define ARGHMM_PRIOR_H

// arghmm includes
#include "expm/matrix_exponential.h"
#include "local_tree.h"
#include "model.h"


namespace arghmm {


void make_tavare_matrix(int k, double n, double *matrix);
double *calc_coal_counts_matrix(int k, double t, double n);

class CoalCountsMatrix
{
public:
    CoalCountsMatrix(int maxlineages) :
        maxlineages(maxlineages),
        tavare(NULL),
        matrix(NULL)
    {        
    }

    ~CoalCountsMatrix()
    {
        if (tavare)
            delete [] tavare;
        if (matrix)
            free(matrix);
    }

    // t -- time interval
    // n -- haploid population size
    void calc(double t, double n)
    {
        if (!tavare)
            tavare = new double [maxlineages*maxlineages];
        make_tavare_matrix(maxlineages, n/t, tavare);
     
        if (matrix)
            free(matrix);
        matrix = expm11(maxlineages, tavare);
    }

    // get probability of going from a to b lineages
    inline double get(int a, int b)
    {
        return matrix[maxlineages*(a-1) + b-1];
    }

    int maxlineages;
    double *tavare;
    double *matrix;
};


extern "C" {

// The probabiluty of going from 'a' lineages to 'b' lineages in time 't'
// with population size 'n'
double prob_coal_counts_matrix(int a, int b, double t, double n);

} // extern "C"



} // namespace arghmm

#endif // ARGHMM_PRIOR_H



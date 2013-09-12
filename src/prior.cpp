
// c/c++ includes
#include "stdlib.h"

// arghmm includes
#include "prior.h"


namespace argweaver {


static inline double H(int i, int j, double n)
{
    double p = i*(i-1)/2.0/n;
    if (j == i - 1)
        return p;
    else if (j == i)
        return -p;
    else
        return 0.0;
}


void make_tavare_matrix(int k, double n, double *matrix)
{
    for (int i=0; i<k; i++)
        for (int j=0; j<k; j++)
            matrix[k*i+j] = H(i+1, j+1, n);
}


double *calc_coal_counts_matrix(int k, double t, double n)
{
    double *tavare = new double [k*k];
    make_tavare_matrix(k, n/t, tavare);
    double *matrix = expm11(k, tavare);
    delete [] tavare;
    return matrix;
}


extern "C" {

// The probability of going from 'a' lineages to 'b' lineages in time 't'
// with population size 'n'
double prob_coal_counts_matrix(int a, int b, double t, double n)
{
    double *matrix = calc_coal_counts_matrix(a, t, n);
    double p = matrix[(a-1)*a + b-1];
    free(matrix);
    return p;
}

} // extern "C"

} // namespace argweaver


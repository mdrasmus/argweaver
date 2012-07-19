#ifndef ARGHMM_COMMON_H
#define ARGHMM_COMMON_H


// headers c++ 
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <algorithm>

#include "t2exp.h"


namespace arghmm {

using namespace std;

// constants
#ifndef INFINITY
#   define INFINITY 1e1000
#endif

#ifndef t2exp
#   define t2exp exp
#endif

// indexing a matrix stored as a single array (row-major)
// m: number of columns
// i: row index
// j: column index
#define matind(m, i, j) ((m)*(i) + (j))


// find value in array
template<class T>
inline int find_array(const T* list, int size, const T value)
{
    for (int i=0; i<size; i++)
        if (list[i] == value)
            return i;
    return -1;
}


template<class ListT, class T>
inline int find_vector(const ListT &list, const T value)
{
    const unsigned int size = list.size();
    for (unsigned int i=0; i<size; i++)
        if (list[i] == value)
            return i;
    return -1;
}



//=============================================================================
// sorting

template <class KeyType, class ValueType>
struct RankSortCmp
{
    RankSortCmp(ValueType *values): values(values) {}
    
    bool operator()(KeyType i, KeyType j)
    { return values[i] < values[j]; }
    
    ValueType *values;
};

template <class KeyType, class ValueType>
void ranksort(KeyType *keys, ValueType *values, int size)
{
    RankSortCmp<KeyType, ValueType> cmp(values);
    sort(keys, keys + size, cmp);
}


void invertPerm(int *perm, int *inv, int size);

template <class T>
void permute(T* array, int *perm, int size)
{
    T *tmp = new T [size];
    
    // transfer permutation to temp array
    for (int i=0; i<size; i++)
        tmp[i] = array[perm[i]];
    
    // copy permutation back to original array
    for (int i=0; i<size; i++)
        array[i] = tmp[i];

    delete [] tmp;
}


//=============================================================================
// Math

inline double frand()
{ return rand() / double(RAND_MAX); }

inline double frand(double max)
{ return rand() / double(RAND_MAX) * max; }

inline double frand(double min, double max)
{ return min + (rand() / double(RAND_MAX) * (max-min)); }

inline int irand(int max)
{
    const int i = int(rand() / float(RAND_MAX) * max); 
    return (i == max) ? max - 1 : i;
}

inline int irand(int min, int max)
{
    const int i = min + int(rand() / float(RAND_MAX) * (max - min)); 
    return (i == max) ? max - 1 : i;
}

inline double expovariate(double lambda)
{ return -log(frand()) / lambda; }


// shuffle array
template<class T>
inline void shuffle(T *list, int size)
{
    for (int i=1; i<size; i++) {
        int j = irand(i);
        T tmp = list[i];
        list[i] = list[j];
        list[j] = tmp;
    }
}

/*
inline double logadd(double lna, double lnb)
{
    if (lna == 1.0)
         return lnb;
    if (lnb == 1.0)
        return lna;
    double diff = lna - lnb;
    if (diff < 500.0)
        return log(exp(diff) + 1.0) + lnb;
    else
         return lna;
}
*/ 


// computes log(a + b) given log(a) and log(b)
inline double logadd(double lna, double lnb) {
    if (lna == -INFINITY)
        return lnb;
    if (lnb == -INFINITY)
        return lna;
    return max(lna, lnb) + log1p(t2exp(-fabs(lna - lnb)));
}


// subtracting numbers in log-space
// NOTE: must have lna > lnb
inline double logsub(double lna, double lnb)
{
    double diff = lna - lnb;
    if (diff < 500) {
        diff = exp(diff) - 1.0;
        if (diff == 0.0)
            return -INFINITY;
        else
            return log(diff) + lnb;
    } else
        return lna;
}


inline double logsum(const double *vals, int nvals, const double threshold=-15)
{
    if (nvals == 0)
        return 1.0;

    double maxval = vals[0];
    int maxi = 0;

    // find maxval
    for (int i=1; i<nvals; i++) {
        if (vals[i] > maxval) {
            maxval = vals[i];
            maxi = i;
        }
    }
    
    // NOTE: for i = maxi, exp(vals[i] - maxval) = 1.0
    // inorder to discount for this value, we start the expsum at 0.0 
    // instead of 1.0
    double expsum = 0.0;
    for (int i=0; i<nvals; i++)
        if (vals[i] - maxval > threshold)
            expsum += t2exp(vals[i] - maxval);
  
    return maxval + log(expsum);        
}


template <class T>
T ipow(T val, int expo)
{
    T result = 1.0;
    unsigned int e = expo;

    if ((int)e < 0) {
	e = -e;
	val = 1.0 / val;
    }

    while (true) {
	if (e & 1)
	    result *= val;
	if ((e >>= 1) == 0)
	    break;
	val *= val;
    }

    return result;
}


template <class T>
T **new_matrix(int nrows, int ncols)
{
    T **mat = new T* [nrows];
    T *block = new T [nrows * ncols];
    for (int i=0; i<nrows; i++)
        mat[i] = &block[i*ncols];
    return mat;
}


template <class T>
void delete_matrix(T **mat, int nrows)
{
    delete [] mat[0];
    delete [] mat;
}


template <class T>
T **new_matrix_rows(int nrows, int ncols)
{
    T **mat = new T* [nrows];
    for (int i=0; i<nrows; i++)
        mat[i] = new T [ncols];
    return mat;
}


template <class T>
void delete_matrix_rows(T **mat, int nrows)
{
    for (int i=0; i<nrows; i++)
        delete [] mat[i];
    delete [] mat;
}



int choose(int n, int k);

double fchoose(int n, int k);


inline int sample(const double *weights, int nweights)
{
    // find total weight
    double total = 0.0;
    for (int i=0; i<nweights; i++)
        total += weights[i];

    // make random sample
    double pick = frand(total);

    // find choosen item index
    double x = 0.0;
    for (int i=0; i<nweights; i++) {
        x += weights[i];
        if (x >= pick)
            return i;
    }
    return nweights - 1;
}


template <class T>
inline T max_array(const T* lst, int size)
{
    T top = lst[0];
    for (int i=1; i<size; i++)
        if (lst[i] > top)
            top = lst[i];
    return top;
}


// test whether two floats are approximately equal
inline bool fequal(double f1, double f2, double rel=.0001, double eabs=1e-12)
{
    if (f1 == f2)
        return true;
    
    double err;
    double diff = fabs(f1 - f2);
    if (f2 == 0.0)
        err = f1;
    else if (f1 == 0)
        err = f2;
    else
        err = diff / fabs(f2);

    if (diff < eabs)
        return true;

    return err < rel;
}


} // namespace arghmm

#endif // ARGHMM_COMMON_H

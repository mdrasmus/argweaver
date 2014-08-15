#include "IntervalIterator.h"
#include "logging.h"

#include <assert.h>
#include <iterator>
#include <list>
#include <iostream>
#include <math.h>
#include <algorithm>


namespace argweaver {

double compute_mean(const vector<double> &scores) {
    double rv=0.0;
    if (scores.size() ==0) {
        printError("Error: trying to get mean with no scores\n");
    }
    for (unsigned int i=0; i < scores.size(); i++) {
        rv += scores[i];
    }
    rv /= (double)scores.size();
    return rv;
}


double compute_stdev(const vector<double> &scores, double mean) {
    double rv=0;
    if (scores.size() <= 1)
        printError("Error: trying to get stdev with %i scores\n",
                (int)scores.size());
    for (unsigned int i=0; i < scores.size(); i++) {
        double diff=scores[i]-mean;
        rv += diff*diff;
    }
    return sqrt(rv/((double)(scores.size()-1)));
}



vector<double> compute_quantiles(vector<double> scores, vector<double> q) {
    vector<double> result(q.size());
    vector<double> p(scores.size());
    unsigned int i;
    std::sort(scores.begin(), scores.end());
    for (i=0; i < scores.size(); i++)
        p[i] = (double)i/scores.size();
    for (i=0; i < q.size(); i++) {
        if (q[i] < 0 || q[i] > 1) {
            printError("Error: quantiles expects values between 0 and 1\n");
            abort();
        }
        int pos = q[i]*scores.size();
        if (pos == (int)scores.size()) pos--;
        if (fabs(q[i]-p[pos]) < 0.00001 && pos > 0)
            result[i] = (scores[pos]+scores[pos-1])/2.0;
        else result[i] = scores[pos];
    }
    return result;
}

}

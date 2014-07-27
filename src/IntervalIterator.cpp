#include "IntervalIterator.h"

#include <assert.h>
#include <iterator>
#include <list>
#include <iostream>
#include <math.h>
#include <algorithm>


namespace argweaver {

double compute_mean(vector<double> scores) {
    double rv=0.0;
    if (scores.size() ==0) {
        fprintf(stderr, "Error: trying to get mean with no scores\n");
    }
    for (unsigned int i=0; i < scores.size(); i++) {
        rv += scores[i];
    }
    rv /= (double)scores.size();
    return rv;
}

    /*
template<class scoreT> scoreT Interval<scoreT>::mean() {
    meanval = mean(scores);
    have_mean = 1;
    return meanval;
    }*/

double compute_stdev(vector<double> scores, double mean) {
    double rv=0;
    if (scores.size() <= 1)
        fprintf(stderr, "Error: trying to get stdev with %i scores\n",
                (int)scores.size());
    for (unsigned int i=0; i < scores.size(); i++) {
        double diff=scores[i]-mean;
        rv += diff*diff;
    }
    return sqrt(rv/((double)(scores.size()-1)));
}

    /*template<class scoreT> scoreT Interval<scoreT>::stdev() {
    if (!have_mean)
        this->mean();
    return stdev(scores, meanval);
    }*/


vector<double> compute_quantiles(vector<double> scores, vector<double> q) {
    vector<double> result(q.size());
    vector<double> p(scores.size());
    unsigned int i;
    std::sort(scores.begin(), scores.end());
    for (i=0; i < scores.size(); i++)
        p[i] = (double)i/scores.size();
    for (i=0; i < q.size(); i++) {
        if (q[i] < 0 || q[i] > 1) {
            fprintf(stderr, "Error: quantiles expects values between 0 and 1\n");
            exit(1);
        }
        int pos = q[i]*scores.size();
        if (pos == (int)scores.size()) pos--;
        if (fabs(q[i]-p[pos]) < 0.00001 && pos > 0)
            result[i] = (scores[pos]+scores[pos-1])/2.0;
        else result[i] = scores[pos];
    }
    return result;
}
    /*
template<class scoreT> vector<scoreT> Interval<scoreT>::quantiles(vector<double> q) {
    return quantiles(scores, q);
    }*/


    /*Interval<scoreT> IntervalIterator<scoreT>::next() {
    Interval<scoreT> rv("", -1, -1);
    if (combined.size() > 0) {
        rv = combined.front();
        combined.pop_front();
    }
    return rv;
    }*/

    /*template<class scoreT>
void IntervalIterator<scoreT>::pushNext(string chr, int start, int end) {
    Interval<scoreT> newCombined(chr, start, end);
    typename std::list<Interval<scoreT> >::iterator curr_it, next_it;
    curr_it = intervals.begin();
    next_it = intervals.begin();
    next_it++;
    while (curr_it != intervals.end()  &&
           curr_it->chrom == chr && curr_it->start == start) {
        newCombined.add_score(curr_it->get_score(0));
        if (curr_it->end < end) {
            fprintf(stderr, "Error\n");
        }
        assert(curr_it->end >= end);
        if (curr_it->end == end) {
            intervals.erase(curr_it);
        } else {
            curr_it->start = end;
        }
        if (next_it == intervals.end()) break;
        curr_it = next_it;
        next_it++;
    }
    combined.push_back(newCombined);
    }*/


    /*template<class scoreT> void IntervalIterator<scoreT>::finish() {
    std::set<int>::iterator it=bounds.begin();
    int startCoord, endCoord;

    if (bounds.size() ==0) return;
    startCoord = *it;
    for (++it; it != bounds.end(); it++) {
        endCoord = *it;
        pushNext(chrom, startCoord, endCoord);
        startCoord = endCoord;
    }
    bounds.clear();
    }*/

    /*template<class scoreT> void IntervalIterator<scoreT>::append(string chr, int start, int end,
                                      scoreT score) {
    Interval<scoreT> newint(chr, start, end, score);
    std::set<int>::iterator it, prev_it;
    int startCoord, endCoord;

    if (intervals.size() > 0 && chrom != chr) {
        this->finish();
    }
    chrom = chr;
    bounds.insert(start);
    bounds.insert(end);
    intervals.push_back(newint);

    it = prev_it = bounds.begin();
    startCoord = *it;
    for (++it; it != bounds.end(); it++) {
        endCoord = *it;
        if (endCoord >= start) break;
        pushNext(chrom, startCoord, endCoord);
        bounds.erase(prev_it);
        prev_it = it;
        startCoord = endCoord;
    }
    }*/

    /*void IntervalIterator<scoreT>::printResults() {
    Interval<scoreT> summary = this->next();
    while (summary.start != summary.end) {
        cout << summary.chrom << "\t" << summary.start << "\t"
             << summary.end;
        for (int i=1; i <= summarize; i++) {
            if (getNumSample==i) {
                printf("\t%i", summary.num_score());
            } else if (getMean==i) {
                printf("\t%f", summary.mean());
            } else if (getStdev==i) {
                printf("\t%f", summary.stdev());
            } else if (getQuantiles==i) {
                vector<double> q = summary.quantiles(quantiles);
                for (unsigned int j=0; j < quantiles.size(); j++) {
                    printf("\t%f", q[j]);
                }
            }
        }
        cout << "\n";
        summary = results->next();
    }
    }*/

}

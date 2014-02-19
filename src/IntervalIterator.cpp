#include "IntervalIterator.h"

#include <assert.h>
#include <iterator>
#include <list>
#include <iostream>
#include <math.h>
#include <algorithm>

using namespace std;
namespace spidir {

double Interval::mean() {
    if (scores.size() ==0) {
        fprintf(stderr, "Error: trying to get mean with no scores\n");
    }
    meanval=0;
    for (unsigned int i=0; i < scores.size(); i++) {
        meanval += scores[i];
    }
    meanval /= (double)scores.size();
    have_mean=1;
    return meanval;
}

double Interval::stdev() {
    double stdev=0;
    if (scores.size() <= 1)
        fprintf(stderr, "Error: trying to get stdev with %i scores\n",
                (int)scores.size());
    if (!have_mean)
        this->mean();
    for (unsigned int i=0; i < scores.size(); i++) {
        double diff=scores[i]-meanval;
        stdev += diff*diff;
    }
    return sqrt(stdev/((double)(scores.size()-1)));
}
    /*
vector<double> Interval::quantiles(vector<double> q) {
    vector<double> result(q.size());
    vector<double> p(scores.size());
    unsigned int i;
    double dist = 1.0/scores.size();
    std::sort(scores.begin(), scores.end());
    for (i=0; i < scores.size(); i++)
        p[i] = (0.5+i)/scores.size();
    for (i=0; i < q.size(); i++) {
        if (q[i] < 0 || q[i] > 1) {
            fprintf(stderr, "Error: quantiles expects values between 0 and 1\n");
            exit(1);
        }
        if (q[i] <= p[0])
            result[i]=scores[0];
        else if (q[i] >= p[scores.size()-1])
            result[i]=scores[scores.size()-1];
        else {
            int pos1=(int)(q[i]*scores.size());
            if (q[i] < p[pos1]) pos1--;
            int pos2=pos1+1;
            if (q[i] < p[pos1] ||
                q[i] > p[pos2]) {
                fprintf(stderr, "Error: %f %f %f\n", p[pos1], q[i], p[pos2]);
                assert(q[i] >= p[pos1] && q[i] <= p[pos2]);
            }
            double d1 = 1.0-(q[i]-p[pos1])/dist;
            double d2 = 1.0-(p[pos2]-q[i])/dist;
            result[i] = d1*scores[pos1] + d2*scores[pos2];
        }
    }
    return result;
    }*/


vector<double> Interval::quantiles(vector<double> q) {
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

double Interval::quantile(double q) {
    vector<double> qv(1);
    qv[0]=q;
    return this->quantiles(qv)[0];
}

Interval IntervalIterator::next() {
    Interval rv("", -1, -1);
    if (combined.size() > 0) {
        rv = combined.front();
        combined.pop_front();
    }
    return rv;
}

void IntervalIterator::pushNext(string chr, int start, int end) {
    Interval newCombined(chr, start, end);
    std::list<Interval>::iterator curr_it, next_it;
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
}


void IntervalIterator::finish() {
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
}

void IntervalIterator::append(string chr, int start, int end, double score) {
    Interval newint(chr, start, end, score);
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
}
}

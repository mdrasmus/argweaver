#ifndef ARGWEAVER_INTERVALS_H
#define ARGWEAVER_INTERVALS_H

#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <list>
#include <vector>
#include <set>
#include <iterator>
#include <assert.h>

using namespace std;
namespace spidir {

class Interval {
 public:
 Interval(string chrom, int start, int end):
    chrom(chrom), start(start), end(end), have_mean(0)
    {
        scores.clear();
    }
 Interval(string chrom, int start, int end, double score):
    chrom(chrom), start(start), end(end), have_mean(1), meanval(score)
    {
        scores.clear();
        scores.push_back(score);
    }
    void add_score(double score) {
        scores.push_back(score);
        have_mean=0;
    }
    int num_score() {
        return scores.size();
    }
    double get_score(int i) {
        if (i < 0 || i >= (int)scores.size()) {
            fprintf(stderr, "Error in get_score; index out of range\n");
            exit(1);
        }
        return scores[i];
    }

    string chrom;
    int start;
    int end;
    double mean();
    double stdev();
    double quantile(double q);
    vector<double> quantiles(vector<double> q);
 private:
    int have_mean;
    double meanval;
    vector<double> scores;
};

class IntervalIterator
{
 public:
 IntervalIterator()
    {
        intervals.clear();
        combined.clear();
        bounds.clear();
        chrom="";
    }
 ~IntervalIterator()
     {
     }

 Interval next();

 void append(string chr, int start, int end, double score);
 void finish();

 private:
 void pushNext(string chr, int start, int end);
 list<Interval> intervals;
 list<Interval> combined;
 set<int> bounds;
 string chrom;
};
}
#endif  //ARGWEAVER_INTERVALS_H

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

namespace argweaver {

using namespace std;

double compute_mean(vector<double>scores);
double compute_stdev(vector<double> scores, double mean);
vector<double> compute_quantiles(vector<double> scores, vector<double> q);


template <class scoreT>
class Interval {
public:
    Interval(string chrom, int start, int end):
        chrom(chrom), start(start), end(end), have_mean(false)
    {
        scores.clear();
    }
    Interval(string chrom, int start, int end, scoreT score):
        chrom(chrom), start(start), end(end), have_mean(true), meanval(score)
    {
        scores.clear();
        scores.push_back(score);
    }
    void add_score(scoreT score) {
        scores.push_back(score);
        have_mean = false;
    }
    int num_score() {
        return scores.size();
    }
    scoreT get_score(int i) {
        if (i < 0 || i >= (int)scores.size()) {
            fprintf(stderr, "Error in get_score; index out of range\n");
            exit(1);
        }
        return scores[i];
    }
    vector<scoreT> get_scores() {
        return scores;
    }

    scoreT mean() {
        meanval = compute_mean(scores);
        have_mean = true;
        return meanval;
    }
    scoreT stdev() {
        if (!have_mean)
            this->mean();
        return compute_stdev(scores, meanval);
    }
    vector<scoreT> quantiles(vector<double> q) {
        return compute_quantiles(scores, q);
    }

    string chrom;
    int start;
    int end;

private:
    bool have_mean;
    scoreT meanval;
    vector<scoreT> scores;
};


template <class scoreT>
class IntervalIterator
{
public:
    IntervalIterator()
    {
        intervals.clear();
        combined.clear();
        bounds.clear();
        chrom = "";
    }
    ~IntervalIterator()
    {
    }

    Interval<scoreT> next() {
        Interval<scoreT> rv("", -1, -1);
        if (combined.size() > 0) {
            rv = combined.front();
            combined.pop_front();
        }
        return rv;
    }

    void append(string chr, int start, int end, scoreT score) {
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
    }

    void finish() {
        std::set<int>::iterator it=bounds.begin();
        int startCoord, endCoord;

        if (bounds.size() == 0) return;
        startCoord = *it;
        for (++it; it != bounds.end(); it++) {
            endCoord = *it;
            pushNext(chrom, startCoord, endCoord);
            startCoord = endCoord;
        }
        bounds.clear();
    }
    // void printResults();

private:
    void pushNext(string chr, int start, int end) {
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
    }

    list<Interval<scoreT> > intervals;
    list<Interval<scoreT> > combined;
    set<int> bounds;
    string chrom;
};

} // namespace argweaver
#endif  //ARGWEAVER_INTERVALS_H

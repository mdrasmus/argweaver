
#include "model.h"



namespace argweaver {


void get_coal_time_steps(const double *times, int ntimes,
                          double *coal_time_steps)
{
    // get midpoints
    double times2[2*ntimes+1];
    for (int i=0; i<ntimes-1; i++) {
        times2[2*i] = times[i];
        times2[2*i+1] = sqrt((times[i+1]+1.0)*(times[i]+1.0));
    }
    times2[2*ntimes] = times[ntimes-1];

    for (int i=0; i<2*ntimes-2; i++)
        coal_time_steps[i] = times2[min(i+1, 2*ntimes)] - times2[i];
    coal_time_steps[2*ntimes-2] = INFINITY;
}



// returns true if regions in track are flush with one another
template <class T>
bool check_map(const Track<T> &track, int start, int end)
{
    // check that track is not empty
    typename Track<T>::const_iterator it = track.begin();
    if (it == track.end()) {
        printError("map is empty");
        return false;
    }

    // check that start and end cover desired range
    if (it->start > start || track.back().end < end) {
        printError("map does not cover entire region");
        return false;
    }

    int last = it->end;
    ++it;
    for (; it != track.end(); ++it) {
        if (it->start != last) {
            printError("map is not complete at %s:%d",
                       it->chrom.c_str(), it->start);
            return false;
        }
        last = it->end;
    }

    return true;
}


// returns true if regions in track are flush with one another
template <class T>
bool complete_map(Track<T> &track, string chrom, int start, int end, const T &default_value)
{
    // check for empty track
    if (track.size() == 0) {
        track.append(chrom, start, end, default_value);
        return true;
    }

    // check that start and end cover desired range
    typename Track<T>::iterator it = track.begin();
    if (it->start > start)
        track.insert(it, RegionValue<T>(chrom, start, it->start, default_value));
    if (track.back().end < end)
        track.append(chrom, track.back().end, end, default_value);

    it = track.begin();
    int last = it->end;
    ++it;
    for (; it != track.end(); ++it) {
        if (it->start > last) {
            it = track.insert(it, RegionValue<T>(chrom, last, it->start, default_value));
        } else if (it->start < last) {
            printError("map contains over laps %s:%d-%d",
                       chrom.c_str(), it->start, last);
            return false;
        }
        last = it->end;
    }

    return true;
}



// Initializes mutation and recombination maps for use
bool ArgModel::setup_maps(string chrom, int start, int end) {

    // check maps
    if (!complete_map(mutmap, chrom, start, end, mu)) {
        printError("mutation map has errors");
        return false;
    }
    if (!complete_map(recombmap, chrom, start, end, rho)) {
        printError("recombination map has errors");
        return false;
    }

    /*
    // setup default maps
    if (mutmap.size() == 0)
        mutmap.append(chrom, start, end, mu);
    if (recombmap.size() == 0)
        recombmap.append(chrom, start, end, rho);

    // check maps
    if (!check_map(mutmap, start, end)) {
        printError("mutation map has errors");
        return false;
    }
    if (!check_map(recombmap, start, end)) {
        printError("recombination map has errors");
        return false;
    }
    */

    // create new mut and recomb maps that share common boundaries
    int pos = start, pos2;
    unsigned int i = 0;
    unsigned int j = 0;
    Track<double> mutmap2;
    Track<double> recombmap2;
    while (i < mutmap.size() || j < recombmap.size()) {
        //printf("recomb[%d] = (%d, %d, %e), mut[%d] = (%d, %d, %e)\n",
        //       j, recombmap[j].start, recombmap[j].end, recombmap[j].value,
        //       i, mutmap[i].start, mutmap[i].end, mutmap[i].value);

        if (mutmap[i].end < recombmap[j].end) {
            // mutation region ends first
            pos2 = mutmap[i].end;
            mutmap2.append(chrom, pos, pos2, mutmap[i].value);
            recombmap2.append(chrom, pos, pos2, recombmap[j].value);
            pos = pos2;
            i++;
        } else if (mutmap[i].end > recombmap[j].end) {
            // recombination region ends first
            pos2 = recombmap[j].end;
            mutmap2.append(chrom, pos, pos2, mutmap[i].value);
            recombmap2.append(chrom, pos, pos2, recombmap[j].value);
            pos = pos2;
            j++;
        } else {
            // mutation and region region end together
            pos2 = recombmap[j].end;
            mutmap2.append(chrom, pos, pos2, mutmap[i].value);
            recombmap2.append(chrom, pos, pos2, recombmap[j].value);
            pos = pos2;
            i++;
            j++;
        }
    }

    // copy over new maps
    mutmap.clear();
    recombmap.clear();
    mutmap.insert(mutmap.begin(), mutmap2.begin(), mutmap2.end());
    recombmap.insert(recombmap.begin(), recombmap2.begin(), recombmap2.end());

    return true;
}


} // namespace argweaver


#include "model.h"



namespace arghmm {


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



// Initializes mutation and recombination maps for use
void ArgModel::setup_maps(string chrom, int start, int end) {
    // setup default maps
    if (mutmap.size() == 0)
        mutmap.append(chrom, start, end, mu);
    if (recombmap.size() == 0)
        recombmap.append(chrom, start, end, rho);

    // create new mut and recomb maps that share common boundaries
    int pos = start, pos2;
    unsigned int i = 0;
    unsigned int j = 0;
    Track<double> mutmap2;
    Track<double> recombmap2;
    while (i < mutmap.size() || j < recombmap.size()) {
        printf("recomb[%d] = (%d, %d, %e), mut[%d] = (%d, %d, %e)\n", 
               j, recombmap[j].start, recombmap[j].end, recombmap[j].value,
               i, mutmap[i].start, mutmap[i].end, mutmap[i].value);

        if (mutmap[i].end < recombmap[j].end) {
            pos2 = mutmap[i].end;
            mutmap2.append(chrom, pos, pos2, mutmap[i].value);
            recombmap2.append(chrom, pos, pos2, recombmap[j].value);
            pos = pos2;
            i++;
        } else if (mutmap[i].end > recombmap[j].end) {
            pos2 = recombmap[j].end;
            mutmap2.append(chrom, pos, pos2, mutmap[i].value);
            recombmap2.append(chrom, pos, pos2, recombmap[j].value);
            pos = pos2;
            j++;
        } else {
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
}


} // namespace arghmm

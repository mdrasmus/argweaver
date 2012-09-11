
#include "model.h"



namespace arghmm {

void ArgModel::setup_maps(string chrom, int start, int end) {
    // setup default maps
    if (mutmap.size() == 0)
        mutmap.append(chrom, start, end, mu);
    if (recombmap.size() == 0)
        recombmap.append(chrom, start, end, rho);

    int pos = start, pos2;
    unsigned int i = 0;
    unsigned int j = 0;
    Track<double> mutmap2;
    Track<double> recombmap2;
    while (i < mutmap.size() || j < recombmap.size()) {
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

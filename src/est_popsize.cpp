
#include <algorithm>
#include <vector>

#include "est_popsize.h"

namespace argweaver {

using namespace std;


void est_popsize_arg(const ArgModel *model, const LocalTrees *trees,
                     double *popsizes)
{

}


void est_popsize_trees2(const ArgModel *model, const LocalTree *const *trees,
                       int ntrees, double *popsizes)
{
    assert(ntrees > 0);

    const int ntimes = model->ntimes;
    const int nleaves = trees[0]->get_num_leaves();
    LineageCounts lineages(ntimes);

    int total_ncoals[ntimes];
    int total_pairs[ntimes];
    int total_ncoals_pairs[ntimes];

    fill(total_ncoals, total_ncoals + ntimes, 0);
    fill(total_pairs, total_pairs + ntimes, 0);
    fill(total_ncoals_pairs, total_ncoals_pairs + ntimes, 0);

    printf("ntrees %d\n", ntrees);

    // count lineages
    for (int i=0; i<ntrees; i++) {
        lineages.count(trees[i]);

        for (int j=0; j<ntimes-1; j++) {
            int start = (j == 0) ? nleaves : lineages.nbranches[j-1];
            int end = lineages.nbranches[j];
            int ncoals = start - end;
            int pairs = start * (start-1) / 2;

            total_ncoals[j] += ncoals;
            total_pairs[j] += pairs;
            total_ncoals_pairs[j] += ncoals * pairs;
        }
    }

    for (int j=0; j<ntimes-1; j++) {
        if (total_ncoals[j] == 0)
            popsizes[j] = 0.0;
        else
            popsizes[j] = .5 * model->time_steps[j] * (
                total_ncoals_pairs[j] + total_pairs[j] - total_ncoals[j]) /
                double(total_ncoals[j]);
        printf("> %d %d %d\n", total_ncoals_pairs[j], total_pairs[j],
               total_ncoals[j]);
        printf("popsize %f\n", popsizes[j]);
    }
}



void est_popsize_trees(const ArgModel *model, const LocalTree *const *trees,
                       int ntrees, double *popsizes)
{
    assert(ntrees > 0);

    const int ntimes = model->ntimes;
    const int nleaves = trees[0]->get_num_leaves();
    LineageCounts lineages(ntimes);

    int total_ncoals[ntimes];
    int total_pairs[ntimes];

    fill(total_ncoals, total_ncoals + ntimes, 0);
    fill(total_pairs, total_pairs + ntimes, 0);

    printf("ntrees %d\n", ntrees);

    // count lineages
    for (int i=0; i<ntrees; i++) {
        lineages.count(trees[i]);

        for (int j=0; j<ntimes-1; j++) {
            int start = (j == 0) ? nleaves : lineages.nbranches[j-1];
            int end = lineages.nbranches[j];
            int ncoals = start - end;
            int pairs = start * (start-1) / 2;

            total_ncoals[j] += ncoals;
            total_pairs[j] += pairs;
        }
    }

    for (int j=0; j<ntimes-1; j++) {
        if (total_ncoals[j] == 0)
            popsizes[j] = 0.0;
        else
            popsizes[j] = .5 * model->time_steps[j] * total_pairs[j] /
                double(total_ncoals[j]);
        printf("> %d %d\n", total_pairs[j], total_ncoals[j]);
        printf("popsize %f\n", popsizes[j]);
    }
}


void est_popsize_trees(const ArgModel *model, const LocalTrees *trees,
                       int step, double *popsizes)
{
    vector<LocalTree*> indep_trees;

    int end = trees->start_coord;
    int pos = end;
    for (LocalTrees::const_iterator it=trees->begin();
         it != trees->end(); ++it)
    {
        int start = end;
        end += it->blocklen;

        while (start <= pos && pos < end) {
            // record tree
            indep_trees.push_back(it->tree);
            pos += step;
        }
    }

    est_popsize_trees(model, &indep_trees[0], indep_trees.size(), popsizes);
}


//=============================================================================
// C-interface

extern "C" {

// estimate population sizes
void arghmm_est_popsizes_trees(LocalTrees *trees, double *times, int ntimes,
                               int step, double *popsizes)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, NULL, 0.0, 0.0);
    est_popsize_trees(&model, trees, step, popsizes);
}


} // extern C

} // namespace argweaver


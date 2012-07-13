// c++ includes
#include <list>
#include <vector>
#include <string.h>

// arghmm includes
#include "common.h"
#include "emit.h"
#include "local_tree.h"
#include "sequences.h"


using namespace std;

namespace arghmm {


// calculate the probability of the sequences given an ARG
double calc_arg_likelihood(ArgModel *model, Sequences *sequences, 
                           LocalTrees *trees)
{
    double lnl = 0.0;
    int nseqs = sequences->get_nseqs();
    double minlen = model->times[1];
    const double log25 = -1.3862943611198906; // log(.25)
    const double *times = model->times;

    if (trees->nnodes < 3)
        return lnl += log25 * sequences->length();

    // get sequences for trees
    char *seqs[nseqs];
    for (int j=0; j<nseqs; j++)
        seqs[j] = sequences->seqs[trees->seqids[j]];

    int end = trees->start_coord;
    for (LocalTrees::iterator it=trees->begin(); it != trees->end(); ++it) {
        int start = end;
        end = start + it->blocklen;
        LocalTree *tree = it->tree;
        const LocalNode *nodes = tree->nodes;
        const int blocklen = end - start;

        int root = tree->root;
        int root1 = nodes[root].child[0];
        int root2 = nodes[root].child[1];

        // initialize per branch substitution counts
        int counts[tree->nnodes];
        for (int j=0; j<tree->nnodes; j++)
            counts[j] = 0;

        // process block
        for (int i=start; i<end; i++) {
            char ancestral[tree->nnodes];
            parsimony_ancestral_seq(tree, seqs, nseqs, i, ancestral);
            
            // count substituions per branch
            for (int j=0; j<tree->nnodes; j++) {
                if (j == root || j == root2)
                    continue;
                if (j == root1 && 
                    ancestral[root1] != ancestral[root2])
                    counts[j]++;
                else if (ancestral[j] != ancestral[nodes[j].parent])
                    counts[j]++;
            }
        }

        // prior probability of root sequence        
        lnl += log25 * blocklen; 
        
        // compute likelihood of seeing each of those substituions
        for (int j=0; j<tree->nnodes; j++) {
            if (j == root || j == root2)
                continue;
            
            int parent = nodes[j].parent;
            double t;
            if (j == root1)
                // wrap branch
                t = 2.0 * times[nodes[root].age] 
                    - times[nodes[root1].age] 
                    - times[nodes[root2].age];
            else
                t = times[nodes[parent].age] - times[nodes[j].age];
            t = max(t, minlen);
            lnl += (counts[j] - blocklen) * model->mu * t +
                counts[j] * log(1/3. - 1/3. * exp(-model->mu * t));
        }
    }

    return lnl;
}



double calc_spr_prob(const ArgModel *model, const LocalTree *tree, 
                     const Spr &spr, LineageCounts &lineages)
{
    double lnl = 0.0;

    const LocalNode *nodes = tree->nodes;
    const int root_age = nodes[tree->root].age;

    // get tree length
    const double treelen = get_treelen(
        tree, model->times, model->ntimes, false);
    const double treelen_b = treelen + model->time_steps[nodes[tree->root].age];

    // get lineage counts
    lineages.count(tree);
    lineages.nrecombs[root_age]--;

    assert(spr.recomb_node != tree->root);
            
    // probability of recombination location in tree
    int k = spr.recomb_time;
    lnl += log(lineages.nbranches[k] * model->time_steps[k] /
               (lineages.nrecombs[k] * treelen_b));

    // probability of re-coalescence
    int j = spr.coal_time;
    int broken_age = nodes[nodes[spr.recomb_node].parent].age;
    int ncoals_j = lineages.ncoals[j] 
        - int(j <= broken_age) - int(j == broken_age);
    int nbranches_j = lineages.nbranches[j] - int(j < broken_age);

    lnl -= log(ncoals_j);
    if (j < model->ntimes - 2)
        lnl += log((1.0 - exp(- model->time_steps[j] * nbranches_j / 
                              (2.0 * model->popsizes[j]))));

    double sum = 0.0;
    for (int m=k; m<j; m++) {
        int nbranches_m = lineages.nbranches[m] - int(m < broken_age);
        
        sum += model->time_steps[m] * nbranches_m / 
            (2.0 * model->popsizes[m]);
    }
    lnl -= sum;
    
    return lnl;
}



// calculate the probability of an ARG given the model parameters
double calc_arg_prior(ArgModel *model, LocalTrees *trees)
{
    double lnl = 0.0;
    LineageCounts lineages(model->ntimes);

    // TODO: add tree prior

    int end = trees->start_coord;
    for (LocalTrees::iterator it=trees->begin(); it != trees->end();) {
        end += it->blocklen;
        LocalTree *tree = it->tree;
        int blocklen = it->blocklen;
        double treelen = get_treelen(tree, model->times, model->ntimes, false);

        // calculate probability P(blocklen | T_{i-1})
        double recomb_rate = max(model->rho * treelen, model->rho);

        if (end < trees->end_coord) {
            // not last block
            // probability of recombining after blocklen
            lnl += log(recomb_rate) - recomb_rate * blocklen;
            
            // get SPR move information
            ++it;
            Spr *spr = &it->spr;
            lnl += calc_spr_prob(model, tree, *spr, lineages);

        } else {
            // last block
            // probability of not recombining after blocklen
            lnl += - recomb_rate * blocklen;
            ++it;
        }
    }

    return lnl;    
}


// calculate the probability of the sequences given an ARG
double calc_arg_joint_prob(ArgModel *model, Sequences *sequences, 
                           LocalTrees *trees)
{
    return calc_arg_likelihood(model, sequences, trees) +
        calc_arg_prior(model, trees);
}



//=============================================================================
// C interface

extern "C" {

double arghmm_likelihood(LocalTrees *trees,
                         double *times, int ntimes,
                         double mu, 
                         char **seqs, int nseqs, int seqlen)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, NULL, 0.0, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    return calc_arg_likelihood(&model, &sequences, trees);
}


double arghmm_prior_prob(LocalTrees *trees,
                         double *times, int ntimes, double *popsizes,
                         double rho)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, 0.0);
    return calc_arg_prior(&model, trees);
}


double arghmm_joint_prob(LocalTrees *trees,
                         double *times, int ntimes, double *popsizes,
                         double mu, double rho,
                         char **seqs, int nseqs, int seqlen)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    return calc_arg_joint_prob(&model, &sequences, trees);
}

/*

double arghmm_likelihood(
    int **ptrees, int **ages, int **sprs, int *blocklens,
    int ntrees, int nnodes, 
    double *times, int ntimes,
    double mu, 
    char **seqs, int nseqs, int seqlen)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, NULL, 0.0, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    LocalTrees trees(ptrees, ages, sprs, blocklens, ntrees, nnodes);

    return calc_arg_likelihood(&model, &sequences, &trees);
}


double arghmm_prior_prob(
    int **ptrees, int **ages, int **sprs, int *blocklens,
    int ntrees, int nnodes, 
    double *times, int ntimes, double *popsizes,
    double rho)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, 0.0);
    LocalTrees trees(ptrees, ages, sprs, blocklens, ntrees, nnodes);

    return calc_arg_prior(&model, &trees);
}


double arghmm_joint_prob(
    int **ptrees, int **ages, int **sprs, int *blocklens,
    int ntrees, int nnodes, 
    double *times, int ntimes, double *popsizes,
    double mu, double rho,
    char **seqs, int nseqs, int seqlen)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    LocalTrees trees(ptrees, ages, sprs, blocklens, ntrees, nnodes);

    return calc_arg_joint_prob(&model, &sequences, &trees);
}

*/


} // extern "C"


} // namespace arghmm

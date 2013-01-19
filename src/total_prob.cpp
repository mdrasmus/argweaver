// c++ includes
#include <list>
#include <vector>
#include <string.h>

// arghmm includes
#include "common.h"
#include "emit.h"
#include "local_tree.h"
#include "sequences.h"
#include "trans.h"


using namespace std;

namespace arghmm {


double calc_arg_likelihood(const ArgModel *model, const Sequences *sequences, 
                           const LocalTrees *trees)
{
    double lnl = 0.0;
    int nseqs = sequences->get_num_seqs();
    
    // special case for truck genealogies
    if (trees->nnodes < 3)
        return lnl += log(.25) * sequences->length();
    
    // get sequences for trees
    char *seqs[nseqs];
    for (int j=0; j<nseqs; j++)
        seqs[j] = sequences->seqs[trees->seqids[j]];

    int end = trees->start_coord;
    for (LocalTrees::const_iterator it=trees->begin(); it!=trees->end(); ++it) {
        int start = end;
        end = start + it->blocklen;
        LocalTree *tree = it->tree;

        lnl += likelihood_tree(tree, model, seqs, nseqs, start, end);
    }

    return lnl;
}


// NOTE: trees should be uncompressed and sequences compressed
double calc_arg_likelihood(const ArgModel *model, const Sequences *sequences, 
                           const LocalTrees *trees, 
                           const SitesMapping* sites_mapping)
{
    if (!sites_mapping)
        return calc_arg_likelihood(model, sequences, trees);

    double lnl = 0.0;
    int nseqs = sequences->get_num_seqs();
    const char default_char = 'A';
    
    // special case for truck genealogies
    if (trees->nnodes < 3)
        return lnl += log(.25) * sequences->length();
    
    int end = trees->start_coord;
    for (LocalTrees::const_iterator it=trees->begin(); it!=trees->end(); ++it) {
        int start = end;
        int blocklen = it->blocklen;
        end = start + blocklen;
        LocalTree *tree = it->tree;

        // get sequences for trees
        char *seqs[nseqs];
        char *matrix = new char [blocklen*nseqs];
        for (int j=0; j<nseqs; j++)
            seqs[j] = &matrix[j*blocklen];

        // find first site within this block
        unsigned int i2 = 0;
        
        // copy sites into new alignment
        for (int i=start; i<end; i++) {
            while (i2 < sites_mapping->all_sites.size() && 
                   sites_mapping->all_sites[i2] < i)
                i2++;
            if (i == sites_mapping->all_sites[i2]) {
                // copy site
                for (int j=0; j<nseqs; j++)
                    seqs[j][i-start] = sequences->seqs[trees->seqids[j]][i2];
            } else {
                // copy non-variant site
                for (int j=0; j<nseqs; j++)
                    seqs[j][i-start] = default_char;
            }
        }

        lnl += likelihood_tree(tree, model, seqs, nseqs, 0, end-start);

        delete [] matrix;
    }

    return lnl;
}

//=============================================================================
// ARG prior


// The probabiluty of going from 'a' lineages to 'b' lineages in time 't'
// with population size 'n'
double prob_coal_counts(int a, int b, double t, double n)
{
    double C = 1.0;
    
    for (int y=0; y<b; y++)
        C *= (b+y)*(a-y)/double(a+y);

    double s = exp(-b*(b-1)*t/2.0/n) * C;

    for (int k=b+1; k<a+1; k++) {
        const double k1 = double(k - 1);
        C *= double(b+k1)*(a-k1)/(a+k1)/(b-k);
        s += exp(-k*k1*t/2.0/n) * (2*k-1) / double(k1+b) * C;
    }
    
    for (int i=1; i<=b; i++)
        s /= i;

    return s;
}


double calc_tree_prior(const ArgModel *model, const LocalTree *tree,
                       LineageCounts &lineages)
{
    lineages.count(tree);
    int nleaves = tree->get_num_leaves();
    double lnl = 0.0;

    // get effective population sizes
    // uses harmonic mean to combine time steps
    double times[model->ntimes];
    double popsizes[model->ntimes];
    times[0] = model->coal_time_steps[0];
    popsizes[0] = model->popsizes[0];

    for (int i=1; i<model->ntimes-1; i++) {
        double t1 = model->coal_time_steps[2*i-1];
        double t2 = model->coal_time_steps[2*i];
        double n1 = model->popsizes[i-1];
        double n2 = model->popsizes[i];

        times[i] = t1 + t2;
        popsizes[i] = (t1 + t2) / (t1/n1 + t2/n2);
    }
    
    for (int i=0; i<model->ntimes-1; i++) {
        int a = (i == 0 ? nleaves : lineages.nbranches[i-1]);
        int b = lineages.nbranches[i];
        lnl += log(prob_coal_counts(a, b, times[i], 2.0*popsizes[i]));
    }

    // TODO: top prior
    
    return lnl;
}


double calc_tree_prior_approx(const ArgModel *model, const LocalTree *tree,
                              LineageCounts &lineages)
{
    lineages.count(tree);
    double lnl = 0.0;
    
    for (int i=0; i<tree->nnodes; i++) {
        if (tree->nodes[i].is_leaf())
            continue;
        int time = tree->nodes[i].age;

        // remove lineage counts
        for (int j=0; j<time; j++)
            lineages.nbranches[j]--;

        lnl += log(calc_state_priors(time, lineages.nbranches, lineages.ncoals,
                                     0, model->popsizes, model->coal_time_steps,
                                     model->ntimes));
        lnl += log(lineages.ncoals[time]);
    }

    if (isnan(lnl))
        lnl = 0.0;
    
    return lnl;
}


void calc_coal_rates_full_tree(const ArgModel *model, const LocalTree *tree, 
                      const Spr &spr, LineageCounts &lineages,
                      double *coal_rates)
{
    int broken_age = tree->nodes[tree->nodes[spr.recomb_node].parent].age;

    for (int i=0; i<2*model->ntimes; i++) {
        int nbranches = lineages.nbranches[i/2] - int(i/2 < broken_age);
        coal_rates[i] = model->coal_time_steps[i] * nbranches / 
            (2.0 * model->popsizes[i/2]);
    }
}


double calc_spr_prob(const ArgModel *model, const LocalTree *tree, 
                     const Spr &spr, LineageCounts &lineages, 
                     double treelen)
{
    assert(spr.recomb_node != tree->root);
    const LocalNode *nodes = tree->nodes;
    const int root_age = nodes[tree->root].age;

    // get tree length, if it is not already given
    if (treelen < 0)
        treelen = get_treelen(tree, model->times, model->ntimes, false);
    const double treelen_b = treelen + model->time_steps[nodes[tree->root].age];

    // get lineage counts
    lineages.count(tree);
    lineages.nrecombs[root_age]--;

    
    double lnl = 0.0;
            
    // probability of recombination location in tree
    int k = spr.recomb_time;
    lnl += log(lineages.nbranches[k] * model->time_steps[k] /
               (lineages.nrecombs[k] * treelen_b));

    // probability of re-coalescence
    double coal_rates[2*model->ntimes];
    calc_coal_rates_full_tree(model, tree, spr, lineages, coal_rates);
    int j = spr.coal_time;
    int broken_age = nodes[nodes[spr.recomb_node].parent].age;

    // probability of recoalescence on choosen branch
    int ncoals_j = lineages.ncoals[j] 
        - int(j <= broken_age) - int(j == broken_age);
    lnl -= log(ncoals_j);

    // probability of recoalescing in choosen time interval
    if (j < model->ntimes - 2)
        lnl += log(1.0 - exp(- coal_rates[2*j] - 
                             (j>k ? coal_rates[2*j-1] : 0.0)));
    
    // probability of not coalescing before time interval j
    for (int m=2*k; m<2*j-1; m++)
        lnl -= coal_rates[m];
    
    assert(!isinf(lnl));
    return lnl;
}



// calculate the probability of an ARG given the model parameters
double calc_arg_prior(const ArgModel *model, const LocalTrees *trees)
{
    double lnl = 0.0;
    LineageCounts lineages(model->ntimes);
    
    // first tree prior
    lnl += calc_tree_prior(model, trees->front().tree, lineages);

    int end = trees->start_coord;
    for (LocalTrees::const_iterator it=trees->begin(); it != trees->end();) {
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
            const Spr *spr = &it->spr;
            lnl += calc_spr_prob(model, tree, *spr, lineages, treelen);

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
double calc_arg_joint_prob(const ArgModel *model, const Sequences *sequences, 
                           const LocalTrees *trees)
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


double arghmm_likelihood_parsimony(LocalTrees *trees,
                                   double *times, int ntimes,
                                   double mu, 
                                   char **seqs, int nseqs, int seqlen)
{
    /*
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, NULL, 0.0, mu);
    Sequences sequences(seqs, nseqs, seqlen);
    return calc_arg_likelihood_parsimony(&model, &sequences, trees);
    */
    abort();
    return 0.0;
}


double arghmm_prior_prob(LocalTrees *trees,
                         double *times, int ntimes, double *popsizes,
                         double rho)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, rho, 0.0);
    return calc_arg_prior(&model, trees);
}


double arghmm_tree_prior_prob(LocalTrees *trees,
                              double *times, int ntimes, double *popsizes)
{
    // setup model, local trees, sequences
    ArgModel model(ntimes, times, popsizes, 0.0, 0.0);
    LineageCounts lineages(ntimes);

    //printf("%f\n", calc_tree_prior_approx(&model, trees->front().tree,
    //                                      lineages));

    return calc_tree_prior(&model, trees->front().tree, lineages);
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



} // extern "C"


} // namespace arghmm


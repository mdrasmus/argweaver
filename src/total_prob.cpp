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
    
    for (int i=0; i<model->ntimes-1; i++) {
        int a = (i == 0 ? nleaves : lineages.nbranches[i-1]);
        int b = lineages.nbranches[i];

        lnl += log(prob_coal_counts(a, b, model->time_steps[i],
                                    2.0*model->popsizes[i]));
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
double calc_arg_prior(const ArgModel *model, const LocalTrees *trees)
{
    double lnl = 0.0;
    LineageCounts lineages(model->ntimes);

    // TODO: fix this before re-enabling
    // first tree prior
    //lnl += calc_tree_prior(model, trees->front().tree, lineages);


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



//=============================================================================
// old code


/*
// calculate the probability of the sequences given an ARG
double calc_arg_likelihood_parsimony(
    const ArgModel *model, const Sequences *sequences, const LocalTrees *trees)
{
    double lnl = 0.0;
    int nseqs = sequences->get_num_seqs();
    double minlen = model->get_mintime();
    const double log25 = log(.25);
    const double *times = model->times;

    if (trees->nnodes < 3)
        return lnl += log25 * sequences->length();

    // get sequences for trees
    char *seqs[nseqs];
    for (int j=0; j<nseqs; j++)
        seqs[j] = sequences->seqs[trees->seqids[j]];

    int end = trees->start_coord;
    for (LocalTrees::const_iterator it=trees->begin(); it!=trees->end(); ++it) {
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
*/

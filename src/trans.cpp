//=============================================================================
// transitions

#include "common.h"
#include "trans.h"

using namespace spidir;

namespace arghmm {


// This notation matches my math notes
void calc_transition_probs(LocalTree *tree, ArgModel *model,
                           const States &states, LineageCounts *lineages,
                           double **transprob)
{
    // get model parameters
    const int ntimes = model->ntimes;
    const double *times = model->times;
    const double *time_steps = model->time_steps;
    const double *popsizes = model->popsizes;
    const double rho = model->rho;
    const int nstates = states.size();
    
    // get tree informations
    LocalNode *nodes = tree->nodes;
    const int *nbranches = lineages->nbranches;
    const int *nrecombs = lineages->nrecombs;
    const int *ncoals = lineages->ncoals;

    // find root node
    int root = tree->root;
    const int root_age_index = nodes[root].age;
    const double root_age = times[root_age_index];
    const double treelen = get_treelen(tree, times, ntimes, false);
    
    // temporary variables
    double B[ntimes];
    double D[ntimes];
    double E[ntimes];
    double G[ntimes];
    double norecombs[ntimes];
    
    // calculate base case (time=0)
    double treelen_b = treelen + time_steps[root_age_index];
    double C = 0.0;
    B[0] = (nbranches[0] + 1.0) * time_steps[0] / (nrecombs[0] + 1.0);
    D[0] = (1.0 - exp(-max(rho * treelen, rho))) / treelen_b;
    E[0] = (1.0 - exp(-time_steps[0] * nbranches[0]
                      / (2.0 * popsizes[0]))) / ncoals[0];
    G[0] = time_steps[0] / (nrecombs[0] + 1.0);
    norecombs[0] = exp(-max(rho * treelen, rho));

    // calculate all other time points (time>0)
    for (int b=1; b<ntimes-1; b++) {
        // get tree length
        double treelen2 = treelen + times[b];
        double treelen2_b;
        if (b > root_age_index) {
            // add wrapped branch
            treelen2 += times[b] - root_age;

            // add basal branch
            treelen2_b = treelen2 + time_steps[b];
        } else {
            // add basal branch
            treelen2_b = treelen2 + time_steps[root_age_index];
        }

        // due to normalization we do not need exp(-rho * treelen)
        const int l = b - 1;
        C += time_steps[l] * nbranches[l] / (2.0 * popsizes[l]);
        double eC = exp(C);

        B[b] = B[b-1] + (nbranches[b] + 1.0) * time_steps[b] / 
            (nrecombs[b] + 1.0) * eC;
        D[b] = (1.0 - exp(-max(rho * treelen2, rho))) / treelen2_b;
        E[b] = (1.0 - exp(-time_steps[b] * nbranches[b] / 
                          (2.0 * popsizes[b]))) / eC / ncoals[b];
        G[b] = eC * time_steps[b] / (nrecombs[b] + 1.0);
        norecombs[b] = exp(-max(rho * treelen2, rho));
    }
    E[ntimes-2] = 1.0 / exp(C) / ncoals[ntimes-2];

    
    // calculate full state transition matrix
    for (int i=0; i<nstates; i++) {
        const int node1 = states[i].node;
        const int a = states[i].time;
        const int c = nodes[node1].age;
        
        for (int j=0; j<nstates; j++) {
            const int node2 = states[j].node;
            const int b = states[j].time;
            const double I = float(a <= b);
            
            if (node1 != node2)
                transprob[i][j] = D[a] * E[b] * (B[min(a,b)] - I * G[a]);
            else {
                double Bc = (c > 0 ? B[c-1] : 0.0);
                transprob[i][j] = D[a] * E[b] * 
                    (2 * B[min(a,b)] - 2 * I * G[a] - Bc);
                if (a == b)
                    transprob[i][j] += norecombs[a];
            }
        }

        // normalize and convert to log scale
        double sum = 0.0;
        for (int j=0; j<nstates; j++)
            sum += transprob[i][j];
        for (int j=0; j<nstates; j++)
            transprob[i][j] = log(transprob[i][j] / sum);
    }
}



// This notation matches my math notes
void calc_transition_probs_compress(LocalTree *tree, ArgModel *model,
    const States &states, LineageCounts *lineages, TransMatrixCompress *matrix)
{
    // get model parameters
    const int nstates = states.size();
    const int ntimes = model->ntimes;
    const double *times = model->times;
    const double *time_steps = model->time_steps;
    const double *popsizes = model->popsizes;
    const double rho = model->rho;
    
    const int *nbranches = lineages->nbranches;
    const int *nrecombs = lineages->nrecombs;
    const int *ncoals = lineages->ncoals;

    // get matrix fields
    double *B = matrix->B;
    double *D = matrix->D;
    double *E = matrix->E;
    double *G = matrix->G;
    double *norecombs = matrix->norecombs;
    double *sums = matrix->sums;

    // find root node
    int root = tree->root;
    const int root_age_index = tree->nodes[root].age;
    const double root_age = times[root_age_index];
    const double treelen = get_treelen(tree, times, ntimes, false);
    
    
    // base cases (time=0)
    double treelen_b = treelen + time_steps[root_age_index];
    double C = 0.0;
    B[0] = (nbranches[0] + 1.0) * time_steps[0] / (nrecombs[0] + 1.0);
    D[0] = (1.0 - exp(-max(rho * treelen, rho))) / treelen_b;
    E[0] = (1.0 - exp(-time_steps[0] * nbranches[0]
                      / (2.0 * popsizes[0]))) / ncoals[0];
    G[0] = time_steps[0] / (nrecombs[0] + 1.0);
    norecombs[0] = exp(-max(rho * treelen, rho));
    
    // calculate all other time points (time>0)
    for (int b=1; b<ntimes-1; b++) {
        // get tree length
        double treelen2 = treelen + times[b];
        double treelen2_b;
        if (b > root_age_index) {
            // add wrapped branch
            treelen2 += times[b] - root_age;

            // add basal branch
            treelen2_b = treelen2 + time_steps[b];
        } else {
            // add basal branch
            treelen2_b = treelen2 + time_steps[root_age_index];
        }

        // due to normalization we do not need exp(-rho * treelen)
        const int l = b - 1;
        C = C + time_steps[l] * nbranches[l] / (2.0 * popsizes[l]);
        double eC = exp(C);
        B[b] = B[b-1] + (nbranches[b] + 1.0) * time_steps[b] / (nrecombs[b] + 1.0)*eC;
        D[b] = (1.0 - exp(-rho * treelen2)) / treelen2_b;
        E[b] = (1.0 - exp(-time_steps[b] * nbranches[b] / 
                          (2.0 * popsizes[b]))) / eC / ncoals[b];
        G[b] = eC * time_steps[b] / (nrecombs[b] + 1.0);
        norecombs[b] = exp(-rho * treelen2);
    }
    E[ntimes-2] = exp(-C) / ncoals[ntimes-2];

    // get max time
    int maxtime = 0;
    for (int k=0; k<nstates; k++)
        if (maxtime < states[k].time)
            maxtime = states[k].time;

    // get branch ages
    double ages1[tree->nnodes];
    double ages2[tree->nnodes];
    for (int i=0; i<tree->nnodes; i++) {
        ages1[i] = tree->nodes[i].age;
        ages2[i] = (tree->root == i) ? 
            maxtime : tree->nodes[tree->nodes[i].parent].age;
    }

    // get normalization factors for each row of the transition matrix    
    for (int i=0; i<nstates; i++) {
        const int node1 = states[i].node;
        const int a = states[i].time;
        const int c = tree->nodes[node1].age;
        double sum = 0.0;

        for (int b=0; b<ntimes-1; b++) {
            const double I = float(a <= b);
            sum += E[b] * (B[min(a,b)] - I * G[a]) * lineages->ncoals[b];
        }

        for (int b=ages1[node1]; b<=ages2[node1]; b++) {
            const double I = float(a <= b);
            const double Bc = (c > 0 ? B[c-1] : 0.0);
            sum += E[b] * (B[min(a,b)] - I * G[a] - Bc);
        }

        sum *= D[a];
        sum += norecombs[a];
        sums[i] = sum;
    }
}



void calc_transition_probs(LocalTree *tree, ArgModel *model,
                           const States &states, LineageCounts *lineages,
                           TransMatrixCompress *matrix, double **transprob)
{
    // get tree and model information
    const int nstates = states.size();
    const LocalNode *nodes = tree->nodes;

    // get matrix fields
    const double *B = matrix->B;
    const double *D = matrix->D;
    const double *E = matrix->E;
    const double *G = matrix->G;
    
    // calculate full state transition matrix
    for (int i=0; i<nstates; i++) {
        const int node1 = states[i].node;
        const int a = states[i].time;
        const int c = nodes[node1].age;
        
        for (int j=0; j<nstates; j++) {
            const int node2 = states[j].node;
            const int b = states[j].time;
            const double I = float(a <= b);
            
            if (node1 != node2)
                transprob[i][j] = D[a] * E[b] * (B[min(a,b)] - I * G[a]);
            else {
                double Bc = (c > 0 ? B[c-1] : 0.0);
                transprob[i][j] = D[a] * E[b] * 
                    (2 * B[min(a,b)] - 2 * I * G[a] - Bc);
                if (a == b)
                    transprob[i][j] += matrix->norecombs[a];
            }
        }

        // normalize and convert to log scale
        double sum = 0.0;
        for (int j=0; j<nstates; j++)
            sum += transprob[i][j];
        for (int j=0; j<nstates; j++)
            transprob[i][j] = log(transprob[i][j] / sum);
    }
}




void get_deterministic_transitions(
    LocalTree *tree, LocalTree *last_tree, const Spr &spr, int *mapping,
    const States &states1, const States &states2,
    int ntimes, int *next_states)
{
    // recomb_node in tree and last_tree
    // coal_node in last_tree
    
    const int nnodes = tree->nnodes;
    const LocalNode *nodes = tree->nodes;
    const LocalNode *last_nodes = last_tree->nodes;
    const int nstates1 = states1.size();

    // make state lookup
    NodeStateLookup state2_lookup(states2, nnodes);
    

    for (int i=0; i<nstates1; i++) {
        const int node1 = states1[i].node;
        const int time1 = states1[i].time;

        if (node1 == spr.coal_node && time1 == spr.coal_time) {
            // not a deterministic case
            next_states[i] = -1;
        
        } else if (node1 != spr.recomb_node) {
            // SPR only removes a subset of descendents, if any
            // trace up from remaining leaf to find correct new state

            int node2;
            const LocalNode *node = &last_nodes[node1];
            bool disrupt = false;
            
            if (node->child[0] == -1) {
                // SPR can't disrupt leaf branch
                node2 = node1;

            } else {
                int child1 = node->child[0];
                int child2 = node->child[1];
                
                if (spr.recomb_node == child1) {
                    // right child is not disrupted
                    node2 = mapping[child2];
                    disrupt = true;
                } else if (spr.recomb_node == child2) {
                    // left child is not disrupted
                    node2 = mapping[child1];
                    disrupt = true;
                } else {
                    // node is not disrupted
                    node2 = mapping[node1];
                }
            }

            // optionally walk up
            if ((spr.coal_node == node1 && spr.coal_time < time1) || 
                (mapping[spr.coal_node] == node2 && spr.coal_time < time1) ||
                (disrupt && mapping[spr.coal_node] == node2 && 
                 spr.coal_time <= time1))
            {
                // XXX: need to walk up for one more case
                // need to change < time1 to <== time for disrupt case
                // coal occurs under us
                node2 = nodes[node2].parent;
            }
            
            assert(nodes[node2].age <= time1);
            int p = nodes[node2].parent;
            if (p != -1)
                assert(time1 <= nodes[p].age);
            
            // set next state
            next_states[i] = state2_lookup.lookup(node2, time1);

        } else {
            // SPR is on same branch as new chromosome
            if (spr.recomb_time >= time1) {
                // we move with SPR subtree
                // TODO: we could probabilistically have subtree move
                // out from underneath.
                next_states[i] = state2_lookup.lookup(
                    mapping[spr.recomb_node], time1);

            } else {
                // SPR should not be able to coal back onto same branch
                // this would be a self cycle
                assert(spr.coal_node != node1);
                
                // SPR subtree moves out from underneath us
                // therefore therefore the new chromosome coalesces with
                // the branch above the subtree

                // search up for parent
                int parent = last_nodes[spr.recomb_node].parent;
                int time2 = last_nodes[parent].age;
                int node2;

                // find other child
                const int *c = last_nodes[parent].child;
                int other = (c[1] == spr.recomb_node ? c[0] : c[1]);

                // find new state in tree
                node2 = (other == spr.coal_node ? 
                         nodes[mapping[other]].parent : mapping[other]);
                next_states[i] = state2_lookup.lookup(node2, time2);
            }
        }
    }
    
}



void calc_transition_probs_switch(
    LocalTree *tree, LocalTree *last_tree, const Spr &spr, int *mapping,
    const States &states1, const States &states2,
    ArgModel *model, LineageCounts *lineages, double **transprob)
{
    // get tree information
    const LocalNode *nodes = tree->nodes;
    const LocalNode *last_nodes = last_tree->nodes;
    const int *nbranches = lineages->nbranches;
    const int *ncoals = lineages->ncoals;
    const int nstates1 = states1.size();
    const int nstates2 = states2.size();
    
    // get model parameters
    const int ntimes = model->ntimes;
    const double *time_steps = model->time_steps;
    const double *popsizes = model->popsizes;
    

    // get deterministic transitions
    int determ[nstates1];
    get_deterministic_transitions(tree, last_tree, spr, mapping,
                                  states1, states2, ntimes, determ);
    

    for (int i=0; i<nstates1; i++) {
        const int node1 = states1[i].node;
        const int time1 = states1[i].time;
        
        if (node1 != spr.coal_node || time1 != spr.coal_time) {
            // deterministic transition case
            assert (determ[i] != -1);
            for (int j=0; j<nstates2; j++)
                transprob[i][j] = -INFINITY;
            transprob[i][determ[i]] = 0.0;
            
        } else {
            // probabilistic transition case
            
            // determine if node1 is still here or not
            int node3;
            int last_parent = last_nodes[spr.recomb_node].parent;
            if (last_parent == node1) {
                // recomb breaks node1 branch, we need to use the other child
                const int *c = last_nodes[last_parent].child;
                node3 = mapping[c[1] == spr.recomb_node ? c[0] : c[1]];
            } else {
                node3 = mapping[node1];
            }

            // find parent of recomb_branch and node1
            int last_parent_age = last_nodes[last_parent].age;

            
            int parent = nodes[mapping[spr.recomb_node]].parent;
            assert(parent == nodes[node3].parent);

            for (int j=0; j<nstates2; j++) {
                const int node2 = states2[j].node;
                const int time2 = states2[j].time;
                
                transprob[i][j] = 0.0;
                if (!((node2 == mapping[spr.recomb_node] 
                       && time2 >= spr.recomb_time) ||
                      (node2 == node3 && time2 == time1) ||
                      (node2 == parent && time2 == time1)))
                    // not a probabilistic transition
                    continue;

                // get lineage counts
                // remove recombination branch and add new branch
                int kbn = nbranches[time2];
                int kcn = ncoals[time2] + 1;
                if (time2 < nodes[parent].age) {
                    kbn -= 1;
                    kcn -= 1;
                }
                if (time2 < time1)
                    kbn += 1;

                double sum = 0.0;
                for (int m=spr.recomb_time; m<time2; m++) {
                    const int nlineages = nbranches[m] + 1
                        - (m < last_parent_age ? 1 : 0);
                    sum += time_steps[m] * nlineages / (2.0 * popsizes[m]);
                }
                transprob[i][j] =
                    (1.0 - exp(- time_steps[time2] * kbn /  
                               (2.0 * popsizes[time2]))) / kcn * exp(-sum);
            }

            // normalize row to ensure they add up to one
            double sum = 0.0;
            for (int j=0; j<nstates2; j++)
                sum += transprob[i][j];
            for (int j=0; j<nstates2; j++) {
                double x = transprob[i][j];
                if (sum > 0.0 and x > 0.0)
                    transprob[i][j] = log(x / sum);
                else
                    transprob[i][j] = -INFINITY;
            }
        }
    }
}



void calc_transition_probs_switch_compress(
    LocalTree *tree, LocalTree *last_tree, const Spr &spr, int *mapping,
    const States &states1, const States &states2,
    ArgModel *model, LineageCounts *lineages, 
    TransMatrixSwitchCompress *transmat_switch_compress)
{
    // get tree information
    const LocalNode *nodes = tree->nodes;
    const LocalNode *last_nodes = last_tree->nodes;
    const int *nbranches = lineages->nbranches;
    const int *ncoals = lineages->ncoals;
    const int nstates1 = states1.size();
    const int nstates2 = states2.size();
    
    // get model parameters
    const int ntimes = model->ntimes;
    const double *time_steps = model->time_steps;
    const double *popsizes = model->popsizes;
    
    
    // get deterministic transitions
    get_deterministic_transitions(tree, last_tree, spr, mapping,
         states1, states2, ntimes, transmat_switch_compress->determ);
    
    // find probabilitistic transition source state
    int src = -1;
    for (int i=0; i<nstates1; i++) {
        if (states1[i].node == spr.coal_node && 
            states1[i].time == spr.coal_time) {
            src = i;
            break;
        }
    }
    assert(src != -1);

    transmat_switch_compress->probsrc = src;    
    int node1 = states1[src].node;
    int time1 = states1[src].time;
    
    // determine if node1 is still here or not
    int node3;
    int last_parent = last_nodes[spr.recomb_node].parent;
    if (last_parent == node1) {
        // recomb breaks node1 branch, we need to use the other child
        const int *c = last_nodes[last_parent].child;
        node3 = mapping[c[1] == spr.recomb_node ? c[0] : c[1]];
    } else {
        node3 = mapping[node1];
    }

    // find parent of recomb_branch and node1
    int last_parent_age = last_nodes[last_parent].age;

            
    int parent = nodes[mapping[spr.recomb_node]].parent;
    assert(parent == nodes[node3].parent);
    
    double *probrow = transmat_switch_compress->probrow;
    for (int j=0; j<nstates2; j++) {
        const int node2 = states2[j].node;
        const int time2 = states2[j].time;
        
        probrow[j] = 0.0;
        if (!((node2 == mapping[spr.recomb_node] 
               && time2 >= spr.recomb_time) ||
              (node2 == node3 && time2 == time1) ||
              (node2 == parent && time2 == time1)))
            // not a probabilistic transition
            continue;
        
        // get lineage counts
        // remove recombination branch and add new branch
        int kbn = nbranches[time2];
        int kcn = ncoals[time2] + 1;
        if (time2 < nodes[parent].age) {
            kbn -= 1;
            kcn -= 1;
        }
        if (time2 < time1)
            kbn += 1;

        double sum = 0.0;
        for (int m=spr.recomb_time; m<time2; m++) {
            const int nlineages = nbranches[m] + 1
                - (m < last_parent_age ? 1 : 0);
            sum += time_steps[m] * nlineages / (2.0 * popsizes[m]);
        }
        probrow[j] = (1.0 - exp(- time_steps[time2] * kbn /  
                                (2.0 * popsizes[time2]))) / kcn * exp(-sum);
    }

    // normalize row to ensure they add up to one
    double sum = 0.0;
    for (int j=0; j<nstates2; j++)
        sum += probrow[j];
    for (int j=0; j<nstates2; j++) {
        double x = probrow[j];
        if (sum > 0.0 and x > 0.0)
            probrow[j] = log(x / sum);
        else
            probrow[j] = -INFINITY;
    }
}


void calc_transition_switch_probs(TransMatrixSwitchCompress *matrix, 
                                  double **transprob)
{
    for (int i=0; i<matrix->nstates1; i++) {
        for (int j=0; j<matrix->nstates2; j++) {
            transprob[i][j] = matrix->get_transition_prob(i, j);
        }
    }
}


void calc_state_priors(const States &states, LineageCounts *lineages, 
                       ArgModel *model, double *priors)
{
    const int nstates = states.size();
    const double *time_steps = model->time_steps;
    const double *popsizes = model->popsizes;
    const int *nbranches = lineages->nbranches;
    const int *ncoals = lineages->ncoals;
    
    for (int i=0; i<nstates; i++) {
        int b = states[i].time;

        double sum = 0.0;
        for (int m=0; m<b; m++)
            sum += time_steps[m] * nbranches[m] / (2.0 * popsizes[m]);
        
        priors[i] = log((1.0 - exp(- time_steps[b] * nbranches[b] /
                          (2.0 * popsizes[b]))) / ncoals[b] * exp(-sum)); 
    }
}


//=============================================================================
// C interface
extern "C" {

double **new_transition_probs(int nnodes, int *ptree, 
                              int *ages, double treelen,
                              intstate *istates, int nstates,
                              int ntimes, double *times, double *time_steps,
                              int *nbranches, int *nrecombs, int *ncoals, 
                              double *popsizes, double rho)
{

    // setup model, local tree, states
    ArgModel model(ntimes, times, popsizes, rho, 0.0);
    LocalTree tree(ptree, nnodes, ages);
    LineageCounts lineages(ntimes);
    lineages.count(&tree);
    States states;
    make_states(istates, nstates, states);

    double **transprob = new_matrix<double>(nstates, nstates);
    calc_transition_probs(&tree, &model, states, &lineages, transprob);
    return transprob;
}


double **new_transition_probs_switch(
    int *ptree, int *last_ptree, int nnodes, 
    int recomb_node, int recomb_time, int coal_node, int coal_time,
    int *ages_index, int *last_ages_index,
    double treelen, double last_treelen,
    intstate *istates1, int nstates1,
    intstate *istates2, int nstates2,
    
    int ntimes, double *times, double *time_steps,
    int *nbranches, int *nrecombs, int *ncoals, 
    double *popsizes, double rho)
{
    // setup model
    ArgModel model(ntimes, times, popsizes, rho, 0.0);
    
    // setup local trees
    LocalTree tree(ptree, nnodes, ages_index);
    LocalTree last_tree(last_ptree, nnodes, last_ages_index);
    Spr spr(recomb_node, recomb_time, coal_node, coal_time);
    int mapping[nnodes];
    make_node_mapping(mapping, nnodes, last_ptree, recomb_node);
    LineageCounts lineages(ntimes);
    lineages.count(&last_tree);
    
    // setup states
    States states1, states2;
    make_states(istates1, nstates1, states1);
    make_states(istates2, nstates2, states2);
    
    double **transprob = new_matrix<double>(nstates1, nstates2);
    calc_transition_probs_switch(&tree, &last_tree, spr, mapping,
        states1, states2, &model, &lineages, transprob);
    return transprob;
}


void delete_transition_probs(double **transmat, int nstates)
{
    delete_matrix<double>(transmat, nstates);
}

} // extern C



} // namespace arghmm

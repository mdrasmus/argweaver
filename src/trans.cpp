//=============================================================================
// transitions

#include "common.h"
#include "matrices.h"
#include "total_prob.h"
#include "thread.h"
#include "trans.h"


namespace arghmm {



// Calculate transition probability within a local block
void calc_transition_probs(const LocalTree *tree, const ArgModel *model,
    const States &states, const LineageCounts *lineages, TransMatrix *matrix)
{
    // get model parameters
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
    G[0] = time_steps[0] * ((nbranches[0] + 1.0) / (nrecombs[0] + 1.0) -
                            (nbranches[0] / (nrecombs[0] + 1.0 + 
                                             int(0 < root_age_index))));
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

        const int l = b - 1;
        C = C + time_steps[l] * nbranches[l] / (2.0 * popsizes[l]);
        const double eC = exp(C);
        B[b] = B[b-1] + (nbranches[b] + 1.0) * time_steps[b] / (nrecombs[b] + 1.0)*eC;
        D[b] = (1.0 - exp(-rho * treelen2)) / treelen2_b;
        E[b] = (1.0 - exp(-time_steps[b] * nbranches[b] / 
                          (2.0 * popsizes[b]))) / eC / ncoals[b];
        G[b] = eC * time_steps[b] * 
            ((nbranches[b] + 1.0) / (nrecombs[b] + 1.0) -
             nbranches[b] / (nrecombs[b] + 1.0 + int(b < root_age_index)));
        norecombs[b] = exp(-max(rho * treelen2, rho));
    }
    E[ntimes-2] = exp(-C) / ncoals[ntimes-2];
}


// Calculate transition probability within a local block
void calc_transition_probs_internal(const LocalTree *tree, 
    const ArgModel *model, const States &states, const LineageCounts *lineages,
    TransMatrix *matrix)
{
    matrix->internal = true;

    // get model parameters
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
    matrix->internal = true;

    // find root node
    const int *c = tree->nodes[tree->root].child;
    const int subtree_root = c[0];
    const int maintree_root = c[1];
    const int root_age_index = tree->nodes[maintree_root].age;
    const double root_age = times[root_age_index];
    const double subtree_age = times[tree->nodes[subtree_root].age];
    const double treelen = get_treelen_internal(tree, times, ntimes);
    
    // base cases (time=0)
    double treelen_b = treelen + time_steps[root_age_index];
    double C = 0.0;
    B[0] = (nbranches[0] + 1.0) * time_steps[0] / (nrecombs[0] + 1.0);
    D[0] = (1.0 - exp(-max(rho * treelen, rho))) / treelen_b;
    E[0] = (1.0 - exp(-time_steps[0] * nbranches[0]
                      / (2.0 * popsizes[0]))) / ncoals[0];
    G[0] = time_steps[0] * ((nbranches[0] + 1.0) / (nrecombs[0] + 1.0) -
                            (nbranches[0] / (nrecombs[0] + 1.0 + 
                                             int(0 < root_age_index))));
    norecombs[0] = exp(-max(rho * treelen, rho));
    
    // calculate all other time points (time>0)
    for (int b=1; b<ntimes-1; b++) {
        // get tree length
        double treelen2 = treelen + times[b] - subtree_age;
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
        
        const int l = b - 1;
        C = C + time_steps[l] * nbranches[l] / (2.0 * popsizes[l]);
        const double eC = exp(C);
        B[b] = B[b-1] + (nbranches[b] + 1.0) * time_steps[b] / (nrecombs[b] + 1.0)*eC;
        D[b] = (1.0 - exp(-rho * treelen2)) / treelen2_b;
        E[b] = (1.0 - exp(-time_steps[b] * nbranches[b] / 
                          (2.0 * popsizes[b]))) / eC / ncoals[b];
        G[b] = eC * time_steps[b] * 
            ((nbranches[b] + 1.0) / (nrecombs[b] + 1.0) -
             nbranches[b] / (nrecombs[b] + 1.0 + int(b < root_age_index)));
        norecombs[b] = exp(-max(rho * treelen2, rho));
    }
    E[ntimes-2] = exp(-C) / ncoals[ntimes-2];
}


// Copies transition probability to dense probability matrix
void get_transition_probs(const LocalTree *tree, const ArgModel *model,
                          const States &states, const LineageCounts *lineages,
                          const TransMatrix *matrix, double **transprob)
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
                    (2*B[min(a,b)] - 2*I*G[a] - Bc);
                if (a == b)
                    transprob[i][j] += matrix->norecombs[a];
            }
            
            transprob[i][j] = log(transprob[i][j]);
        }
    }
}



// Calculate transition probability within a local block
// Stores transition probabilities in dense matrix
void calc_transition_probs(const LocalTree *tree, const ArgModel *model,
                           const States &states, const LineageCounts *lineages,
                           double **transprob)
{
    TransMatrix matrix(model->ntimes, states.size());
    calc_transition_probs(tree, model, states, lineages, &matrix);
    get_transition_probs(tree, model, states, lineages, &matrix, transprob);
}


//=============================================================================
// functions for switch matrix calculation


void get_deterministic_transitions(
    const LocalTree *tree, const LocalTree *last_tree, 
    const Spr &spr, const int *mapping,
    const States &states1, const States &states2,
    int ntimes, int *next_states, bool internal)
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

        if ((node1 == spr.coal_node && time1 == spr.coal_time) ||
            (node1 == spr.recomb_node && time1 == spr.recomb_time)) {
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
                const int child1 = node->child[0];
                const int child2 = node->child[1];
                
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
                // coal occurs under us
                node2 = nodes[node2].parent;
            }
            
            if (internal && nodes[node2].age > time1) {
                // this is a deadend
                next_states[i] = -1;
                continue;
            }

            assert(nodes[node2].age <= time1);
            const int p = nodes[node2].parent;
            if (p != -1) {
                if (internal && time1 > nodes[p].age) {
                    // this is a deadend
                    next_states[i] = -1;
                    continue;
                }
                
                assert(time1 <= nodes[p].age);
            }
            
            // set next state
            next_states[i] = state2_lookup.lookup(node2, time1);

        } else {
            // SPR is on same branch as new chromosome
            if (spr.recomb_time > time1) {
                // we move with SPR subtree
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
                const int parent = last_nodes[spr.recomb_node].parent;
                const int time2 = last_nodes[parent].age;

                // find other child
                const int *c = last_nodes[parent].child;
                const int other = (c[1] == spr.recomb_node ? c[0] : c[1]);

                // find new state in tree
                const int node2 = (other == spr.coal_node ? 
                    nodes[mapping[other]].parent : mapping[other]);
                next_states[i] = state2_lookup.lookup(node2, time2);
            }
        }
    }
}


void get_recomb_transition_switch(
    const LocalTree *tree, const LocalTree *last_tree, 
    const Spr &spr, const int *mapping,
    const States &states1, const States &states2, int next_states[2])
{
    const int nnodes = tree->nnodes;
    const LocalNode *nodes = tree->nodes;
    const LocalNode *last_nodes = last_tree->nodes;

    // make state lookup
    NodeStateLookup state2_lookup(states2, nnodes);
    
    // SPR subtree moves out from underneath us
    // therefore therefore the new chromosome coalesces with
    // the branch above the subtree

    // search up for parent
    int parent = last_nodes[spr.recomb_node].parent;
    int time2 = last_nodes[parent].age;
    int node2;

    // find other child
    int other = last_tree->get_sibling(spr.recomb_node);

    // find new state in tree
    node2 = (other == spr.coal_node ? 
             nodes[mapping[other]].parent : mapping[other]);

    // stay case
    next_states[0] = state2_lookup.lookup(mapping[spr.recomb_node], 
                                          spr.recomb_time);
    // escape case
    next_states[1] = state2_lookup.lookup(node2, time2);    
}


double calc_recomb_recoal(
    const LocalTree *last_tree, const ArgModel *model, 
    const LineageCounts *lineages, 
    const Spr &spr, const State state1, 
    const int recomb_parent_age, double last_treelen,
    const bool internal=false)
{
    const int *nbranches = lineages->nbranches;
    const int *ncoals = lineages->ncoals;

    // get times
    int a = state1.time;
    const int k = spr.recomb_time;
    const int j = spr.coal_time;

    double last_treelen_b;


    int root_age;
    if (internal) {
        int subtree_root = last_tree->nodes[last_tree->root].child[0];
        int maintree_root = last_tree->nodes[last_tree->root].child[1];
        root_age = last_tree->nodes[maintree_root].age;

        assert(spr.recomb_node != subtree_root);

        // detect sprs onto subtree root branch
        if (spr.coal_node == subtree_root) {
            if (a < spr.coal_time)
                return 0.0;
            if (spr.recomb_node == maintree_root) {
                if (state1.node != maintree_root) {
                    return 0.0;
                }
            }
        }

        // detect sprs from anywhere in maintree to anywhere in subtree
        // state1 must be topologically above recomb point
        int ptr = spr.coal_node;
        int ptr2 = -1, ptr3 = -1;
        while (ptr != last_tree->root) {
            ptr2 = ptr;
            ptr = last_tree->nodes[ptr].parent;
        }        
        ptr = spr.recomb_node;
        while (ptr != last_tree->root) {
            ptr3 = ptr;
            ptr = last_tree->nodes[ptr].parent;
        }
        if (ptr2 == subtree_root && ptr3 == maintree_root &&
            state1.time == spr.recomb_time) {
            // check is needed, ensure state1 is not topologically below recomb
            ptr = last_tree->nodes[state1.node].parent;
            while (ptr != last_tree->root) {
                if (ptr == spr.recomb_node)
                    return 0.0;
                ptr = last_tree->nodes[ptr].parent;
            }
        }
        
        
        last_treelen += model->times[a] - 
            model->times[last_tree->nodes[subtree_root].age];

        if (a > root_age) {
            // add wrapped branch
            last_treelen += model->times[a] - model->times[root_age];

            // add basal branch
            last_treelen_b = last_treelen + model->time_steps[a];
        } else {
            // add basal branch
            last_treelen_b = last_treelen + model->time_steps[root_age];
        }
    } else {
        root_age = last_tree->nodes[last_tree->root].age;
        last_treelen = get_treelen_branch(
            last_tree, model->times, model->ntimes,
            state1.node, state1.time, last_treelen, false);
        last_treelen_b = last_treelen + get_basal_branch(
            last_tree, model->times, model->ntimes,
            state1.node, state1.time);
    }


    // recomb prob
    //int root_age = last_tree->nodes[last_tree->root].age;
    int nbranches_k = lineages->nbranches[k] + int(k < a);
    int nrecombs_k = lineages->nrecombs[k] + int(k <= a) + 
        int(k == a) - int(k >= max(root_age, a));
    double recomb_prob = nbranches_k * model->time_steps[k] /
        (nrecombs_k * last_treelen_b) * 
        (1.0 - exp(-max(model->rho * last_treelen, model->rho)));
  

    // coal prob
    double sum = 0.0;
    for (int m=k; m<j; m++) {
        int nbranches_m = lineages->nbranches[m] - int(m < recomb_parent_age) 
            + int(m < a);
        sum += model->time_steps[m] * nbranches_m / (2.0 * model->popsizes[m]);
    }
    int nbranches_j = nbranches[j] - int(j < recomb_parent_age) + int(j < a);
    int ncoals_j = ncoals[j] - int(j <= recomb_parent_age)
        - int(j == recomb_parent_age) + int(j <= a) + int(j == a);

    if (internal) {
        int subtree_root = last_tree->nodes[last_tree->root].child[0];
        int maintree_root = last_tree->nodes[last_tree->root].child[1];
        if (spr.recomb_node == maintree_root) {
            // special cases for properly calculating ncoals_j and nbranches_j
            if (spr.coal_time >= last_tree->nodes[subtree_root].age) {
                nbranches_j = 1;
                ncoals_j++;
            }
        }
    }

    
    double p = recomb_prob * exp(-sum) / ncoals_j;

    if (j < model->ntimes - 2)
        p *= 1.0 - exp(-model->time_steps[j] * nbranches_j / 
                       (2.0*model->popsizes[j]));

    if (ncoals_j <= 0 || nbranches_j <= 0 || 
        nrecombs_k <= 0 || nbranches_k <= 0) {
        printf("counts %d %d %d %d %e\n", 
               ncoals_j, nbranches_j, nrecombs_k, nbranches_k, p);
        assert(false);
    }

    assert(!isnan(p) && p>0);
    return p;
}



void calc_transition_probs_switch(
    const LocalTree *tree, const LocalTree *last_tree, 
    const Spr &spr, const int *mapping,
    const States &states1, const States &states2,
    const ArgModel *model, const LineageCounts *lineages, 
    TransMatrixSwitch *transmat_switch)
{
    const int nstates1 = states1.size();
    const int nstates2 = states2.size();
    int recomb_parent_age;        
    
    
    assert(assert_spr(last_tree, tree, &spr, mapping));

    double last_treelen = get_treelen(
            last_tree, model->times, model->ntimes, false);
    

    // get deterministic transitions
    get_deterministic_transitions(tree, last_tree, spr, mapping,
         states1, states2, model->ntimes, transmat_switch->determ);
    for (int i=0; i<nstates1; i++) {
        int j = transmat_switch->determ[i];
        if (j >= 0) {
            if (states1[i].node == spr.recomb_node && 
                states1[i].time > spr.recomb_time)
                recomb_parent_age = states1[i].time;
            else
                recomb_parent_age = last_tree->nodes[
                    last_tree->nodes[spr.recomb_node].parent].age;

            //printf("det %d %d\n", i, j);
            transmat_switch->determprob[i] = log(calc_recomb_recoal(
              last_tree, model, lineages, spr, 
              states1[i], recomb_parent_age, last_treelen));
        }
    }
    
    
    // find probabilitistic transition source states
    int recoalsrc = -1;
    int recombsrc = -1;
    for (int i=0; i<nstates1; i++) {
        if (states1[i].node == spr.recomb_node && 
            states1[i].time == spr.recomb_time) {
            recombsrc = i;
        } else if (states1[i].node == spr.coal_node && 
            states1[i].time == spr.coal_time) {
            recoalsrc = i;
        }
    }
    assert(recoalsrc != -1);
    assert(recombsrc != -1);
    transmat_switch->recoalsrc = recoalsrc;
    transmat_switch->recombsrc = recombsrc;


    
    // compute recomb case
    // [0] = stay, [1] = escape
    int recomb_next_states[2];
    get_recomb_transition_switch(tree, last_tree, spr, mapping,
                                 states1, states2, recomb_next_states);
    for (int j=0; j<nstates2; j++)
        transmat_switch->recombrow[j] = -INFINITY;
    
    {
        // stay case (recomb above)
        int j = recomb_next_states[0];
        recomb_parent_age = last_tree->nodes[
            last_tree->nodes[spr.recomb_node].parent].age;
        transmat_switch->recombrow[j] = log(calc_recomb_recoal(
            last_tree, model, lineages, spr, states1[recombsrc], 
            recomb_parent_age, last_treelen));
        
        // escape case (recomb below)
        j = recomb_next_states[1];
        recomb_parent_age = states1[recombsrc].time;
        transmat_switch->recombrow[j] = log(calc_recomb_recoal(
            last_tree, model, lineages, spr, states1[recombsrc], 
            recomb_parent_age, last_treelen));
    }
    

    // compute recoal case
    int node1 = states1[recoalsrc].node;
    int time1 = states1[recoalsrc].time;
    
    // determine if node1 is still here or not
    int node3;
    int last_parent = last_tree->nodes[spr.recomb_node].parent;
    if (last_parent == node1) {
        // recomb breaks node1 branch, we need to use the other child
        const int *c = last_tree->nodes[last_parent].child;
        node3 = mapping[c[1] == spr.recomb_node ? c[0] : c[1]];
    } else {
        node3 = mapping[node1];
    }
            
    int parent = tree->nodes[mapping[spr.recomb_node]].parent;
    assert(parent == tree->nodes[node3].parent);
    
    for (int j=0; j<nstates2; j++) {
        const int node2 = states2[j].node;
        const int time2 = states2[j].time;
                
        transmat_switch->recoalrow[j] = -INFINITY;
        if (!((node2 == mapping[spr.recomb_node] 
               && time2 >= spr.recomb_time) ||
              (node2 == node3 && time2 == time1) ||
              (node2 == parent && time2 == time1)))
            // not a probabilistic transition
            continue;

        recomb_parent_age = last_tree->nodes[last_tree->nodes[spr.recomb_node].parent].age;
        Spr spr2 = spr;
        spr2.coal_time = time2;
        transmat_switch->recoalrow[j] = log(calc_recomb_recoal(
            last_tree, model, lineages, spr2, states1[recoalsrc],
            recomb_parent_age, last_treelen));
    }
}


void calc_transition_probs_switch_internal(
    const LocalTree *tree, const LocalTree *last_tree, 
    const Spr &spr, const int *mapping,
    const States &states1, const States &states2,
    const ArgModel *model, const LineageCounts *lineages, 
    TransMatrixSwitch *transmat_switch)
{
    const int nstates1 = states1.size();
    const int nstates2 = states2.size();
    int recomb_parent_age;
    
    // remove from the top case
    if (nstates1 == 0) {
        // switching between two completely specified blocks
        if (nstates2 == 0) {
            transmat_switch->determ[0] = 0;
            transmat_switch->determprob[0] = 0.0;
            transmat_switch->recoalsrc = -1;
            transmat_switch->recombsrc = -1;
            return;
        }

        assert(spr.coal_node == last_tree->root);
        const int maintree_root = tree->nodes[tree->root].child[1];
        for (int j=0; j<nstates2; j++) {
            if (states2[j].node == maintree_root &&
                states2[j].time == spr.coal_time) {
                transmat_switch->determ[0] = j;
                transmat_switch->determprob[0] = 0.0;
                transmat_switch->recoalsrc = -1;
                transmat_switch->recombsrc = -1;
                return;
            }
        }

        assert(false);
    }

    double last_treelen = get_treelen_internal(
            last_tree, model->times, model->ntimes);


    // fall off the top case
    if (nstates2 == 0) {
        fill(transmat_switch->determ, transmat_switch->determ + nstates1, 0.0);

        for (int i=0; i<nstates1; i++) {
            if (states1[i].node == spr.recomb_node && 
                states1[i].time > spr.recomb_time)
                recomb_parent_age = states1[i].time;
            else
                recomb_parent_age = last_tree->nodes[
                    last_tree->nodes[spr.recomb_node].parent].age;
                
            transmat_switch->determprob[i] = log(calc_recomb_recoal(
                last_tree, model, lineages, spr, 
                states1[i], recomb_parent_age, last_treelen, true));
        }

        transmat_switch->recoalsrc = -1;
        transmat_switch->recombsrc = -1;
        return;
    }

    
    // get deterministic transitions
    get_deterministic_transitions(tree, last_tree, spr, mapping,
        states1, states2, model->ntimes, transmat_switch->determ,
        true);
    for (int i=0; i<nstates1; i++) {
        int j = transmat_switch->determ[i];
        if (j >= 0) {
            if (states1[i].node == spr.recomb_node && 
                states1[i].time > spr.recomb_time)
                recomb_parent_age = states1[i].time;
            else
                recomb_parent_age = last_tree->nodes[
                    last_tree->nodes[spr.recomb_node].parent].age;

            transmat_switch->determprob[i] = log(calc_recomb_recoal(
              last_tree, model, lineages, spr, 
              states1[i], recomb_parent_age, last_treelen, true));
        }
    }

    
    // find probabilitistic transition source states
    int recoalsrc = -1;
    int recombsrc = -1;
    for (int i=0; i<nstates1; i++) {
        if (states1[i].node == spr.recomb_node && 
            states1[i].time == spr.recomb_time) {
            recombsrc = i;
        } else if (states1[i].node == spr.coal_node && 
            states1[i].time == spr.coal_time) {
            recoalsrc = i;
        }
    }
    transmat_switch->recoalsrc = recoalsrc;
    transmat_switch->recombsrc = recombsrc;


    
    // compute recomb case
    // [0] = stay, [1] = escape
    if (recombsrc != -1) {
        int recomb_next_states[2];
        get_recomb_transition_switch(tree, last_tree, spr, mapping,
                                     states1, states2, recomb_next_states);
        for (int j=0; j<nstates2; j++)
            transmat_switch->recombrow[j] = -INFINITY;
    
        // stay case (recomb above)
        int j = recomb_next_states[0];
        if (j != -1) {
            recomb_parent_age = last_tree->nodes[
                last_tree->nodes[spr.recomb_node].parent].age;
            transmat_switch->recombrow[j] = log(calc_recomb_recoal(
                last_tree, model, lineages, spr, states1[recombsrc],
                recomb_parent_age, last_treelen, true));
            assert(!isnan(transmat_switch->recombrow[j]));
        }

        // escape case (recomb below)
        j = recomb_next_states[1];
        if (j != -1) {
            recomb_parent_age = states1[recombsrc].time;
            transmat_switch->recombrow[j] = log(calc_recomb_recoal(
                last_tree, model, lineages, spr, states1[recombsrc],
                recomb_parent_age, last_treelen, true));
            assert(!isnan(transmat_switch->recombrow[j]));
        }
    }
    

    // compute recoal case
    if (recoalsrc == -1) {
        for (int j=0; j<nstates2; j++)
            transmat_switch->recoalrow[j] = -INFINITY;
        return;
    }

    int node1 = states1[recoalsrc].node;
    int time1 = states1[recoalsrc].time;
    
    // determine if node1 is still here or not
    int node3;
    int last_parent = last_tree->nodes[spr.recomb_node].parent;
    if (last_parent == node1) {
        // recomb breaks node1 branch, we need to use the other child
        const int *c = last_tree->nodes[last_parent].child;
        node3 = mapping[c[1] == spr.recomb_node ? c[0] : c[1]];
    } else {
        node3 = mapping[node1];
    }
            
    int parent = tree->nodes[mapping[spr.recomb_node]].parent;
    assert(parent == tree->nodes[node3].parent);
    
    for (int j=0; j<nstates2; j++) {
        const int node2 = states2[j].node;
        const int time2 = states2[j].time;
                
        transmat_switch->recoalrow[j] = -INFINITY;
        if (!((node2 == mapping[spr.recomb_node] 
               && time2 >= spr.recomb_time) ||
              (node2 == node3 && time2 == time1) ||
              (node2 == parent && time2 == time1)))
            // not a probabilistic transition
            continue;

        recomb_parent_age = last_tree->nodes[last_tree->nodes[spr.recomb_node].parent].age;
        Spr spr2 = spr;
        spr2.coal_time = time2;
        transmat_switch->recoalrow[j] = log(calc_recomb_recoal(
            last_tree, model, lineages, spr2, states1[recoalsrc],
            recomb_parent_age, last_treelen, true));
    }
}



void get_transition_probs_switch(const TransMatrixSwitch *matrix, 
                                 double **transprob)
{
    for (int i=0; i<matrix->nstates1; i++) {
        for (int j=0; j<matrix->nstates2; j++) {
            transprob[i][j] = matrix->get_log(i, j);
        }
    }
}


void calc_transition_probs_switch(
    const LocalTree *tree, const LocalTree *last_tree, 
    const Spr &spr, const int *mapping,
    const States &states1, const States &states2,
    const ArgModel *model, const LineageCounts *lineages, double **transprob)
{
    TransMatrixSwitch transmat_switch(states1.size(), states2.size());
    calc_transition_probs_switch(tree, last_tree, spr, mapping,
                                 states1, states2, model, lineages, 
                                 &transmat_switch);
    get_transition_probs_switch(&transmat_switch, transprob);
}


//=============================================================================
// prior for state space

void calc_state_priors(const States &states, const LineageCounts *lineages, 
                       const ArgModel *model, double *priors,
                       const int minage)
{
    const int nstates = states.size();
    const double *time_steps = model->time_steps;
    const double *popsizes = model->popsizes;
    const int *nbranches = lineages->nbranches;
    const int *ncoals = lineages->ncoals;
    
    // special case
    if (nstates == 0) {
        priors[0] = 1.0;
        return;
    }

    for (int i=0; i<nstates; i++) {
        int b = states[i].time;

        if (b < minage) {
            priors[i] = 0.0;
            continue;
        }

        double sum = 0.0;
        for (int m=minage; m<b; m++)
            sum += time_steps[m] * nbranches[m] / (2.0 * popsizes[m]);

        priors[i] = (1.0 - exp(- time_steps[b] * nbranches[b] /
                               (2.0 * popsizes[b]))) / ncoals[b] * exp(-sum);
    }
}


void calc_state_priors_log(const States &states, const LineageCounts *lineages, 
                       const ArgModel *model, double *priors,
                       const int minage)
{
    const int nstates = max((int) states.size(), 1);
    calc_state_priors(states, lineages, model, priors, minage);
    for (int i=0; i<nstates; i++)
        priors[i] = log(priors[i]);
}



//=============================================================================
// assert transition matrices

double calc_recomb_prob(const LocalTree *tree, const ArgModel *model)
{
    double treelen = get_treelen(tree, model->times, model->ntimes, false);
    return 1.0 - exp(-max(model->rho * treelen, model->rho));
}


bool assert_transmat(const LocalTree *tree, const ArgModel *model,
                     const TransMatrix *matrix)
{
    const int nleaves = tree->get_num_leaves();
    const int nnodes = tree->nnodes;
    const int nnodes2 = nnodes + 2;

    
    LineageCounts lineages(model->ntimes);
    States states;
    get_coal_states(tree, model->ntimes, states);
    int nstates = states.size();
    int displace[nnodes2];
    int newleaf = nleaves;

    // get local tree we can modify
    LocalTree tree2(tree->nnodes, tree->capacity + 2);
    tree2.copy(*tree);
    
    
    for (int i=0; i<nstates; i++) {
        for (int j=0; j<nstates; j++) {
            State state1 = states[i];
            State state2 = states[j];
            Spr spr;
            int sister_age = tree2.nodes[state1.node].age;
            
            add_tree_branch(&tree2, state1.node, state1.time);

            double p = -INFINITY;
            
            // recomb could be on new branch
            // recoal is state2
            spr.recomb_node = newleaf;
            spr.coal_node = state2.node;
            spr.coal_time = state2.time;

            // sum over possible recombination times
            for (int rtime=0; rtime<=min(state1.time, state2.time); rtime++) {
                spr.recomb_time = rtime;
                p = logadd(p, calc_spr_prob(model, &tree2, spr, lineages));
            }


            if (state1.node == state2.node) {
                // recomb could be state1.node
                // recoal is on new branch or parent of state1.node
                spr.recomb_node = state1.node;

                if (state2.time < state1.time) {
                    spr.coal_node = newleaf;
                } else if (state2.time >= state1.time) {
                    spr.coal_node = tree2.nodes[newleaf].parent;
                }
                spr.coal_time = state2.time;
                
                // sum over possible recombination times
                for (int rtime=sister_age; 
                     rtime<=min(state1.time, state2.time); rtime++) {
                    spr.recomb_time = rtime;
                    p = logadd(p, calc_spr_prob(model, &tree2, spr, lineages));
                }
            }

            // calculate probability of recombination
            double recomb_prob = calc_recomb_prob(&tree2, model);
            p += log(recomb_prob);

            if (i == j)
                p = logadd(p, log(1.0 - recomb_prob));

            
            double p2 = matrix->get_log(tree, states, i, j);

            // compare probabilities
            printf("> (%d,%d)->(%d,%d): %e = %e; %e\n", 
                       state1.node, state1.time,
                       state2.node, state2.time, p, p2, p - p2);
            if (!fequal(p, p2, 1e-4, 1e-9)) {
                return false;
            }

            remove_tree_branch(&tree2, newleaf, displace);
        }
    }

    return true;
}



bool assert_transmat_switch(const LocalTree *tree, const Spr &_spr,
                            const ArgModel *model, 
                            const TransMatrixSwitch *matrix)
{
    const int nleaves = tree->get_num_leaves();
    const int nnodes = tree->nnodes;
    const int ntimes = model->ntimes;

    // get local tree we can modify
    LocalTree tree1(tree->nnodes, tree->capacity + 2);
    tree1.copy(*tree);
    
    // build next tree from last tree
    LocalTree tree2(nnodes);
    tree2.copy(tree1);
    apply_spr(&tree2, _spr);

    // define mapping
    int mapping[nnodes];
    for (int i=0; i<nnodes; i++)
        mapping[i] = i;
    mapping[tree1.nodes[_spr.recomb_node].parent] = -1;

    // get states
    States states1, states2;
    get_coal_states(&tree1, ntimes, states1);
    get_coal_states(&tree2, ntimes, states2);
    int nstates1 = states1.size();
    int nstates2 = states2.size();

    // get lineages
    LineageCounts lineages(ntimes);

    // add/remove branch data
    int displace[nnodes + 2];
    int newleaf = nleaves;
    int displaced = nnodes;
    int newcoal = nnodes + 1;


    for (int i=0; i<nstates1; i++) {
        for (int j=0; j<nstates2; j++) {
            State state1 = states1[i];
            State state2 = states2[j];
            double p = -INFINITY;
            Spr spr = _spr;

            double p2 = matrix->get_log(i, j);
            if (p2 == -INFINITY)
                continue;

            // add new branch and modify spr
            add_tree_branch(&tree1, state1.node, state1.time);
            add_tree_branch(&tree2, state2.node, state2.time);        
            add_spr_branch(&tree2, &tree1, state2, state1,
                           &spr, mapping, newleaf, displaced, newcoal);
            
            // calculate probability of recombination
            double recomb_prob = calc_recomb_prob(&tree1, model);
            p = log(recomb_prob);
            p += calc_spr_prob(model, &tree1, spr, lineages);

            remove_tree_branch(&tree1, newleaf, displace);
            remove_tree_branch(&tree2, newleaf, displace);
            
            
            // compare probabilities
            printf("> %d,%d (%d,%d)->(%d,%d): %e = %e; %e\n", 
                   i, j,
                   state1.node, state1.time,
                   state2.node, state2.time, p, p2, p - p2);
            if (!fequal(p, p2, 1e-4, 1e-9)) {
                return false;
            }
        }
    }

    return true;
}




bool assert_transmat_internal(const LocalTree *tree, const ArgModel *model,
                              const TransMatrix *matrix)
{
    LineageCounts lineages(model->ntimes);
    States states;
    get_coal_states_internal(tree, model->ntimes, states);
    int nstates = states.size();
    int subtree_root = tree->nodes[tree->root].child[0];

    // get local tree we can modify
    LocalTree tree2(tree->nnodes);
    tree2.copy(*tree);
    
    
    for (int i=0; i<nstates; i++) {
        for (int j=0; j<nstates; j++) {
            State state1 = states[i];
            State state2 = states[j];
            Spr spr;
            int sister_age = tree2.nodes[state1.node].age;
            
            Spr add_spr(subtree_root, tree->nodes[subtree_root].age,
                        state1.node, state1.time);
            apply_spr(&tree2, add_spr);

            // calculate probability of recombination
            double recomb_prob = calc_recomb_prob(&tree2, model);

            double p = -INFINITY;

            // recomb could be on new branch
            // recoal is state2
            spr.recomb_node = subtree_root;
            spr.coal_node = state2.node;
            spr.coal_time = state2.time;
            assert(tree2[subtree_root].age <= spr.coal_time);

            // sum over possible recombination times
            for (int rtime=tree2[subtree_root].age; 
                 rtime<=min(state1.time, state2.time); rtime++) {
                spr.recomb_time = rtime;
                p = logadd(p, calc_spr_prob(model, &tree2, spr, lineages));
            }


            if (state1.node == state2.node) {
                // recomb could be state1.node
                // recoal is on new branch or parent of state1.node
                spr.recomb_node = state1.node;

                if (state2.time < state1.time) {
                    spr.coal_node = subtree_root;
                } else if (state2.time >= state1.time) {
                    spr.coal_node = tree2.nodes[subtree_root].parent;
                }
                spr.coal_time = state2.time;
                
                // sum over possible recombination times
                for (int rtime=sister_age; 
                     rtime<=min(state1.time, state2.time); rtime++) {
                    spr.recomb_time = rtime;
                    p = logadd(p, calc_spr_prob(model, &tree2, spr, lineages));
                }
            }

            p += log(recomb_prob);

            if (i == j)
                p = logadd(p, log(1.0 - recomb_prob));

            
            double p2 = matrix->get_log(tree, states, i, j);

            printf("> (%d,%d)->(%d,%d): %e = %e; %e\n", 
                   state1.node, state1.time,
                   state2.node, state2.time, p, p2, p - p2);
            if (!fequal(p, p2, 1e-4, 1e-9)) {
                return false;
            }

            Spr remove_spr(subtree_root, tree2[subtree_root].age, 
                           tree2.root, model->ntimes+1);
            apply_spr(&tree2, remove_spr);
        }
    }

    return true;
}


bool assert_transmat_switch_internal(const LocalTree *last_tree, 
                                     const LocalTree *tree, 
                                     const Spr &_spr, const int *_mapping, 
                                     ArgModel *model, 
                                     TransMatrixSwitch *transmat_switch)
{
    LineageCounts lineages(model->ntimes);
    States states, last_states;
    get_coal_states_internal(last_tree, model->ntimes, last_states);
    get_coal_states_internal(tree, model->ntimes, states);
    int nstates1 = last_states.size();
    int nstates2 = states.size();
    int last_subtree_root = last_tree->nodes[last_tree->root].child[0];
    int subtree_root = tree->nodes[tree->root].child[0];    
    
    for (int i=0; i<nstates1; i++) {
        for (int j=0; j<nstates2; j++) {
            State state1 = last_states[i];
            State state2 = states[j];

            double p2 = transmat_switch->get_log(i, j);
            if (p2 == -INFINITY)
                continue;

            // get local tree we can modify
            LocalTree last_tree2(last_tree->nnodes);
            last_tree2.copy(*last_tree);
            LocalTree tree2(tree->nnodes);
            tree2.copy(*tree);
            Spr spr = _spr;
            int mapping[tree->nnodes];
            for (int k=0; k<tree->nnodes; k++)
                mapping[k] = _mapping[k];


            // add branches
            Spr add_spr(last_subtree_root, last_tree2[last_subtree_root].age,
                        state1.node, state1.time);
            apply_spr(&last_tree2, add_spr);
            Spr add_spr2(subtree_root, tree2[subtree_root].age,
                         state2.node, state2.time);
            apply_spr(&tree2, add_spr2);
            add_spr_branch(&tree2, &last_tree2, state2, state1,
                           &spr, mapping, subtree_root, last_subtree_root);

            double recomb_prob = calc_recomb_prob(&last_tree2, model);
            double p = log(recomb_prob);
            p += calc_spr_prob(model, &last_tree2, spr, lineages);
            
            printf("> (%d,%d)->(%d,%d): %e = %e; %e\n", 
                   state1.node, state1.time,
                   state2.node, state2.time, p, p2, p - p2);
            if (!fequal(p, p2, 1e-4, 1e-9)) {
                return false;
            }
        }
    }    

    return true;
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
    make_node_mapping(last_ptree, nnodes, recomb_node, mapping);
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


bool arghmm_assert_transmat(int nnodes, int *ptree, int *ages, 
                            int ntimes, double *times, 
                            double *popsizes, double rho)
{
    ArgModel model(ntimes, times, popsizes, rho, 0.0);
    LocalTree tree(ptree, nnodes, ages);

    States states;
    get_coal_states(&tree, ntimes, states);

    LineageCounts lineages(ntimes);
    lineages.count(&tree);

    TransMatrix transmat(ntimes, states.size());
    calc_transition_probs(&tree, &model, states, &lineages, &transmat);
    
    return assert_transmat(&tree, &model, &transmat);
}



bool arghmm_assert_transmat_switch(
    int nnodes, int *ptree, int *ages, 
    int recomb_node, int recomb_time, int coal_node, int coal_time,
    int ntimes, double *times, 
    double *popsizes, double rho)
{
    ArgModel model(ntimes, times, popsizes, rho, 0.0);
    LocalTree last_tree(ptree, nnodes, ages);

    // build next tree from last tree
    LocalTree tree(nnodes);
    tree.copy(last_tree);
    Spr spr(recomb_node, recomb_time, coal_node, coal_time);
    apply_spr(&tree, spr);

    // define mapping
    int mapping[nnodes];
    for (int i=0; i<nnodes; i++)
        mapping[i] = i;
    mapping[last_tree.nodes[spr.recomb_node].parent] = -1;

    // get states
    States states, last_states;
    get_coal_states(&last_tree, ntimes, last_states);
    get_coal_states(&tree, ntimes, states);

    // get lineages
    LineageCounts lineages(ntimes);
    lineages.count(&last_tree);

    TransMatrixSwitch transmat_switch(last_states.size(), states.size());
    calc_transition_probs_switch(&tree, &last_tree, 
        spr, mapping, last_states, states, &model, &lineages, 
        &transmat_switch);

    return assert_transmat_switch(&last_tree, spr, &model, &transmat_switch);
}


bool arghmm_assert_transmat_internal(int nnodes, int *ptree, int *ages, 
                                     int ntimes, double *times, 
                                     double *popsizes, double rho)
{
    ArgModel model(ntimes, times, popsizes, rho, 0.0);
    LocalTree tree(ptree, nnodes, ages);

    // randomly choose branch to remove
    int node = irand(nnodes);
    Spr remove_spr(node, tree[node].age, tree.root, ntimes+1);
    apply_spr(&tree, remove_spr);
    int *c = tree[tree.root].child;
    if (c[0] != node) {
        int tmp = c[0];
        c[0] = c[1];
        c[1] = tmp;
    }


    States states;
    get_coal_states_internal(&tree, ntimes, states);

    LineageCounts lineages(ntimes);
    lineages.count(&tree, true);

    TransMatrix transmat(ntimes, states.size());
    calc_transition_probs_internal(&tree, &model, states, &lineages, &transmat);
    
    return assert_transmat_internal(&tree, &model, &transmat);
}


bool arghmm_assert_transmat_switch_internal(
    LocalTrees *trees, int ntimes, double *times, double *popsizes, double rho)
{
    ArgModel model(ntimes, times, popsizes, rho, 0.0);

    const int maxtime = model.ntimes + 1;
    int *removal_path = new int [trees->get_num_trees()];
    int node = irand(trees->nnodes);
    sample_arg_removal_path(trees, node, removal_path);
    remove_arg_thread_path(trees, removal_path, maxtime);

    // build matrices
    ArgHmmMatrices matrices;
    ArgHmmMatrixIter matrix_iter(&model, NULL, trees);
    matrix_iter.set_internal(true);

    LocalTree *last_tree = NULL;
    for (matrix_iter.begin(); matrix_iter.more(); matrix_iter.next()) {
        matrix_iter.get_matrices(&matrices);
        LocalTrees::iterator it = matrix_iter.get_tree_iter();
        LocalTree *tree = it->tree;

        if (matrices.transmat_switch) {
            if (!assert_transmat_switch_internal(last_tree, tree, it->spr,
                                                 it->mapping, &model, 
                                                 matrices.transmat_switch))
                return false;
        }
        
        last_tree = tree;
    }

    return true;
}



} // extern C

} // namespace arghmm




/*
//=============================================================================
// OLD CODE


double calc_recomb_recoal(
    const LocalTree *last_tree, const ArgModel *model, 
    const LineageCounts *lineages, 
    const Spr &spr, const State state1, const State state2, 
    const int recomb_parent_age)
{
    const int *nbranches = lineages->nbranches;
    const int *ncoals = lineages->ncoals;

    // get times
    const int a = state1.time;
    const int k = spr.recomb_time;
    const int j = spr.coal_time;
    
    double last_treelen = get_treelen_branch(
        last_tree, model->times, model->ntimes,
        state1.node, state1.time, -1, false);
    double last_treelen_b = last_treelen + get_basal_branch(
        last_tree, model->times, model->ntimes,
        state1.node, state1.time);
    
    // recomb prob
    int nbranches_k = lineages->nbranches[k] + int(k < a);
    int nrecombs_k = lineages->nrecombs[k] + int(k < a) + int(k == a);
    double recomb_prob = nbranches_k * model->time_steps[k] /
        (nrecombs_k * last_treelen_b) * 
        (1.0 - exp(-max(model->rho * last_treelen, model->rho)));

    // coal prob
    double sum = 0.0;
    for (int m=k; m<j; m++) {
        int nbranches_m = lineages->nbranches[m] 
            - int(m < recomb_parent_age) + int(m < a);
        sum += model->time_steps[m] * nbranches_m / (2.0 * model->popsizes[m]);
    }
    int nbranches_j = nbranches[j] - int(j < recomb_parent_age) + int(j < a);
    int ncoals_j = ncoals[j] - int(j <= recomb_parent_age)
        - int(j == recomb_parent_age) + int(j <= a) + int (j == a);
    
    double p = recomb_prob * exp(-sum) / ncoals_j;
    if (j < model->ntimes - 1)
        p *= 1.0 - exp(-model->time_steps[j] * nbranches_j / 
                       (2.0*model->popsizes[j]));

    if (ncoals_j <= 0 || nbranches_j <= 0 || 
        nrecombs_k <= 0 || nbranches_k <= 0) {
        printf("counts %d %d %d %d %e\n", 
               ncoals_j, nbranches_j, nrecombs_k, nbranches_k, p);
        assert(false);
    }

    return p;
}
*/

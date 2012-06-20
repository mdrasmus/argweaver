
#include "local_tree.h"
#include "matrices.h"

namespace arghmm {

using namespace std;


//=============================================================================
// Sample recombinations


double recomb_prob_unnormalized(ArgModel *model, LocalTree *tree, 
                                const LineageCounts &lineages, 
                                const State &last_state, 
                                const State &state,
                                NodePoint &recomb)
{

    const int k = recomb.time;
    const int j = state.time;

    int nbranches_k = lineages.nbranches[k]
        + int(k < last_state.time);
    int nrecombs_k = lineages.nrecombs[k]
        + int(k <= last_state.time)
        + int(k == last_state.time);

    double sum = 0.0;
    int recomb_parent_age = (recomb.node == -1 || 
                             (recomb.node == last_state.node &&
                              recomb.time < last_state.time)) ? 
        last_state.time :
        tree->nodes[tree->nodes[recomb.node].parent].age;

    for (int m=k; m<j; m++) {
        int nbranches_m = lineages.nbranches[m] 
            + int(m < last_state.time) 
            - int(m < recomb_parent_age);
        sum += (model->time_steps[m] * nbranches_m 
                / (2.0 * model->popsizes[m]));
    }

    return (nbranches_k * model->time_steps[k] / nrecombs_k) * exp(- sum);
}


void sample_recombinations(
    LocalTrees *trees, ArgModel *model, ArgHmmMatrixList *matrix_list,
    int *thread_path, vector<int> &recomb_pos, vector<NodePoint> &recombs)
{
    States states;
    LineageCounts lineages(model->ntimes);
    const int new_node = -1;
    vector <NodePoint> candidates;
    vector <double> probs;


    // loop through local blocks
    int blocki = 0;
    int end = trees->start_coord;
    for (LocalTrees::iterator it=trees->begin(); 
         it != trees->end(); ++it, blocki++) {

        // get local block information
        int start = end + 1;  // don't allow new recomb at start
        end = start - 1 + it->blocklen;
        ArgHmmMatrices matrices = matrix_list->matrices[blocki];
        LocalTree *tree = it->tree;
        lineages.count(tree);
        get_coal_states(tree, model->ntimes, states);
        int next_recomb = -1;

        
        // loop through positions in block
        for (int i=start; i<end; i++) {
            
            if (thread_path[i] == thread_path[i-1]) {
                // no change in state, recombination is optional
                if (i > next_recomb) {
                    // sample the next recomb pos
                    int last_state = thread_path[i-1];

                    double rate;
                    if (matrices.transmat_compress) {
                        TransMatrixCompress *m = matrices.transmat_compress;
                        int a = states[last_state].time;
                        double self_trans = exp(m->get_transition_prob(
                            tree, states, last_state, last_state));
                        rate = 1.0 - (m->norecombs[a] / self_trans);
                        
                    } else {
                        // TODO: need to compute row sum
                        assert(false);
                        
                        double treelen_b = get_treelen(
                            tree, model->times, model->ntimes);
                        double treelen = treelen_b - get_basal_branch(
                            tree, model->times, model->ntimes, -1, -1);

                        double treelen2_b = get_treelen_branch(
                            tree, model->times, model->ntimes,
                            states[last_state].node,
                            states[last_state].time, treelen_b);
                        double treelen2 = treelen2_b - get_basal_branch(
                            tree, model->times, model->ntimes, 
                            states[last_state].node, states[last_state].time);
                        double self_trans = matrices.transmat[last_state][last_state];
                        rate = max(1.0 - exp(-model->rho * (treelen2 - treelen)
                                             - self_trans), model->rho);
                    }
                        

                    // NOTE: the min prevents large floats from overflowing
                    // when cast to int
                    next_recomb = int(min(float(end), i + expovariate(rate)));
                }

                if (i < next_recomb)
                    continue;
            }

            // sample recombination
            next_recomb = -1;
            State state = states[thread_path[i]];
            State last_state = states[thread_path[i-1]];
            
            // there must be a recombination
            // either because state changed or we choose to recombine
            // find candidates
            candidates.clear();
            int end_time = min(state.time, last_state.time);
            if (state.node == last_state.node) {
                // y = v, k in [0, min(timei, last_timei)]
                // y = node, k in Sr(node)
                for (int k=tree->nodes[state.node].age; k<=end_time; k++)
                    candidates.push_back(NodePoint(state.node, k));
            }
            
            for (int k=0; k<=end_time; k++)
                candidates.push_back(NodePoint(new_node, k));

            // compute probability of each candidate
            probs.clear();
            for (vector<NodePoint>::iterator it=candidates.begin(); 
                 it != candidates.end(); ++it) {
                probs.push_back(recomb_prob_unnormalized(
                    model, tree, lineages, last_state, state, *it));
            }

            // sample recombination
            recomb_pos.push_back(i);
            recombs.push_back(candidates[sample(&probs[0], probs.size())]);

            assert(recombs[recombs.size()-1].time <= min(state.time,
                                                         last_state.time));
        }
    }
}



void max_recombinations(
    LocalTrees *trees, ArgModel *model, ArgHmmMatrixList *matrix_list,
    int *thread_path, vector<int> &recomb_pos, vector<NodePoint> &recombs)
{
    States states;
    LineageCounts lineages(model->ntimes);
    const int new_node = -1;
    vector <NodePoint> candidates;
    vector <double> probs;


    // loop through local blocks
    int blocki = 0;
    int end = trees->start_coord;
    for (LocalTrees::iterator it=trees->begin(); 
         it != trees->end(); ++it, blocki++) {

        // get local block information
        int start = end + 1;  // don't allow new recomb at start
        end = start - 1 + it->blocklen;
        ArgHmmMatrices matrices = matrix_list->matrices[blocki];
        LocalTree *tree = it->tree;
        lineages.count(tree);
        get_coal_states(tree, model->ntimes, states);
        
        // loop through positions in block
        for (int i=start; i<end; i++) {
            // its never worth sampling a recombination if you don't have to
            if (thread_path[i] == thread_path[i-1])
                continue;
            
            State state = states[thread_path[i]];
            State last_state = states[thread_path[i-1]];
            
            // find candidate recombinations
            candidates.clear();
            int end_time = min(state.time, last_state.time);
            if (state.node == last_state.node) {
                // y = node, k in [node.age, min(timei, last_timei)]
                for (int k=tree->nodes[state.node].age; k<=end_time; k++)
                    candidates.push_back(NodePoint(state.node, k));
            }
            
            // y = v, k in [0, min(timei, last_timei)]
            for (int k=0; k<=end_time; k++)
                candidates.push_back(NodePoint(new_node, k));

            // compute probability of each candidate
            probs.clear();
            for (vector<NodePoint>::iterator it=candidates.begin(); 
                 it != candidates.end(); ++it) {
                probs.push_back(recomb_prob_unnormalized(
                    model, tree, lineages, last_state, state, *it));
            }

            // find max prob recombination
            recomb_pos.push_back(i);
            int argmax = 0;
            double maxprob = probs[0];
            for (unsigned int m=1; m<probs.size(); m++) {
                if (probs[m] > maxprob) {
                    argmax = m;
                    maxprob = probs[m];
                }
            }
            recombs.push_back(candidates[argmax]);

            assert(recombs[recombs.size()-1].time <= min(state.time,
                                                         last_state.time));
        }
    }
}


} // namespace arghmm

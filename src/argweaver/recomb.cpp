
#include "local_tree.h"
#include "matrices.h"

namespace argweaver {

using namespace std;


//=============================================================================
// Sample recombinations


double recomb_prob_unnormalized(const ArgModel *model, const LocalTree *tree,
                                const LineageCounts &lineages,
                                const State &last_state,
                                const State &state,
                                const NodePoint &recomb, bool internal)
{
    const int k = recomb.time;
    const int j = state.time;

    int root_time;
    int recomb_parent_age;
    if (internal) {
        int subtree_root = tree->nodes[tree->root].child[0];
        int maintree_root = tree->nodes[tree->root].child[1];
        root_time = max(tree->nodes[maintree_root].age, last_state.time);
        recomb_parent_age = (recomb.node == subtree_root ||
                             tree->nodes[recomb.node].parent == -1 ||
                             recomb.node == last_state.node) ?
            last_state.time :
            tree->nodes[tree->nodes[recomb.node].parent].age;
    } else {
        root_time = max(tree->nodes[tree->root].age, last_state.time);
        recomb_parent_age = (recomb.node == -1 ||
                             tree->nodes[recomb.node].parent == -1 ||
                             recomb.node == last_state.node) ?
            last_state.time :
            tree->nodes[tree->nodes[recomb.node].parent].age;
    }

    int nbranches_k = lineages.nbranches[k]
        + int(k < last_state.time);
    int nrecombs_k = lineages.nrecombs[k]
        + int(k <= last_state.time)
        + int(k == last_state.time)
        - int(k == root_time);



    // The commented-out section would be correct if we were rounding
    // recombs to nearest time point in the same way as coals (keep because
    // hope to make this change eventually)
    /*double precomb = nbranches_k*model->coal_time_steps[2*k]/nrecombs_k;
    if (k > 0) {
        int m=k-1;
        int nbranches_m = lineages.nbranches[m]
            + int(m < last_state.time)
            - int(m < recomb_parent_age);
        int nrecombs_m = lineages.nrecombs[m]
            + int(m <= last_state.time)
            + int(m == last_state.time);
        precomb += nbranches_m*model->coal_time_steps[2*k-1]/nrecombs_m;
    }*/
    double precomb = nbranches_k * model->time_steps[k] / nrecombs_k;

    // probability of not coalescing before time j-1
    double sum = 0.0;
    for (int m=k; m<j-1; m++) {
        int nbranches_m = lineages.nbranches[m]
            + int(m < last_state.time)
            - int(m < recomb_parent_age);
        sum += (model->time_steps[m] * nbranches_m
                / (2.0 * model->popsizes[m]));
    }

    // probability of coalescing at time j
    double pcoal = 1.0;
    int nbranches_j = lineages.nbranches[j]
        + int(j < last_state.time)
        - int(j < recomb_parent_age);
    assert(nbranches_j > 0);
    if (k == j) {
        // in this case coalesce event must happen between time
        // j,j+1/2 (coal_time_step[2j])
        // NOTE: if k == ntimes - 2, coalescence is probability 1.0
        if (j < model->ntimes - 2) {
            pcoal = 1.0 - exp(-model->coal_time_steps[2*j] * nbranches_j
                              / (2.0 * model->popsizes[j]));
        }
    } else {
        // otherwise it could happen anytime between j-1/2 and j+1/2
        int m = j - 1;
        int nbranches_m = lineages.nbranches[m]
            + int(m < last_state.time)
            - int(m < recomb_parent_age);
        // have to not coalesce in the first half interval before j
        sum += (model->coal_time_steps[2*m] * nbranches_m
                / (2.0 * model->popsizes[m]));

        // NOTE: if k == ntimes - 2, coalescence is probability 1.0
        if (j < model->ntimes - 2) {
            pcoal = 1.0 - exp(
                -model->coal_time_steps[2*j-1] * nbranches_m
                / (2.0 * model->popsizes[m])
                - model->coal_time_steps[2*j] * nbranches_j
                / (2.0 * model->popsizes[j]));
        }
    }

    // probability of recoalescing on a choosen branch
    int ncoals_j = lineages.ncoals[j]
        - int(j <= recomb_parent_age) - int(j == recomb_parent_age)
        + int(j <= last_state.time) + int(j == last_state.time);
    pcoal /= ncoals_j;

    return precomb * exp(- sum) * pcoal;
}


// Returns the possible recombination events that are compatiable with
// the transition (last_state -> state).
void get_possible_recomb(const LocalTree *tree,
                         const State last_state, const State state,
                         bool internal, vector<NodePoint> &candidates)
{
    // represents the branch above the new node in the tree.
    const int new_node = -1;

    int end_time = min(state.time, last_state.time);
    if (state.node == last_state.node) {
        // y = v, k in [0, min(timei, last_timei)]
        // y = node, k in Sr(node)
        for (int k=tree->nodes[state.node].age; k<=end_time; k++)
            candidates.push_back(NodePoint(state.node, k));
    }

    if (internal) {
        const int subtree_root = tree->nodes[tree->root].child[0];
        const int subtree_root_age = tree->nodes[subtree_root].age;
        for (int k=subtree_root_age; k<=end_time; k++)
            candidates.push_back(NodePoint(subtree_root, k));
    } else {
        for (int k=0; k<=end_time; k++)
            candidates.push_back(NodePoint(new_node, k));
    }
}


void sample_recombinations(
    const LocalTrees *trees, const ArgModel *model,
    ArgHmmMatrixIter *matrix_iter,
    int *thread_path, vector<int> &recomb_pos, vector<NodePoint> &recombs,
    bool internal)
{
    States states;
    LineageCounts lineages(model->ntimes);
    vector <NodePoint> candidates;
    vector <double> probs;

    // loop through local blocks
    for (matrix_iter->begin(); matrix_iter->more(); matrix_iter->next()) {

        // get local block information
        ArgHmmMatrices &matrices = matrix_iter->ref_matrices();
        LocalTree *tree = matrix_iter->get_tree_spr()->tree;
        lineages.count(tree, internal);
        matrices.states_model.get_coal_states(tree, states);
        int next_recomb = -1;

        // don't sample recombination if there is no state space
        if (internal && states.size() == 0)
            continue;

        int start = matrix_iter->get_block_start();
        int end = matrix_iter->get_block_end();
        if (matrices.transmat_switch || start == trees->start_coord) {
            // don't allow new recomb at start if we are switching blocks
            start++;
        }

        // loop through positions in block
        for (int i=start; i<end; i++) {

            if (thread_path[i] == thread_path[i-1]) {
                // no change in state, recombination is optional

                if (i > next_recomb) {
                    // sample the next recomb pos
                    int last_state = thread_path[i-1];
                    TransMatrix *m = matrices.transmat;
                    int a = states[last_state].time;
                    double self_trans = m->get(
                        tree, states, last_state, last_state);
                    double rate = 1.0 - (m->norecombs[a] / self_trans);

                    // NOTE: the min prevents large floats from overflowing
                    // when cast to int
                    next_recomb = int(min(double(end), i + expovariate(rate)));
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
            get_possible_recomb(tree, last_state, state, internal, candidates);

            // compute probability of each candidate
            probs.clear();
            for (vector<NodePoint>::iterator it=candidates.begin();
                 it != candidates.end(); ++it) {
                probs.push_back(recomb_prob_unnormalized(
                    model, tree, lineages, last_state, state, *it, internal));
            }

            // sample recombination
            recomb_pos.push_back(i);
            recombs.push_back(candidates[sample(&probs[0], probs.size())]);

            assert(recombs[recombs.size()-1].time <= min(state.time,
                                                         last_state.time));
        }
    }
}


} // namespace argweaver

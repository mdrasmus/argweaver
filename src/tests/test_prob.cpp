#include "gtest/gtest.h"

#include "argweaver/common.h"
#include "argweaver/local_tree.h"
#include "argweaver/model.h"
#include "argweaver/recomb.h"
#include "argweaver/states.h"
#include "argweaver/thread.h"
#include "argweaver/total_prob.h"


namespace argweaver {

// The probability of all possible SPRs in the DSMC model should add up to one.
TEST(ProbTest, test_calc_spr_prob)
{
    // Setup model.
    int ntimes = 5;
    double times[] = {0, 10, 20, 30, 40};
    double rho = 1e-9;
    double mu = 2.5e-9;
    double popsize = 1e4;
    ArgModel model(ntimes, times, NULL, rho, mu);
    model.set_popsizes(popsize, ntimes);

    // Read tree.
    const char *newick =
        "((0,1)5[&&NHX:age=10],((2,3)6[&&NHX:age=20],4)7[&&NHX:age=20])8[&&NHX:age=30]";
    LocalTree tree;
    parse_local_tree(newick, &tree, times, ntimes);

    LineageCounts lineages(ntimes);
    States states;
    get_coal_states_external(&tree, ntimes, states);

    double prob = -INFINITY;
    for (unsigned int i=0; i<states.size(); i++) {
        const State recomb = states[i];
        for (unsigned int j=0; j<states.size(); j++) {
            const State coal = states[j];

            // Ignore nonsense SPRs
            // - re-coalescence must be older than recombination
            // - re-coalescence must not be the same branch as recombination,
            //   bubbles are not allowed in the model.
            // - recombination cannot occur on the basal branch.
            // - re-coalescence immediately above broken node is
            //   redundant with re-coalescence with sibling at same time.
            if (coal.time < recomb.time ||
                coal.node == recomb.node ||
                recomb.node == tree.root)
                continue;
            const int broken = tree.nodes[recomb.node].parent;
            if (coal.node == broken && coal.time == tree.nodes[broken].age)
                continue;

            const Spr spr(recomb.node, recomb.time, coal.node, coal.time);
            const double p = calc_spr_prob(&model, &tree, spr, lineages);
            prob = logadd(prob, p);
        }
    }
    EXPECT_NEAR(prob, 0.0, 1e-5);
}


// The probability of all recombinations between two given trees should be
// proportional to the SPR probabilities.
TEST(ProbTest, test_sample_recomb_external)
{
    // Setup model.
    int ntimes = 5;
    double times[] = {0, 10, 20, 30, 40};
    double rho = 1e-9;
    double mu = 2.5e-9;
    double popsize = 1e4;
    ArgModel model(ntimes, times, NULL, rho, mu);
    model.set_popsizes(popsize, ntimes);
    bool internal = false;

    // Read tree.
    const char *newick =
        "((0,1)5[&&NHX:age=10],((2,3)6[&&NHX:age=20],4)7[&&NHX:age=20])8[&&NHX:age=30]";
    LocalTree tree;
    parse_local_tree(newick, &tree, times, ntimes);

    // Create a copy of the tree that will explicitly have a new branch and SPR.
    LocalTree tree2(0, tree.capacity + 2);
    tree2.copy(tree);

    // Get thread states.
    LineageCounts lineages(ntimes);
    LineageCounts lineages2(ntimes);
    lineages.count(&tree, internal);
    States states;
    get_coal_states(&tree, ntimes, states, internal);
    int nstates = states.size();

    // add/remove branch data
    int nleaves = tree.get_num_leaves();
    int newleaf = nleaves;
    int displaced = tree.nnodes;

    // Loop over all possible thread transitions (state1 -> state2).
    for (int i=0; i<nstates; i++) {
        for (int j=0; j<nstates; j++) {
            State state1 = states[i];
            State state2 = states[j];

            // Construct tree2 (explicitly add a branch at state1).
            add_tree_branch(&tree2, state1.node, state1.time);
            lineages2.count(&tree2, internal);

            // Get all candidate recombination points.
            vector<NodePoint> candidates;
            get_possible_recomb(&tree, state1, state2, internal, candidates);

            // Compute probability of each candidate recombination point.
            vector<double> probs;
            vector<double> probs2;
            for (vector<NodePoint>::iterator it=candidates.begin();
                 it != candidates.end(); ++it) {
                // Construct SPR representing this recombination.
                Spr spr(it->node, it->time, state2.node, state2.time);

                // Adjust recomb node to account for new nodes/displacement.
                if (spr.recomb_node == -1)
                    // Recombination is on new branch.
                    spr.recomb_node = newleaf;
                else if (spr.recomb_node == nleaves)
                    // Recombination was on displaced branch.
                    spr.recomb_node = displaced;
                if (spr.coal_node == nleaves)
                    /// Recoal was on displaced branch.
                    spr.coal_node = displaced;

                // Adjust re-coalescence.
                if (state1.node == state2.node) {
                    // 1. recomb is on newleaf
                    //    recoal is on state1 or parent of state1
                    // 2. recomb is on state1.node
                    //    recoal is on newleaf or parent of state1.node
                    if (state2.time < state1.time) {
                        spr.coal_node = newleaf;
                    } else if (state2.time >= state1.time) {
                        spr.coal_node = tree2.nodes[newleaf].parent;
                    }
                }

                probs.push_back(recomb_prob_unnormalized(
                    &model, &tree, lineages, state1, state2, *it, internal));
                probs2.push_back(exp(calc_spr_prob(&model, &tree2, spr,
                                                   lineages2)));
            }

            // Revert changes to tree2.
            remove_tree_branch(&tree2, newleaf);

            // Normalize probabilities.
            double total = 0.0, total2 = 0.0;
            for (unsigned int k=0; k<probs.size(); k++)
                total += probs[k];
            for (unsigned int k=0; k<probs2.size(); k++)
                total2 += probs2[k];

            // Assert that normalized probabilities are equal.
            for (unsigned int k=0; k<probs.size(); k++) {
                double p = probs[k] / total;
                double p2 = probs2[k] / total2;
                EXPECT_NEAR(p, p2, 1e-3);
            }
        }
    }
}


// The probability of all recombinations between two given trees should be
// proportional to the SPR probabilities.
TEST(ProbTest, test_sample_recomb_internal)
{
    // Setup model.
    int ntimes = 6;
    double times[] = {0, 10, 20, 30, 40, 50};
    double rho = 1e-9;
    double mu = 2.5e-9;
    double popsize = 1e4;
    ArgModel model(ntimes, times, NULL, rho, mu);
    model.set_popsizes(popsize, ntimes);
    bool internal = true;

    // Read tree.
    const char *newick =
        "((0,1)5[&&NHX:age=10],((2,3)6[&&NHX:age=20],4)7[&&NHX:age=20])8[&&NHX:age=30]";
    LocalTree tree;
    parse_local_tree(newick, &tree, times, ntimes);

    // Make tree a partial tree by setting root age above valid range.
    tree.nodes[tree.root].age = model.get_removed_root_time();

    // Create a copy of the tree that will explicitly have a new branch and SPR.
    LocalTree tree2(tree);
    int subtree_root = tree2.nodes[tree2.root].child[0];

    // Get thread states.
    LineageCounts lineages(ntimes);
    LineageCounts lineages2(ntimes);
    lineages.count(&tree, internal);
    States states;
    get_coal_states(&tree, ntimes, states, internal);
    int nstates = states.size();

    // Loop over all possible thread transitions (state1 -> state2).
    for (int i=0; i<nstates; i++) {
        for (int j=0; j<nstates; j++) {
            State state1 = states[i];
            State state2 = states[j];

            // Construct tree2 explicitly using state1.
            Spr add_spr(subtree_root, tree2.nodes[subtree_root].age,
                        state1.node, state1.time);
            apply_spr(&tree2, add_spr);
            lineages2.count(&tree2);

            // Get all candidate recombination points.
            vector<NodePoint> candidates;
            get_possible_recomb(&tree, state1, state2, internal, candidates);

            // Compute probability of each candidate recombination point.
            vector<double> probs;
            vector<double> probs2;
            for (vector<NodePoint>::iterator it=candidates.begin();
                 it != candidates.end(); ++it) {
                // Construct SPR representing this recombination.
                Spr spr(it->node, it->time, state2.node, state2.time);

                // Adjust re-coalescence.
                if (state1.node == state2.node) {
                    // 1. recomb is on subtree_root
                    //    recoal is on state1 or parent of state1
                    // 2. recomb is on state1.node
                    //    recoal is on subtree_root or parent of state1.node
                    if (state2.time < state1.time) {
                        spr.coal_node = subtree_root;
                    } else if (state2.time >= state1.time) {
                        spr.coal_node = tree2.nodes[subtree_root].parent;
                    }
                }

                probs.push_back(recomb_prob_unnormalized(
                    &model, &tree, lineages, state1, state2, *it, internal));
                probs2.push_back(exp(calc_spr_prob(&model, &tree2, spr,
                                                   lineages2)));
            }

            // Revert changes to tree2.
            Spr remove_spr(subtree_root, tree2[subtree_root].age,
                           tree2.root, model.get_removed_root_time());
            apply_spr(&tree2, remove_spr);

            // Normalize probabilities.
            double total = 0.0, total2 = 0.0;
            for (unsigned int k=0; k<probs.size(); k++)
                total += probs[k];
            for (unsigned int k=0; k<probs2.size(); k++)
                total2 += probs2[k];

            // Assert that normalized probabilities are equal.
            for (unsigned int k=0; k<probs.size(); k++) {
                double p = probs[k] / total;
                double p2 = probs2[k] / total2;
                EXPECT_NEAR(p, p2, 1e-3);
            }
        }
    }
}


}  // namespace

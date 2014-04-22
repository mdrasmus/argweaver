#include "gtest/gtest.h"

#include "common.h"
#include "local_tree.h"
#include "model.h"
#include "states.h"
#include "total_prob.h"


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


}  // namespace

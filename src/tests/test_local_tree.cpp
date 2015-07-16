#include "gtest/gtest.h"

#include "argweaver/local_tree.h"


namespace argweaver {


// Parse a local tree from newick.
TEST(LocalTreeTest, parse_local_tree)
{
    const char *newick = "((0[&&NHX:age=0],1[&&NHX:age=0])3[&&NHX:age=10],2[&&NHX:age=0])4[&&NHX:age=20]";
    int ntimes = 5;
    double times[] = {0, 10, 20, 30, 40};

    LocalTree tree;
    bool result = parse_local_tree(newick, &tree, times, ntimes);

    // Assert parse.
    EXPECT_EQ(result, true);

    // Assert ages.
    EXPECT_EQ(tree.nodes[0].age, 0);
    EXPECT_EQ(tree.nodes[1].age, 0);
    EXPECT_EQ(tree.nodes[2].age, 0);
    EXPECT_EQ(tree.nodes[3].age, 1);
    EXPECT_EQ(tree.nodes[4].age, 2);

    // Assert structure.
    EXPECT_EQ(tree.nnodes, 5);
    EXPECT_EQ(tree.root, 4);
    EXPECT_EQ(tree.nodes[0].parent, 3);
    EXPECT_EQ(tree.nodes[1].parent, 3);
    EXPECT_EQ(tree.nodes[2].parent, 4);
    EXPECT_EQ(tree.nodes[3].parent, 4);
}


// Parse a local tree from newick with implied age of 0 for leaves.
TEST(LocalTreeTest, parse_local_tree_implied_leaf_ages)
{
    const char *newick = "((0,1)3[&&NHX:age=10],2)4[&&NHX:age=20]";
    int ntimes = 5;
    double times[] = {0, 10, 20, 30, 40};

    LocalTree tree;
    bool result = parse_local_tree(newick, &tree, times, ntimes);

    // Assert parse.
    EXPECT_EQ(result, true);

    // Assert ages.
    EXPECT_EQ(tree.nodes[0].age, 0);
    EXPECT_EQ(tree.nodes[1].age, 0);
    EXPECT_EQ(tree.nodes[2].age, 0);
    EXPECT_EQ(tree.nodes[3].age, 1);
    EXPECT_EQ(tree.nodes[4].age, 2);
}


// Parse a local tree from newick with extra meta data.
TEST(LocalTreeTest, parse_local_tree_extra_meta_data)
{
    const char *newick = "((0,1[&&NHX:])3[&&NHX:age=10:key=value],2)4[&&NHX:key=value:age=20]";
    int ntimes = 5;
    double times[] = {0, 10, 20, 30, 40};

    LocalTree tree;
    bool result = parse_local_tree(newick, &tree, times, ntimes);

    // Assert parse.
    EXPECT_EQ(result, true);

    // Assert ages.
    EXPECT_EQ(tree.nodes[0].age, 0);
    EXPECT_EQ(tree.nodes[1].age, 0);
    EXPECT_EQ(tree.nodes[2].age, 0);
    EXPECT_EQ(tree.nodes[3].age, 1);
    EXPECT_EQ(tree.nodes[4].age, 2);
}


}  // namespace

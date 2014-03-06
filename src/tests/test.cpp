#include "gtest/gtest.h"

namespace {

TEST(ExampleTest, test1) {
    EXPECT_EQ(1, 1);
}

TEST(ExampleTest, test2) {
    EXPECT_EQ(1, 1);
}


}  // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

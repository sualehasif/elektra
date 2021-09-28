#include "gtest/gtest.h"
#include "tests/test_parallel_euler_tour_tree.h"
#include "tests/test_parallel_skip_list.h"

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

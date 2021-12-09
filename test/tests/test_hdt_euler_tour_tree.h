#pragma once

#include "../../elektra/parallel_euler_tour_tree/hdt_euler_tour_tree.h"

#include "gtest/gtest.h"
#include <parlay/sequence.h>

namespace elektra::testing {
namespace hdt_euler_tour_tree {

using EulerTourTree = parallel_euler_tour_tree::HdtEulerTourTree;

TEST(HdtEulerTourTree, ComponentSize) {
  const int n = 5;
  EulerTourTree ett = EulerTourTree{n};

  // 0 - 1   4
  // | \
  // 2   3 - 5
  parlay::sequence<std::pair<int, int>> links{{
    {0, 1},
    {0, 2},
    {0, 3},
    {5, 3},
  }};
  ett.BatchLink(links, false);
  EXPECT_EQ(ett.ComponentSize(3), 5);
  EXPECT_EQ(ett.ComponentSize(4), 1);

  ett.Cut(0, 3);
  EXPECT_EQ(ett.ComponentSize(2), 3);
  EXPECT_EQ(ett.ComponentSize(3), 2);
}

}  // namespace hdt_element
}  // namespace elektra::testing

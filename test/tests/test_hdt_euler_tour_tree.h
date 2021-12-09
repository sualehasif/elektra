#pragma once

#include "../../elektra/parallel_euler_tour_tree/hdt_euler_tour_tree.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <parlay/sequence.h>

namespace elektra::testing {
namespace hdt_euler_tour_tree {

using EulerTourTree = parallel_euler_tour_tree::HdtEulerTourTree;

using ::testing::IsEmpty;
using ::testing::Pair;
using ::testing::UnorderedElementsAre;

TEST(HdtEulerTourTree, ComponentSize) {
  const int n = 6;
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

TEST(HdtEulerTourTree, GetAndClearLevelIEdges) {
  const int n = 7;
  EulerTourTree ett = EulerTourTree{n};

  // 0 - 1   4   6
  // | \     |
  // 2   3   5
  parlay::sequence<std::pair<int, int>> links{{
    {0, 1},
    {0, 2},
    {0, 3},
    {4, 5},
  }};
  parlay::sequence<bool> is_level_i_edge{{true, false, true, true}};
  ett.BatchLink(links, is_level_i_edge);

  EXPECT_THAT(ett.GetAndClearLevelIEdges(3), UnorderedElementsAre(Pair(0, 1), Pair(0, 3)));
  EXPECT_THAT(ett.GetAndClearLevelIEdges(3), IsEmpty());
  EXPECT_THAT(ett.GetAndClearLevelIEdges(4), UnorderedElementsAre(Pair(4, 5)));
  EXPECT_THAT(ett.GetAndClearLevelIEdges(4), IsEmpty());
  EXPECT_THAT(ett.GetAndClearLevelIEdges(6), IsEmpty());
}

}  // namespace hdt_element
}  // namespace elektra::testing

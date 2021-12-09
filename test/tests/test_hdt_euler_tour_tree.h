#pragma once

#include "../../elektra/parallel_euler_tour_tree/hdt_euler_tour_tree.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <parlay/sequence.h>

namespace elektra::testing {
namespace hdt_euler_tour_tree {

using EulerTourTree = parallel_euler_tour_tree::HdtEulerTourTree;

using ::testing::_;
using ::testing::InSequence;
using ::testing::IsEmpty;
using ::testing::MockFunction;
using ::testing::Pair;
using ::testing::UnorderedElementsAre;

using MockForEachIncidentVertexCallback = MockFunction<void(uint32_t, uint32_t, uint32_t)>;

TEST(HdtEulerTourTree, ComponentSize_Singleton) {
  EulerTourTree ett = EulerTourTree{1};
  EXPECT_EQ(ett.ComponentSize(0), 1);
}

TEST(HdtEulerTourTree, ComponentSize) {
  const int n = 5;
  EulerTourTree ett = EulerTourTree{n};

  // 0 - 1
  // | \
  // 2   3 - 4
  parlay::sequence<std::pair<int, int>> links{{
    {0, 1},
    {0, 2},
    {0, 3},
    {4, 3},
  }};
  ett.BatchLink(links, false);
  EXPECT_EQ(ett.ComponentSize(3), 5);

  ett.Cut(0, 3);
  EXPECT_EQ(ett.ComponentSize(2), 3);
  EXPECT_EQ(ett.ComponentSize(3), 2);
}

TEST(HdtEulerTourTree, GetAndClearLevelIEdges_Singleton) {
  EulerTourTree ett = EulerTourTree{1};
  EXPECT_THAT(ett.GetAndClearLevelIEdges(0), IsEmpty());
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

TEST(HdtEulerTourTree, FindNonTreeEdges_SingletonWithoutEdges) {
  EulerTourTree ett = EulerTourTree{1};
  const auto finder = ett.CreateNontreeEdgeFinder(0);

  MockForEachIncidentVertexCallback mock_callback;
  EXPECT_CALL(mock_callback, Call(_, _, _)).Times(0);

  finder.ForEachIncidentVertex(mock_callback.AsStdFunction(), 0, 100);
  finder.ForEachIncidentVertex(mock_callback.AsStdFunction(), 20, 40);
}

TEST(HdtEulerTourTree, FindNonTreeEdges_SingletonWithEdges) {
  EulerTourTree ett = EulerTourTree{1};
  ett.UpdateNontreeEdgeCounts({0}, parlay::sequence<int>(1, 50));
  const auto finder = ett.CreateNontreeEdgeFinder(0);

  InSequence s;
  MockForEachIncidentVertexCallback mock_callback;
  EXPECT_CALL(mock_callback, Call(0, 0, 50));
  EXPECT_CALL(mock_callback, Call(0, 20, 40));

  finder.ForEachIncidentVertex(mock_callback.AsStdFunction(), 0, 100);
  finder.ForEachIncidentVertex(mock_callback.AsStdFunction(), 20, 40);
}

}  // namespace hdt_element
}  // namespace elektra::testing

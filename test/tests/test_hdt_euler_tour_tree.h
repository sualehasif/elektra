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
  // Check ComponentSize() works on a singleton vertex.
  EulerTourTree ett = EulerTourTree{1};
  EXPECT_EQ(ett.ComponentSize(0), 1);
  EXPECT_THAT(ett.ComponentVertices(0), UnorderedElementsAre(0));
  EXPECT_THAT(ett.ComponentEdges(0), IsEmpty());
}

TEST(HdtEulerTourTree, ComponentSize) {
  const uint32_t n = 5;
  EulerTourTree ett = EulerTourTree{n};

  // 0 - 1
  // | \
  // 2   3 - 4
  parlay::sequence<std::pair<uint32_t, uint32_t>> links{{{0, 1}, {0, 2}, {0, 3}, {4, 3}}};
  ett.BatchLink(links, false);
  EXPECT_EQ(ett.ComponentSize(3), 5);
  EXPECT_THAT(ett.ComponentVertices(3), UnorderedElementsAre(0, 1, 2, 3, 4));
  EXPECT_THAT(
      ett.ComponentEdges(3),
      UnorderedElementsAre(Pair(0, 1), Pair(0, 2), Pair(0, 3), Pair(3, 4)));

  ett.Cut(0, 3);
  EXPECT_EQ(ett.ComponentSize(2), 3);
  EXPECT_THAT(ett.ComponentVertices(2), UnorderedElementsAre(0, 1, 2));
  EXPECT_THAT(ett.ComponentEdges(2), UnorderedElementsAre(Pair(0, 1), Pair(0, 2)));
  EXPECT_EQ(ett.ComponentSize(3), 2);
  EXPECT_THAT(ett.ComponentVertices(3), UnorderedElementsAre(3, 4));
  EXPECT_THAT(ett.ComponentEdges(3), UnorderedElementsAre(Pair(3, 4)));
}

TEST(HdtEulerTourTree, GetAndClearLevelIEdges_Singleton) {
  // Check GetAndClearLevelIEdges() works on a singleton vertex.
  EulerTourTree ett = EulerTourTree{1};
  EXPECT_THAT(ett.GetAndClearLevelIEdges(0), IsEmpty());
}

TEST(HdtEulerTourTree, GetAndClearLevelIEdges) {
  const uint32_t n = 7;
  EulerTourTree ett = EulerTourTree{n};

  // 0 - 1   4   6
  // | \     |
  // 2   3   5
  parlay::sequence<std::pair<uint32_t, uint32_t>> links{{{0, 1}, {0, 2}, {0, 3}, {4, 5}}};
  parlay::sequence<bool> is_level_i_edge{{true, false, true, true}};
  ett.BatchLink(links, is_level_i_edge);

  EXPECT_THAT(ett.GetAndClearLevelIEdges(3), UnorderedElementsAre(Pair(0, 1), Pair(0, 3)));
  EXPECT_THAT(ett.GetAndClearLevelIEdges(3), IsEmpty());
  EXPECT_THAT(ett.GetAndClearLevelIEdges(4), UnorderedElementsAre(Pair(4, 5)));
  EXPECT_THAT(ett.GetAndClearLevelIEdges(4), IsEmpty());
  EXPECT_THAT(ett.GetAndClearLevelIEdges(6), IsEmpty());
}

TEST(HdtEulerTourTree, FindNonTreeEdges_SingletonWithoutEdges) {
  // Check FindNonTreeEdges() works on a singleton vertex without any incident
  // non-tree edges.
  EulerTourTree ett = EulerTourTree{1};
  const auto finder = ett.CreateNontreeEdgeFinder(0);

  MockForEachIncidentVertexCallback mock_callback;
  EXPECT_CALL(mock_callback, Call(_, _, _)).Times(0);

  EXPECT_EQ(finder.NumEdges(), 0);
  finder.ForEachIncidentVertex(0, 100, mock_callback.AsStdFunction());
  finder.ForEachIncidentVertex(20, 40, mock_callback.AsStdFunction());
}

TEST(HdtEulerTourTree, FindNonTreeEdges_SingletonWithEdges) {
  // Check FindNonTreeEdges() works on a singleton vertex with incident non-tree
  // edges.
  EulerTourTree ett = EulerTourTree{1};
  ett.IncrementNontreeEdgeCounts(parlay::sequence<uint32_t>(1, 0), parlay::sequence<uint32_t>(1, 50));
  const auto finder = ett.CreateNontreeEdgeFinder(0);

  InSequence s;
  MockForEachIncidentVertexCallback mock_callback;
  EXPECT_CALL(mock_callback, Call(0, 0, 50));
  EXPECT_CALL(mock_callback, Call(0, 20, 40));
  EXPECT_CALL(mock_callback, Call(0, 1, 2));

  EXPECT_EQ(finder.NumEdges(), 50);
  finder.ForEachIncidentVertex(0, 100, mock_callback.AsStdFunction());
  finder.ForEachIncidentVertex(20, 40, mock_callback.AsStdFunction());
  finder.ForEachIncidentVertex(1, 1, mock_callback.AsStdFunction());
  finder.ForEachIncidentVertex(1, 2, mock_callback.AsStdFunction());
}

TEST(HdtEulerTourTree, FindNonTreeEdges_NoEdges) {
  EulerTourTree ett = EulerTourTree{5};

  // 0 - 1
  // | \
  // 2   3 - 4
  parlay::sequence<std::pair<uint32_t, uint32_t>> links{{{0, 1}, {0, 2}, {0, 3}, {4, 3}}};
  ett.BatchLink(links, false);

  const auto finder = ett.CreateNontreeEdgeFinder(0);
  MockForEachIncidentVertexCallback mock_callback;
  EXPECT_CALL(mock_callback, Call(_, _, _)).Times(0);

  finder.ForEachIncidentVertex(0, 100, mock_callback.AsStdFunction());
  finder.ForEachIncidentVertex(20, 40, mock_callback.AsStdFunction());
}

TEST(HdtEulerTourTree, FindNonTreeEdges) {
  // Check FindNonTreeEdges() works.
  uint32_t n = 5;
  EulerTourTree ett = EulerTourTree{n};
  // 0 - 1
  // | \
  // 2   3   4
  parlay::sequence<std::pair<uint32_t, uint32_t>> links{{{0, 1}, {0, 2}, {0, 3}}};
  ett.BatchLink(links, false);

  const parlay::sequence<uint32_t> vertices = {0, 1, 2, 3, 4};
  const parlay::sequence<uint32_t> edge_counts = {0, 3, 1, 4, 5};
  ett.IncrementNontreeEdgeCounts(vertices, edge_counts);

  // We don't know how the NontreeEdgeFinder will order the non-tree edges, so
  // we'll explicitly track how many times each of the edges has been visited.
  std::vector<std::vector<uint32_t>> edge_visit_counts(5, std::vector<uint32_t>{});

  const auto visit = [&](uint32_t vertex_id, uint32_t start, uint32_t end) {
    for (size_t i = start; i < end; i++) {
      edge_visit_counts[vertex_id][i]++;
    }
  };
  const auto reset_counts = [&]() {
    for (uint32_t i = 0; i < n; i++) {
      edge_visit_counts[i] = std::vector<uint32_t>(edge_counts[i], 0);
    }
  };
  const auto sum_counts = [&]() {
    uint32_t sum = 0;
    for (uint32_t i = 0; i < n; i++) {
      for (uint32_t j = 0; j < edge_counts[i]; j++) {
        sum += edge_visit_counts[i][j];
      }
    }
    return sum;
  };
  const auto unique_visits = [&]() {
    uint32_t count = 0;
    for (uint32_t i = 0; i < n; i++) {
      for (uint32_t j = 0; j < edge_counts[i]; j++) {
        if (edge_visit_counts[i][j] > 0) {
          count++;
        }
      }
    }
    return count;
  };

  const auto finder = ett.CreateNontreeEdgeFinder(0);
  EXPECT_EQ(finder.NumEdges(), 8);

  // Check that querying disjoint intervals gives unique edges and gives the
  // correct number of edges.
  reset_counts();
  finder.ForEachIncidentVertex(0, 2, visit);
  EXPECT_EQ(sum_counts(), 2);
  finder.ForEachIncidentVertex(2, 6, visit);
  EXPECT_EQ(sum_counts(), 6);
  finder.ForEachIncidentVertex(6, 7, visit);
  EXPECT_EQ(sum_counts(), 7);
  EXPECT_EQ(unique_visits(), 7);

  // Check that querying the same interval gives a deterministic answer
  reset_counts();
  finder.ForEachIncidentVertex(1, 4, visit);
  EXPECT_EQ(sum_counts(), 3);
  finder.ForEachIncidentVertex(1, 4, visit);
  EXPECT_EQ(sum_counts(), 6);
  EXPECT_EQ(unique_visits(), 3);
}

}  // namespace hdt_euler_tour_tree
}  // namespace elektra::testing

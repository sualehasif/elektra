#pragma once

#include <parlay/parallel.h>
#include <parlay/sequence.h>

#include <utility>

#include "../../elektra/parallel_euler_tour_tree/euler_tour_tree.hpp"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace elektra::testing {
namespace {

namespace pett = parallel_euler_tour_tree;

using ::testing::IsEmpty;
using ::testing::Pair;
using ::testing::UnorderedElementsAre;

TEST(ParallelEulerTourTreeTest, IsDisconnectedInitially) {
  const pett::EulerTourTree ett = pett::EulerTourTree{3};
  EXPECT_FALSE(ett.IsConnected(0, 1));
  EXPECT_FALSE(ett.IsConnected(0, 2));
  EXPECT_FALSE(ett.IsConnected(1, 2));
}

// Create a small line graph and test that it's connected. Then cut it in half
// and test that it's disconnected.
TEST(ParallelEulerTourTreeTest, ShortLineGraph) {
  const int n = 7;
  pett::EulerTourTree ett = pett::EulerTourTree{n};

  for (int i = 0; i < n - 1; i++) {
    ett.Link(i, i + 1);
  }
  EXPECT_TRUE(ett.IsConnected(0, n - 1));

  ett.Cut(n / 2, n / 2 + 1);
  EXPECT_TRUE(ett.IsConnected(0, n / 2));
  EXPECT_TRUE(ett.IsConnected(n / 2 + 1, n - 1));
  EXPECT_FALSE(ett.IsConnected(0, n - 1));
}

// Adds edges creating two large star graphs linked by an edge, then removes the
// star graph edges.
//
// This larger test is to check that parallel link and parallel cut works as
// expected.
// TODO(tomtseng): The ETT operates on batches sequentially if the batch is
// smaller than some size S. But we may want to tune S, and we don't want that
// to affect this test. We should let this test manually override S by exposing
// the S as a optionally tunable parameter.
TEST(ParallelEulerTourTreeTest, BigStarGraphs) {
  const int n = 200;
  ASSERT_EQ(n % 2, 0);
  pett::EulerTourTree ett = pett::EulerTourTree{n};

  auto links = parlay::sequence<std::pair<int, int>>::uninitialized(n - 1);
  auto cuts = parlay::sequence<std::pair<int, int>>::uninitialized(n - 2);
  // The centers of the two stars are vertices n / 2 - 1 and n / 2.
  for (int i = 0; i < n / 2 - 1; i++) {
    links[i] = make_pair(i, n / 2 - 1);
    cuts[i] = make_pair(i, n / 2 - 1);
  }
  for (int i = n / 2 + 1; i < n; i++) {
    links[i - 2] = make_pair(n / 2, i);
    cuts[i - 2] = make_pair(n / 2, i);
  }
  links[n - 2] = make_pair(n / 2 - 1, n / 2);

  ett.BatchLink(links);
  EXPECT_TRUE(ett.IsConnected(0, 1));
  EXPECT_TRUE(ett.IsConnected(0, n / 2 - 1));
  EXPECT_TRUE(ett.IsConnected(0, n - 1));
  EXPECT_TRUE(ett.IsConnected(n / 2 - 1, n / 2));

  ett.BatchCut(cuts);
  EXPECT_FALSE(ett.IsConnected(0, 1));
  EXPECT_FALSE(ett.IsConnected(0, n / 2 - 1));
  EXPECT_FALSE(ett.IsConnected(0, n - 1));
  EXPECT_TRUE(ett.IsConnected(n / 2 - 1, n / 2));
}

TEST(ParallelEulerTourTreeTest, ComponentVerticesAndEdges) {
  pett::EulerTourTree ett = pett::EulerTourTree{6};
  ett.Link(0, 2);
  ett.Link(1, 2);
  ett.Link(3, 5);

  EXPECT_THAT(ett.ComponentVertices(0), UnorderedElementsAre(0, 1, 2));
  EXPECT_THAT(ett.ComponentEdges(0), UnorderedElementsAre(Pair(0, 2), Pair(1, 2)));
  EXPECT_THAT(ett.ComponentVertices(3), UnorderedElementsAre(3, 5));
  EXPECT_THAT(ett.ComponentEdges(3), UnorderedElementsAre(Pair(3, 5)));
  EXPECT_THAT(ett.ComponentVertices(4), UnorderedElementsAre(4));
  EXPECT_THAT(ett.ComponentEdges(4), IsEmpty());
}

TEST(ParallelEulerTourTreeTest, GetRepresentative) {
  const int n = 8;
  pett::EulerTourTree ett = pett::EulerTourTree{n};

  ett.Link(0, 1);
  ett.Link(0, 2);
  ett.Link(1, 3);
  ett.Link(1, 4);
  ett.Link(6, 7);

  std::vector<int> reps(n);
  for (int i = 0; i < n; i++) {
    reps[i] = ett.GetRepresentative(i);
  }
  EXPECT_EQ(reps[0], reps[1]);
  EXPECT_EQ(reps[0], reps[2]);
  EXPECT_EQ(reps[0], reps[3]);
  EXPECT_EQ(reps[0], reps[4]);
  EXPECT_NE(reps[0], reps[5]);
  EXPECT_NE(reps[0], reps[6]);
  EXPECT_NE(reps[5], reps[6]);
  EXPECT_EQ(reps[6], reps[7]);
}

}  // namespace
}  // namespace elektra::testing

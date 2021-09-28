#pragma once

#include "../../elektra/parallel_euler_tour_tree/euler_tour_tree.hpp"

#include <utility>

#include "gtest/gtest.h"
#include <parlay/parallel.h>
#include <parlay/sequence.h>

namespace elektra::testing {
namespace {

namespace pett = parallel_euler_tour_tree;

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
// This larger test is to check that parallel link and parallel cut works as expected.
// TODO(tomtseng): The ETT operates on batches sequentially if the batch is smaller
// than some size S. But we may want to tune S, and we don't want that to affect
// this test. We should let this test manually override S by exposing the S as a
// optionally tunable parameter.
TEST(ParallelEulerTourTreeTest, BigStarGraphs) {
  const int n = 200;
  pett::EulerTourTree ett = pett::EulerTourTree{n};

  auto links = parlay::sequence<std::pair<int, int>>::uninitialized(199);
  auto cuts = parlay::sequence<std::pair<int, int>>::uninitialized(198);
  for (int i = 0; i <= 99; i++) {
    links[i] = make_pair(i, 100);
    cuts[i] = make_pair(i, 100);
  }
  for (int i = 100; i <= 197; i++) {
    links[i] = make_pair(101, i + 2);
    cuts[i] = make_pair(101, i + 2);
  }
  links[198] = make_pair(100, 101);

  ett.BatchLink(links, links.size());
  EXPECT_TRUE(ett.IsConnected(0, 1));
  EXPECT_TRUE(ett.IsConnected(0, 100));
  EXPECT_TRUE(ett.IsConnected(0, 199));
  EXPECT_TRUE(ett.IsConnected(100, 101));

  ett.BatchCut(cuts, cuts.size());
  EXPECT_FALSE(ett.IsConnected(0, 1));
  EXPECT_FALSE(ett.IsConnected(0, 100));
  EXPECT_FALSE(ett.IsConnected(0, 199));
  EXPECT_TRUE(ett.IsConnected(100, 101));
}

}  // namespace
}  // namespace elektra::testing

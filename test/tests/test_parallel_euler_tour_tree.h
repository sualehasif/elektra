#include <parlay/parallel.h>
#include <parlay/random.h>
#include <parlay/sequence.h>

#include <random>
#include <string>
#include <utility>

#include "../../elektra/parallel_euler_tour_tree/euler_tour_tree.h"
#include "../../elektra/utilities/simple_forest_connectivity.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

// TODO: parallelize the loops in this based on the original.

namespace elektra::testing {
namespace {
using EulerTourTree = parallel_euler_tour_tree::EulerTourTree;

using ::testing::ElementsAreArray;
using ::testing::IsEmpty;
using ::testing::Pair;
using ::testing::UnorderedElementsAre;

TEST(ParallelEulerTourTreeTest, RandomLinksAndCuts) {
  // Generate `link_attempts_per_round` edges randomly, keeping each one that
  // doesn't add a cycle into the forest. Then call `BatchLink` on all of
  // them.

  constexpr uint32_t num_vertices{100};
  constexpr int num_rounds{2};
  constexpr int link_attempts_per_round{200};
  constexpr int cut_ratio{3};

  EulerTourTree ett{num_vertices};
  SimpleForestConnectivity reference_solution{num_vertices};
  std::unordered_set<std::pair<uint32_t, uint32_t>, HashIntPairStruct> edges{};
  std::mt19937 rng{0};
  std::uniform_int_distribution<std::mt19937::result_type> vert_dist{0, num_vertices - 1};
  std::uniform_int_distribution<std::mt19937::result_type> coin{0, 1};

  for (int i = 0; i < num_rounds; i++) {
    parlay::sequence<std::pair<uint32_t, uint32_t>> links{};
    for (int j = 0; j < link_attempts_per_round; j++) {
      const unsigned long u{vert_dist(rng)}, v{vert_dist(rng)};
      if (!reference_solution.IsConnected(u, v)) {
        reference_solution.Link(u, v);
        edges.emplace(u, v);
        links.push_back(std::make_pair(u, v));
      }
    }
    ett.BatchLink(links);

    for (unsigned u = 0; u < num_vertices; u++) {
      for (unsigned v = 0; v < num_vertices; v++) {
        EXPECT_EQ(reference_solution.IsConnected(u, v), ett.IsConnected(u, v));
      }
    }

    // Call `BatchCut` over each `cut_ratio`-th edge.
    parlay::sequence<std::pair<unsigned, unsigned>> cuts{};
    int cnt{0};
    for (auto e : edges) {
      if (++cnt % cut_ratio == 0) {
        cuts.push_back(e);
      }
    }
    for (size_t j = 0; j < cuts.size(); j++) {
      pair<unsigned, unsigned> cut{cuts[j]};
      edges.erase(cut);
      if (coin(rng) == 1) {
        cuts[j] = make_pair(cut.second, cut.first);
      }
      reference_solution.Cut(cut.first, cut.second);
    }
    ett.BatchCut(cuts);

    for (unsigned u = 0; u < num_vertices; u++) {
      for (unsigned v = 0; v < num_vertices; v++) {
        EXPECT_EQ(reference_solution.IsConnected(u, v), ett.IsConnected(u, v));
      }
    }
  }
}

TEST(ParallelEulerTourTreeTest, IsDisconnectedInitially) {
  const EulerTourTree ett = EulerTourTree{3};
  EXPECT_FALSE(ett.IsConnected(0, 1));
  EXPECT_FALSE(ett.IsConnected(0, 2));
  EXPECT_FALSE(ett.IsConnected(1, 2));
}

// Create a small line graph and test that it's connected. Then cut it in half
// and test that it's disconnected.
TEST(ParallelEulerTourTreeTest, ShortLineGraph) {
  const unsigned n = 7;
  EulerTourTree ett = EulerTourTree{n};

  for (unsigned i = 0; i < n - 1; i++) {
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
  const unsigned n = 200;
  ASSERT_EQ(n % 2, 0);
  EulerTourTree ett = EulerTourTree{n};

  auto links = parlay::sequence<std::pair<unsigned, unsigned>>::uninitialized(n - 1);
  auto cuts = parlay::sequence<std::pair<unsigned, unsigned>>::uninitialized(n - 2);
  // The centers of the two stars are vertices n / 2 - 1 and n / 2.
  for (unsigned i = 0; i < n / 2 - 1; i++) {
    links[i] = make_pair(i, n / 2 - 1);
    cuts[i] = make_pair(i, n / 2 - 1);
  }
  for (unsigned i = n / 2 + 1; i < n; i++) {
    links[i - 2] = make_pair(n / 2, i);
    cuts[i - 2] = make_pair(n / 2, i);
  }
  links[n - 2] = make_pair(n / 2 - 1, n / 2);

  ett.BatchLink(links);
  EXPECT_TRUE(ett.IsConnected(0, 1));
  EXPECT_TRUE(ett.IsConnected(0, n / 2 - 1));
  EXPECT_TRUE(ett.IsConnected(0, n - 1));
  EXPECT_TRUE(ett.IsConnected(n / 2 - 1, n / 2));
  EXPECT_EQ(ett.ComponentSize(5), 200);

  ett.BatchCut(cuts);
  EXPECT_FALSE(ett.IsConnected(0, 1));
  EXPECT_FALSE(ett.IsConnected(0, n / 2 - 1));
  EXPECT_FALSE(ett.IsConnected(0, n - 1));
  EXPECT_TRUE(ett.IsConnected(n / 2 - 1, n / 2));
  EXPECT_EQ(ett.ComponentSize(5), 1);
  EXPECT_EQ(ett.ComponentSize(n / 2 - 1), 2);
}

TEST(ParallelEulerTourTreeTest, ComponentVerticesAndEdges) {
  EulerTourTree ett = EulerTourTree{6};
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

TEST(ParallelEulerTourTreeTest, ComponentEdgesOfLargeTree) {
  // Make sure ComponentEdges still works when the ETT is large (and functions
  // presumably run in parallel
  const unsigned n = 400;
  EulerTourTree ett = EulerTourTree{n};
  std::vector<std::pair<unsigned, unsigned>> expected_edges;
  expected_edges.reserve(n - 1);
  for (unsigned i = 1; i < n / 3; i++) {
    expected_edges.emplace_back(0, i);
    ett.Link(0, i);
  }
  for (unsigned i = n / 3; i < n; i++) {
    expected_edges.emplace_back(i - 1, i);
    ett.Link(i - 1, i);
  }

  auto edges = ett.ComponentEdges(0);
  std::sort(edges.begin(), edges.end());
  EXPECT_THAT(edges, ElementsAreArray(expected_edges));
}

TEST(ParallelEulerTourTreeTest, GetRepresentative) {
  const unsigned n = 8;
  EulerTourTree ett = EulerTourTree{n};

  ett.Link(0, 1);
  ett.Link(0, 2);
  ett.Link(1, 3);
  ett.Link(1, 4);
  ett.Link(6, 7);

  std::vector<unsigned> reps(n);
  for (unsigned i = 0; i < n; i++) {
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

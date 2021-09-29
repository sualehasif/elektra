#include <parlay/parallel.h>
#include <parlay/random.h>
#include <parlay/sequence.h>

#include <random>
#include <string>

#include "../../elektra/parallel_euler_tour_tree/euler_tour_tree.hpp"
#include "../../elektra/utilities/simple_forest_connectivity.h"
#include "gtest/gtest.h"

// TODO: parallelize the loops in this based on the original.

constexpr int num_vertices{100};
constexpr int link_attempts_per_round{200};
constexpr int cut_ratio{3};
constexpr int num_rounds{2};

namespace elektra::testing {
namespace {
using EulerTourTree = parallel_euler_tour_tree::EulerTourTree;

class ParallelETTTest : public ::testing::Test {
 protected:
  ParallelETTTest() {
    std::cout << "Parallel ETT: Testing";
    rng.seed(0);
  }

  ~ParallelETTTest() override {}

  std::mt19937 rng{};

  std::uniform_int_distribution<std::mt19937::result_type> vert_dist{
      0, num_vertices - 1};
  std::uniform_int_distribution<std::mt19937::result_type> coin{0, 1};

  SimpleForestConnectivity reference_solution{num_vertices};
  EulerTourTree ett{num_vertices};
  std::unordered_set<std::pair<int, int>, HashIntPairStruct> edges{};

  parlay::sequence<std::pair<int, int>> ett_input =
      parlay::sequence<std::pair<int, int>>{};

  //   void CheckAllPairsConnectivity(
  //       const SimpleForestConnectivity& reference_solution,
  //       const EulerTourTree& ett) {
  //     for (int u = 0; u < num_vertices; u++) {
  //       for (int v = 0; v < num_vertices; v++) {
  //         if (reference_solution.IsConnected(u, v) != ett.IsConnected(u, v))
  //         {
  //           std::cout << "Error: " << u << " " << v << std::endl;
  //         }

  //         assert(reference_solution.IsConnected(u, v) == ett.IsConnected(u,
  //         v));
  //       }
  //     }
  //   }
};

TEST_F(ParallelETTTest, SingleRoundLinkingElements) {
  // Generate `link_attempts_per_round` edges randomly, keeping each one that
  // doesn't add a cycle into the forest. Then call `BatchLink` on all of
  // them.

  for (int j = 0; j < link_attempts_per_round; j++) {
    const unsigned long u{vert_dist(rng)}, v{vert_dist(rng)};
    if (!reference_solution.IsConnected(u, v)) {
      reference_solution.Link(u, v);
      edges.emplace(u, v);
      ett_input.push_back(std::make_pair(u, v));
    }
  }
  ett.BatchLink(ett_input);

  std::cout << "Edges: ";
  for (auto& edge : edges) {
    std::cout << edge.first << "<->" << edge.second << "  ||   ";
  }
  std::cout << std::endl;

  for (int u = 0; u < num_vertices; u++) {
    for (int v = 0; v < num_vertices; v++) {
      if (reference_solution.IsConnected(u, v) != ett.IsConnected(u, v)) {
        std::cout << "Error: Connected component of " << u << ": ";
        for (auto& c : ett.ConnectedComponent(u)) {
          std::cout << c << " ";
        }
        std::cout << std::endl;
        std::cout << "Error: Connected component of " << v << ": ";
        for (auto& c : ett.ConnectedComponent(v)) {
          std::cout << c << " ";
        }
        std::cout << std::endl;
      }

      EXPECT_EQ(reference_solution.IsConnected(u, v), ett.IsConnected(u, v));
    }
  }
}

TEST_F(ParallelETTTest, MultipleRoundLinkingElements) {
  // Empty ett_input
  ett_input.assign(0, std::make_pair(0, 0));

  EXPECT_EQ(ett_input.size(), 0);

  for (int i = 0; i < num_rounds; i++) {
    // Generate `link_attempts_per_round` edges randomly, keeping each one that
    // doesn't add a cycle into the forest. Then call `BatchLink` on all of
    // them.

    for (int j = 0; j < link_attempts_per_round; j++) {
      const unsigned long u{vert_dist(rng)}, v{vert_dist(rng)};
      if (!reference_solution.IsConnected(u, v)) {
        reference_solution.Link(u, v);
        edges.emplace(u, v);
        ett_input.push_back(std::make_pair(u, v));
      }
    }
    ett.BatchLink(ett_input);

    for (int u = 0; u < num_vertices; u++) {
      for (int v = 0; v < num_vertices; v++) {
        if (reference_solution.IsConnected(u, v) != ett.IsConnected(u, v)) {
          std::cout << "Error: " << u << " " << v << std::endl;
        }

        EXPECT_EQ(reference_solution.IsConnected(u, v), ett.IsConnected(u, v));
      }
    }

    std::cout << "Round " << i << ": " << ett_input.size()
              << " edges inserted\n";
  }
}

}  // namespace
}  // namespace elektra::testing

#include <iostream>
#include <random>

#include "elektra/batch_dynamic_connectivity/connectivity.h"
#include "elektra/utilities/simple_forest_connectivity.h"

using EulerTourTree = parallel_euler_tour_tree::EulerTourTree;

constexpr int num_vertices{500};
constexpr int link_attempts_per_round{400};
constexpr int cut_ratio{3};
constexpr int num_rounds{10};

void CheckAllPairsConnectivity(
    const SimpleForestConnectivity& reference_solution,
    const EulerTourTree& ett) {
  for (int u = 0; u < num_vertices; u++) {
    for (int v = 0; v < num_vertices; v++) {
      if (reference_solution.IsConnected(u, v) != ett.IsConnected(u, v)) {
        std::cout << "Error: " << u << " " << v << std::endl;
      }

      assert(reference_solution.IsConnected(u, v) == ett.IsConnected(u, v));
    }
  }
}

int main() {
  std::cout << "Hello, from elektra!\n";

  std::cout << "Parallel ETT: Test 1";

  std::mt19937 rng{};
  rng.seed(0);
  std::uniform_int_distribution<std::mt19937::result_type> vert_dist{
      0, num_vertices - 1};
  std::uniform_int_distribution<std::mt19937::result_type> coin{0, 1};

  SimpleForestConnectivity reference_solution{num_vertices};
  EulerTourTree ett{num_vertices};
  std::unordered_set<std::pair<int, int>, HashIntPairStruct> edges{};

  auto ett_input = parlay::sequence<std::pair<int, int>>{};

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
    ett.BatchLink(ett_input, ett_input.size());
    CheckAllPairsConnectivity(reference_solution, ett);

    std::cout << "Round " << i << ": " << ett_input.size()
              << " edges inserted\n";
    std::cout << "Passed\n";

    // auto ett_input_cut = parlay::sequence<std::pair<int, int>>{};
    // // Call `BatchCut` over each `cut_ratio`-th edge.
    // int cnt{0};
    // for (auto e : edges) {
    //   if (++cnt % cut_ratio == 0) {
    //     ett_input_cut.push_back(e);
    //   }
    // }
    // for (int j = 0; j < ett_input_cut.size(); j++) {
    //   pair<int, int> cut{ett_input_cut[j]};
    //   edges.erase(cut);
    //   if (coin(rng) == 1) {
    //     ett_input_cut[j] = make_pair(cut.second, cut.first);
    //   }
    //   reference_solution.Cut(cut.first, cut.second);
    // }
    // ett.BatchCut(ett_input_cut, ett_input_cut.size());
    // CheckAllPairsConnectivity(reference_solution, ett);
  }

  std::cout << "Test 1 complete." << std::endl;
}

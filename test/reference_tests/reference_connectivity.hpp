#include <parlay/parallel.h>
#include <parlay/random.h>
#include <parlay/sequence.h>

#include <dynamic_graph/dynamic_connectivity.hpp>
#include <random>
#include <string>

#include "../../elektra/batch_dynamic_connectivity/dynamic_connectivity.h"
#include "gtest/gtest.h"

namespace elektra::testing {
namespace {
using std::make_pair;

using UndirectedEdge = pair<V, V>;
using BatchDynamicConnectivity = bdcty::BatchDynamicConnectivity;

class ParallelConnectivityTest : public ::testing::Test {
 protected:
  ParallelConnectivityTest() {
    std::cout << "Testing on " << parlay::num_workers() << " Parlay workers."
              << std::endl;
    rng.seed(0);
  }

  std::mt19937 rng{};

  //   std::uniform_int_distribution<std::mt19937::result_type> vert_dist{
  //       0, num_vertices - 1};
  //   std::uniform_int_distribution<std::mt19937::result_type> coin{0, 1};

  //   SimpleForestConnectivity reference_solution{num_vertices};
  //   EulerTourTree ett{num_vertices};
  //   std::unordered_set<std::pair<int, int>, HashIntPairStruct> edges{};

  //   parlay::sequence<std::pair<int, int>> ett_input =
  //       parlay::sequence<std::pair<int, int>>{};
};

TEST(ParallelConnectivityTest, SimpleGraph) {
  DynamicConnectivity graph(1);
  EXPECT_TRUE(graph.IsConnected(0, 0));
}

}  // namespace
}  // namespace elektra::testing
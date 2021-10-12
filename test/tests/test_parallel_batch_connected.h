#include <parlay/parallel.h>
#include <parlay/random.h>
#include <parlay/sequence.h>

#include <random>
#include <string>

#include "../../elektra/batch_dynamic_connectivity/dynamic_connectivity.h"
#include "gtest/gtest.h"

namespace elektra::testing {
namespace {

using UndirectedEdge = dynamicGraph::UndirectedEdge;
using BatchDynamicConnectivity =
    batchDynamicConnectivity::BatchDynamicConnectivity;

class ParallelConnectivityTest : public ::testing::Test {
 protected:
  ParallelConnectivityTest() {
    std::cout << "Parallel Connectivity: Testing" << std::endl;
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

TEST_F(ParallelConnectivityTest, SimpleEdgeInsertion1) {
  parlay::sequence<UndirectedEdge> edges;

  edges.push_back(UndirectedEdge(0, 1));
  edges.push_back(UndirectedEdge(1, 2));
  edges.push_back(UndirectedEdge(3, 4));

  BatchDynamicConnectivity x(5, edges);

  parlay::sequence<std::pair<Vertex, Vertex>> queries;
  parlay::sequence<char> expectedOut;

  queries.push_back(std::make_pair(0, 1));

  expectedOut.push_back(true);

  auto result = x.BatchConnected(queries);
  for (int i = 0; i < (int)queries.size(); i++) {
    EXPECT_EQ(result[i], expectedOut[i]);
  }
}

TEST_F(ParallelConnectivityTest, SimpleInsertionAndQuery1) {
  parlay::sequence<UndirectedEdge> edges;

  edges.push_back(UndirectedEdge(0, 1));
  edges.push_back(UndirectedEdge(1, 2));
  edges.push_back(UndirectedEdge(3, 4));

  BatchDynamicConnectivity x(5, edges);

  parlay::sequence<std::pair<Vertex, Vertex>> queries;
  parlay::sequence<char> expectedOut;

  queries.push_back(std::make_pair(0, 1));
  queries.push_back(std::make_pair(1, 2));
  queries.push_back(std::make_pair(0, 2));
  queries.push_back(std::make_pair(3, 4));
  queries.push_back(std::make_pair(0, 3));
  queries.push_back(std::make_pair(0, 4));
  queries.push_back(std::make_pair(1, 3));
  queries.push_back(std::make_pair(1, 4));

  expectedOut.push_back(true);
  expectedOut.push_back(true);
  expectedOut.push_back(true);
  expectedOut.push_back(true);
  expectedOut.push_back(false);
  expectedOut.push_back(false);
  expectedOut.push_back(false);
  expectedOut.push_back(false);

  auto result = x.BatchConnected(queries);
  for (int i = 0; i < (int)queries.size(); i++) {
    EXPECT_EQ(result[i], expectedOut[i]);
  }
}

TEST_F(ParallelConnectivityTest, SimpleInsertionAndQuery2) {
  parlay::sequence<UndirectedEdge> edges;
  // connected component with first four edges
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < i; j++) {
      edges.push_back(UndirectedEdge(i, j));
    }
  }
  edges.push_back(UndirectedEdge(5, 6));
  edges.push_back(UndirectedEdge(6, 7));

  BatchDynamicConnectivity x(10, edges);

  parlay::sequence<std::pair<Vertex, Vertex>> queries;
  parlay::sequence<char> expectedOut;

  for (long i = 0; i < 10; i++) {
    for (long j = 0; j < 10; j++) {
      queries.push_back(std::make_pair(i, j));
      expectedOut.push_back((i == j) || ((i < 5) && (j < 5)) ||
                            ((5 <= i) && (5 <= j) && (i < 8) &&
                             (j < 8)));  //|| ((i < 10) && (j < 10))
    }
  }

  auto result = x.BatchConnected(queries);
  for (int i = 0; i < (int)queries.size(); i++) {
    EXPECT_EQ(result[i], expectedOut[i]);
  }
}

}  // namespace
}  // namespace elektra::testing

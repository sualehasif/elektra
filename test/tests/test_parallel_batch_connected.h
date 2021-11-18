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

TEST_F(ParallelConnectivityTest, SimpleEdgeDeletion1) {
  parlay::sequence<UndirectedEdge> edges;

  edges.push_back(UndirectedEdge(0, 1));
  edges.push_back(UndirectedEdge(1, 2));
  edges.push_back(UndirectedEdge(3, 4));

  BatchDynamicConnectivity x(5, edges);

  parlay::sequence<std::pair<Vertex, Vertex>> queries;
  parlay::sequence<char> expectedOut;

  queries.push_back(std::make_pair(0, 2));

  expectedOut.push_back(true);

  auto result = x.BatchConnected(queries);
  for (int i = 0; i < (int)queries.size(); i++) {
    EXPECT_EQ(result[i], expectedOut[i]);
  }

  // x.PrintStructure();

  // do some deletions
  parlay::sequence<UndirectedEdge> deletions;
  deletions.push_back(UndirectedEdge{1, 2});

  // delete the edge in the graph
  x.BatchDeleteEdges(deletions);

  // x.PrintStructure();

  // do queries again
  queries.clear();
  expectedOut.clear();

  queries.push_back(std::make_pair(0, 2));
  expectedOut.push_back(false);

  result = x.BatchConnected(queries);
  for (int i = 0; i < (int)queries.size(); i++) {
    EXPECT_EQ(result[i], expectedOut[i]);
  }
}

TEST_F(ParallelConnectivityTest, SimpleEdgeDeletion2) {
  parlay::sequence<UndirectedEdge> edges;

  // we are looking at a graph with a triangle with an edge sticking out of it
  edges.push_back(UndirectedEdge(0, 1));
  edges.push_back(UndirectedEdge(1, 2));
  edges.push_back(UndirectedEdge(0, 2));
  edges.push_back(UndirectedEdge(2, 3));

  BatchDynamicConnectivity x(4, edges);

  // x.PrintStructure();

  parlay::sequence<std::pair<Vertex, Vertex>> queries;
  parlay::sequence<char> expectedOut;

  // set up the sequence of queries
  queries.push_back(std::make_pair(0, 2));
  expectedOut.push_back(true);

  queries.push_back(std::make_pair(0, 1));
  expectedOut.push_back(true);

  queries.push_back(std::make_pair(0, 3));
  expectedOut.push_back(true);

  auto result = x.BatchConnected(queries);
  for (int i = 0; i < (int)queries.size(); i++) {
    EXPECT_EQ(result[i], expectedOut[i]);
  }

  // do some deletions
  parlay::sequence<UndirectedEdge> deletions;
  deletions.push_back(UndirectedEdge{1, 2});

  // delete the edge in the graph
  x.BatchDeleteEdges(deletions);

  // do queries again
  queries.clear();
  expectedOut.clear();

  // set up the sequence of queries
  queries.push_back(std::make_pair(0, 2));
  expectedOut.push_back(true);

  queries.push_back(std::make_pair(0, 1));
  expectedOut.push_back(true);

  queries.push_back(std::make_pair(0, 3));
  expectedOut.push_back(true);

  result = x.BatchConnected(queries);
  for (int i = 0; i < (int)queries.size(); i++) {
    EXPECT_EQ(result[i], expectedOut[i]);
  }
}

TEST_F(ParallelConnectivityTest, SimpleEdgeDeletion3) {
  parlay::sequence<UndirectedEdge> edges;

  auto num_vertices = 5;

  // we are looking at a graph:
  //     o 0
  //    / \
  // 2 o---o 1
  //    \ / \
  //   3 o---o 4
  edges.push_back(UndirectedEdge(0, 1));
  edges.push_back(UndirectedEdge(1, 2));
  edges.push_back(UndirectedEdge(0, 2));
  edges.push_back(UndirectedEdge(2, 3));
  edges.push_back(UndirectedEdge(3, 1));
  edges.push_back(UndirectedEdge(1, 4));
  edges.push_back(UndirectedEdge(3, 4));

  BatchDynamicConnectivity x(num_vertices, edges);

  parlay::sequence<std::pair<Vertex, Vertex>> queries;
  parlay::sequence<char> expectedOut;

  // set up the sequence of queries
  queries.push_back(std::make_pair(0, 1));
  expectedOut.push_back(true);

  queries.push_back(std::make_pair(0, 2));
  expectedOut.push_back(true);

  queries.push_back(std::make_pair(0, 4));
  expectedOut.push_back(true);

  auto result = x.BatchConnected(queries);
  for (int i = 0; i < (int)queries.size(); i++) {
    EXPECT_EQ(result[i], expectedOut[i]);
  }

  // do some deletions
  // the graph is now:
  //     o 0
  //    / \
  // 2 o   o 1
  //    \   \
  //   3 o   o 4
  parlay::sequence<UndirectedEdge> deletions;
  deletions.push_back(UndirectedEdge{1, 2});
  deletions.push_back(UndirectedEdge{1, 3});
  deletions.push_back(UndirectedEdge{3, 4});

  // delete the edge in the graph
  x.BatchDeleteEdges(deletions);

  // do queries again
  queries.clear();
  expectedOut.clear();

  // set up the sequence of queries
  queries.push_back(std::make_pair(0, 1));
  expectedOut.push_back(true);

  // set up the sequence of queries
  queries.push_back(std::make_pair(1, 2));
  expectedOut.push_back(true);

  queries.push_back(std::make_pair(0, 2));
  expectedOut.push_back(true);

  queries.push_back(std::make_pair(0, 4));
  expectedOut.push_back(true);

  queries.push_back(std::make_pair(3, 4));
  expectedOut.push_back(true);

  result = x.BatchConnected(queries);
  for (int i = 0; i < (int)queries.size(); i++) {
    EXPECT_EQ(result[i], expectedOut[i]);
  }
}

TEST_F(ParallelConnectivityTest, SimpleEdgeDeletion4) {
  parlay::sequence<UndirectedEdge> edges;

  auto num_vertices = 5;

  // we are looking at a graph:
  //     o 0
  //    / \
  // 2 o---o 1
  //    \ / \
  //   3 o---o 4
  edges.push_back(UndirectedEdge(0, 1));
  edges.push_back(UndirectedEdge(1, 2));
  edges.push_back(UndirectedEdge(0, 2));
  edges.push_back(UndirectedEdge(2, 3));
  edges.push_back(UndirectedEdge(3, 1));
  edges.push_back(UndirectedEdge(1, 4));
  edges.push_back(UndirectedEdge(3, 4));

  BatchDynamicConnectivity x(num_vertices, edges);

  parlay::sequence<std::pair<Vertex, Vertex>> queries;
  parlay::sequence<char> expectedOut;

  // set up the sequence of queries
  queries.push_back(std::make_pair(0, 1));
  expectedOut.push_back(true);

  queries.push_back(std::make_pair(0, 2));
  expectedOut.push_back(true);

  queries.push_back(std::make_pair(0, 4));
  expectedOut.push_back(true);

  auto result = x.BatchConnected(queries);
  for (int i = 0; i < (int)queries.size(); i++) {
    EXPECT_EQ(result[i], expectedOut[i]);
  }

  // Print the graph
  // x.PrintStructure();

  // do some deletions
  //     o 0
  //    /
  // 2 o   o 1
  //    \   \
  //   3 o   o 4
  parlay::sequence<UndirectedEdge> deletions;
  deletions.push_back(UndirectedEdge{0, 1});
  deletions.push_back(UndirectedEdge{1, 2});
  deletions.push_back(UndirectedEdge{1, 3});
  deletions.push_back(UndirectedEdge{3, 4});

  // delete the edge in the graph
  x.BatchDeleteEdges(deletions);

  // Print the graph
  // x.PrintStructure();

  // do queries again
  queries.clear();
  expectedOut.clear();

  // set up the sequence of queries
  queries.push_back(std::make_pair(0, 3));
  expectedOut.push_back(true);

  queries.push_back(std::make_pair(0, 2));
  expectedOut.push_back(true);

  queries.push_back(std::make_pair(2, 3));
  expectedOut.push_back(true);

  queries.push_back(std::make_pair(1, 4));
  expectedOut.push_back(true);

  queries.push_back(std::make_pair(0, 1));
  expectedOut.push_back(false);

  queries.push_back(std::make_pair(1, 2));
  expectedOut.push_back(false);

  queries.push_back(std::make_pair(0, 4));
  expectedOut.push_back(false);

  queries.push_back(std::make_pair(1, 3));
  expectedOut.push_back(false);

  queries.push_back(std::make_pair(3, 4));
  expectedOut.push_back(false);

  result = x.BatchConnected(queries);
  for (int i = 0; i < (int)queries.size(); i++) {
    EXPECT_EQ(result[i], expectedOut[i]);
  }
}

// TEST_F(ParallelConnectivityTest, RandomInsertionAndQuery1) {}
// TEST_F(ParallelConnectivityTest, RandomInsertionAndQuery2) {}

}  // namespace
}  // namespace elektra::testing

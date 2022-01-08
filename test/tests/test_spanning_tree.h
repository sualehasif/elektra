#include <parlay/parallel.h>
#include <parlay/random.h>
#include <parlay/sequence.h>

#include <random>

#include "../../elektra/spanning_tree.h"
#include "gtest/gtest.h"

namespace elektra::testing {
namespace {

class SpanningTree : public ::testing::Test {
 protected:
  using E = std::pair<uintE, uintE>;
  SpanningTree() {
    std::cout << "Testing on " << parlay::num_workers() << " Parlay workers."
              << std::endl;
    rng.seed(0);
  }

  std::mt19937 rng{};
};

// Test that the spanning tree is correct for a small graph.
TEST_F(SpanningTree, MiniGraph) {
  constexpr int num_vertices = 5;
  parlay::sequence<uintE> vertices(num_vertices);

  for (int i = 0; i < num_vertices; ++i) {
    vertices[i] = i;
  }

  parlay::sequence<E> edges = {E(0, 1), E(0, 2), E(1, 2), E(1, 3),
                               E(2, 3), E(2, 4), E(3, 4)};

  parlay::sequence<E> tree;

  elektra::SpanningTree(edges, tree);
  bdcty::PrintEdgeSequence(tree, "Spanning tree");

  bdcty::PrintEdgeSequence(tree, "Spanning tree cut");

  // print the spanning tree
  // std::cout << "Spanning tree: " << std::endl;
  // for (auto& e : tree) {
  //   std::cout << e.first << " " << e.second << std::endl;
  // }

  EXPECT_EQ(tree.size(), num_vertices - 1);

  // initialize a set of all vertices returned by the spanning tree
  std::set<uintE> vertices_in_spanning_tree;
  for (auto &edge : tree) {
    vertices_in_spanning_tree.insert(edge.first);
    vertices_in_spanning_tree.insert(edge.second);
  }

  // check that all vertices are in the spanning tree
  for (int i = 0; i < num_vertices; ++i) {
    EXPECT_EQ(vertices_in_spanning_tree.count(i), 1);
  }
}

// Test that the spanning tree is correct for a small graph.
TEST_F(SpanningTree, MiniGraph2) {
  constexpr int num_vertices = 5;
  parlay::sequence<uintE> vertices(num_vertices);

  for (int i = 0; i < num_vertices; ++i) {
    vertices[i] = i + 5;
  }

  parlay::sequence<E> edges = {E(0, 1), E(0, 2), E(1, 2), E(1, 3),
                               E(2, 3), E(2, 4), E(3, 4)};

  // increase the value of each vertex in edges by 5
  for (auto &e : edges) {
    e.first += 5;
    e.second += 5;
  }

  parlay::sequence<E> tree;

  elektra::SpanningTree(edges, tree);

  // print the spanning tree
  // std::cout << "Spanning tree: " << std::endl;
  // for (auto& e : tree) {
  //   std::cout << e.first << " " << e.second << std::endl;
  // }

  EXPECT_EQ(tree.size(), num_vertices - 1);

  // initialize a set of all vertices returned by the spanning tree
  std::set<uintE> vertices_in_spanning_tree;
  for (auto &edge : tree) {
    vertices_in_spanning_tree.insert(edge.first);
    vertices_in_spanning_tree.insert(edge.second);
  }

  // check that all vertices are in the spanning tree
  for (int i = 0; i < num_vertices; ++i) {
    EXPECT_EQ(vertices_in_spanning_tree.count(vertices[i]), 1);
  }
}

TEST_F(SpanningTree, FixedSmallGraph) {
  constexpr uintE num_vertices = 10;
  constexpr int num_edges = 20;

  parlay::sequence<uintE> vertices(num_vertices);
  parlay::sequence<E> edges(num_edges);

  for (uintE i = 0; i < num_vertices; ++i) {
    vertices[i] = i;
  }

  for (uintE i = 0; i < num_edges; ++i) {
    edges[i] = E(rng() % num_vertices, rng() % num_vertices);
  }

  edges =
      parlay::filter(edges, [&](const E &e) { return e.first != e.second; });

  // push back an edges starting with every vertex so that every vertex is
  // is initalized in the spanning tree
  for (uintE i = 0; i < num_vertices; ++i) {
    auto x = (rng() % num_vertices);
    if (x != i) {
      edges.push_back(E(i, x));
    } else {
      edges.push_back(E(i, i + 1));
    }
  }

  parlay::sequence<E> tree;

  elektra::SpanningTree(edges, tree);

  EXPECT_EQ(tree.size(), num_vertices - 1);

  // initialize a set of all vertices returned by the spanning tree
  std::set<uintE> vertices_in_spanning_tree;
  for (auto &edge : tree) {
    vertices_in_spanning_tree.insert(edge.first);
    vertices_in_spanning_tree.insert(edge.second);
  }

  // check that all vertices are in the spanning tree
  for (uintE i = 0; i < num_vertices; ++i) {
    EXPECT_EQ(vertices_in_spanning_tree.count(i), 1);
  }
}

TEST_F(SpanningTree, BigGraph) {
  constexpr uintE num_vertices = 50;
  constexpr int num_edges = 200;
  rng.seed(162763);

  parlay::sequence<uintE> vertices(num_vertices);
  parlay::sequence<E> edges(num_edges);

  for (uintE i = 0; i < num_vertices; ++i) {
    vertices[i] = i;
  }

  for (uintE i = 0; i < num_edges; ++i) {
    edges[i] = E(rng() % num_vertices, rng() % num_vertices);
  }

  edges =
      parlay::filter(edges, [&](const E &e) { return e.first != e.second; });

  // push back an edges starting with every vertex so that every vertex is
  // is initalized in the spanning tree
  for (uintE i = 0; i < num_vertices; ++i) {
    auto x = (rng() % num_vertices);
    if (x != i) {
      edges.push_back(E(i, x));
    } else {
      edges.push_back(E(i, i + 1));
    }
  }

  parlay::sequence<E> tree;

  elektra::SpanningTree(edges, tree);

  EXPECT_EQ(tree.size(), num_vertices - 1);

  // initialize a set of all vertices returned by the spanning tree
  std::set<uintE> vertices_in_spanning_tree;
  for (auto &edge : tree) {
    vertices_in_spanning_tree.insert(edge.first);
    vertices_in_spanning_tree.insert(edge.second);
  }

  // check that all vertices are in the spanning tree
  for (uintE i = 0; i < num_vertices; ++i) {
    EXPECT_EQ(vertices_in_spanning_tree.count(i), 1);
  }
}

}  // namespace
}  // namespace elektra::testing
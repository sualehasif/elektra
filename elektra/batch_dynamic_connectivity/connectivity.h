#pragma once

#include <parlay/hash_table.h>
#include <parlay/primitives.h>
#include <parlay/sequence.h>

#include <boost/pending/disjoint_sets.hpp>
#include <boost/property_map/property_map.hpp>
#include <cstdint>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "graph.h"
#include "hash_pair.h"
#include "parallel_euler_tour_tree/euler_tour_tree.hpp"
#include "resizable_table.h"

#define INITIAL_SIZE 50

/** This class represents an undirected graph that can undergo efficient edge
 *  insertions, edge deletions, and connectivity queries. Multiple edges between
 *  a pair of vertices are supported.
 */

namespace batchDynamicConnectivity {

using UndirectedEdge = dynamicGraph::UndirectedEdge;
using UndirectedEdgeHash = dynamicGraph::UndirectedEdgeHash;
using BatchDynamicET = parallel_euler_tour_tree::EulerTourTree;

using treeSet = std::unordered_set<UndirectedEdge, UndirectedEdgeHash>;
using vertexSet = parlay::hashtable<parlay::hash_numeric<Vertex>>;

class BatchDynamicConnectivity {
 public:
  /** Initializes an empty graph with a fixed number of vertices.
   *
   *  @param[in] num_vertices Number of vertices in the graph.
   */
  explicit BatchDynamicConnectivity(int num_vertices);

  explicit BatchDynamicConnectivity(int num_vertices,
                                    const parlay::sequence<UndirectedEdge> &se);

  /** Deallocates the data structure. */
  ~BatchDynamicConnectivity() = default;

  /** The default constructor is invalid because the number of vertices in the
   *  graph must be known. */
  BatchDynamicConnectivity() = delete;

  /** Copy constructor not implemented. */
  BatchDynamicConnectivity(const BatchDynamicConnectivity &other) = delete;

  /** Copy assignment not implemented. */
  BatchDynamicConnectivity &operator=(const BatchDynamicConnectivity &other) =
      delete;

  /** Move constructor. */
  BatchDynamicConnectivity(BatchDynamicConnectivity &&other) noexcept;

  /** Move assignment not implemented. */
  BatchDynamicConnectivity &operator=(
      BatchDynamicConnectivity &&other) noexcept;

  /** Returns true if vertices \p u and \p v are connected in the graph.
   *
   *  Efficiency:
   *
   *  @param[in] suv A sequence of pair's of vertices
   */
  parlay::sequence<char> BatchConnected(
      parlay::sequence<std::pair<Vertex, Vertex>> suv) const;

  /** Returns true if edge \p edge is in the graph.
   *
   *  Efficiency: constant on average.
   *
   *  @param[in] edge Edge.
   *  @returns True if \p edge is in the graph, false if it is not.
   */
  // TODO(Sualeh): implement this
  // bool HasEdge(const UndirectedEdge &edge) const;

  /** Returns the number of vertices in `v`'s connected component.
   *
   * Efficiency: logarithmic in the size of the graph.
   *
   * @param[in] v Vertex.
   * @returns The number of vertices in \p v's connected component.
   */
  // TODO(Sualeh): implement this
  // int64_t GetSizeOfConnectedComponent(Vertex v) const;

  /** Adds an edge to the graph.
   *
   *  The edge must not already be in the graph and must not be a self-loop
   * edge.
   *
   *  Efficiency:
   *
   *  @param[in] se A sequence of edges
   */
  void BatchAddEdges(const parlay::sequence<UndirectedEdge> &se);

  /** Deletes an edge from the graph.
   *
   *  An exception will be thrown if the edge is not in the graph.
   *
   *  Efficiency:
   *
   *  @param[in] se A sequence of edges
   */
  void BatchDeleteEdges(parlay::sequence<UndirectedEdge> &se);

  parlay::sequence<Vertex> BatchFindRepr(const parlay::sequence<Vertex> &sv);

  void PrintStructure();
  void PrintLevel(int8_t level);
  void PrintNonTreeEdges();
  void PrintNonTreeEdgesForLevel(int8_t level);

 private:
  const int64_t num_vertices_;

  const int8_t max_level_;
  const detail::EdgeInfo empty_info = {-1, detail::EdgeType::kNonTree};
  const std::tuple<std::pair<Vertex, Vertex>, detail::EdgeInfo> empty_edge =
      std::make_tuple(std::make_pair(-1, -1), empty_info);

  // `spanning_forests_[i]` stores F_i, the spanning forest for the i-th
  // subgraph. In particular, `spanning_forests[0]` is a spanning forest for the
  // whole graph.
  parlay::sequence<BatchDynamicET *> parallel_spanning_forests_;

  // `non_tree_adjacency_lists_[i][v]` contains the vertices connected to vertex
  // v by level-i non-tree edges.
  // `non_tree_adjacency_lists_` is a vector of size `max_level_`. (i.e. the
  // first index represents the levels.)
  // `non_tree_adjacency_lists_[i]` is a vector of size `num_vertices_`. (i.e.
  // the second index represents the vertices.)
  parlay::sequence<parlay::sequence<vertexSet>> non_tree_adjacency_lists_;

  // TODO: use a concurrent map here.
  // All edges in the graph.
  // std::unordered_map<UndirectedEdge, detail::EdgeInfo, UndirectedEdgeHash>
  //     edges_;
  elektra::resizable_table<std::pair<Vertex, Vertex>, detail::EdgeInfo,
                           HashIntPairStruct>
      edges_;

  UndirectedEdge componentSearch(int level, Vertex v);

  treeSet getSpanningTree(const parlay::sequence<UndirectedEdge> &se);

  void replacementSearch(int level, parlay::sequence<int> components,
                         parlay::sequence<pair<int, int>> &promotedEdges);

  template <typename Rank, typename Parent>
  treeSet constructTree(Rank &r, Parent &p,
                        const parlay::sequence<UndirectedEdge> &se);

  parlay::sequence<int> removeDuplicates(parlay::sequence<int> &seq);
};

auto generateInitialVertexLayer(int numVertices, int max_level_) {
  auto vtxLayer =
      parlay::sequence<vertexSet>::from_function(numVertices, [&](int n) {
        return vertexSet(INITIAL_SIZE, parlay::hash_numeric<Vertex>{});
      });

  return vtxLayer;
}

BatchDynamicConnectivity::BatchDynamicConnectivity(int numVertices)
    : num_vertices_(numVertices), max_level_(parlay::log2_up(numVertices)) {
  parallel_spanning_forests_ = parlay::sequence<BatchDynamicET *>(max_level_);

  parlay::parallel_for(0, max_level_, [&](int i) {
    BatchDynamicET *ET = new BatchDynamicET{numVertices};
    parallel_spanning_forests_[i] = ET;
  });

  non_tree_adjacency_lists_ =
      parlay::sequence<parlay::sequence<vertexSet>>(max_level_);

  parlay::parallel_for(0, max_level_, [&](int i) {
    auto vtxLayer = generateInitialVertexLayer(numVertices, max_level_);
    non_tree_adjacency_lists_[i] = vtxLayer;
  });

  // make the base size smaller and update n_elms if we ever update.
  edges_ = elektra::resizable_table<pair<Vertex, Vertex>, detail::EdgeInfo,
                                    HashIntPairStruct>(
      num_vertices_ * num_vertices_, empty_edge, HashIntPairStruct());
}

BatchDynamicConnectivity::BatchDynamicConnectivity(
    int numVertices, const parlay::sequence<UndirectedEdge> &se)
    : num_vertices_(numVertices), max_level_(parlay::log2_up(numVertices)) {
  parallel_spanning_forests_ = parlay::sequence<BatchDynamicET *>(max_level_);

  parlay::parallel_for(0, max_level_, [&](int i) {
    BatchDynamicET *ET = new BatchDynamicET{numVertices};
    parallel_spanning_forests_[i] = ET;
  });

  non_tree_adjacency_lists_ =
      parlay::sequence<parlay::sequence<vertexSet>>(max_level_);

  parlay::parallel_for(0, max_level_, [&](int i) {
    auto vtxLayer = generateInitialVertexLayer(numVertices, max_level_);
    non_tree_adjacency_lists_[i] = vtxLayer;
  });

  edges_ = elektra::resizable_table<pair<Vertex, Vertex>, detail::EdgeInfo,
                                    HashIntPairStruct>(
      num_vertices_ * num_vertices_, empty_edge, HashIntPairStruct());

  BatchAddEdges(se);

  // auto maxLevelEulerTree = parallel_spanning_forests_[max_level_ - 1];
}

// TODO: fix the getRepresentative function
// parlay::sequence<Vertex> BatchDynamicConnectivity::BatchFindRepr(
//     const parlay::sequence<Vertex> &sv) {
//   auto pMaxLevelEulerTree = parallel_spanning_forests_[max_level_ - 1];
//   // get representatives from the highest level foredt
//   return parlay::map(sv, [&](Vertex v) {
//     return (int64_t)pMaxLevelEulerTree->getRepresentative(v);
//   });
// }

// TODO: get the concurrent union find for this
template <typename Rank, typename Parent>
treeSet BatchDynamicConnectivity::constructTree(
    Rank &r, Parent &p, const parlay::sequence<UndirectedEdge> &se) {
  // Given a sequence of edges, returns a set of a
  // spanning forest of the graph formed by them
  boost::disjoint_sets<Rank, Parent> dsu(r, p);
  treeSet tree;

  // BUG POSSIBLE: check whether this causes a problem.
  for (auto v : se) {
    dsu.make_set(v.first);
    dsu.make_set(v.second);
  }

  for (auto v : se) {
    Vertex first = v.first;
    Vertex second = v.second;
    // TODO is there a race condition here if we paralize this? How can we
    // resolve that

    if (dsu.find_set(first) != dsu.find_set(second)) {
      tree.insert(v);
      dsu.link(first, second);
    }
  }
  return tree;
}

// TODO: add parallel DSU structure to implement this
treeSet BatchDynamicConnectivity::getSpanningTree(
    const parlay::sequence<UndirectedEdge> &se) {
  // I am assuming the interface in
  // https://github.com/ParAlg/gbbs/blob/master/gbbs/union_find.h?fbclid=IwAR0U_Nbe1SpQF7mbmN0CEGLyF-5v362oy1q-9eQLvjQz916jhfTH69bMx9s
  // could be worth paralelizing this

  typedef std::map<Vertex, size_t> rank_t;
  typedef std::map<Vertex, Vertex> parent_t;

  rank_t rank_map;
  parent_t parent_map;

  boost::associative_property_map<rank_t> rank_pmap(rank_map);
  boost::associative_property_map<parent_t> parent_pmap(parent_map);

  return constructTree(rank_pmap, parent_pmap, se);
}

parlay::sequence<std::pair<int, int>> edgeBatchToPairArray(
    parlay::sequence<UndirectedEdge> &se) {
  // turns a sequence of edges to an array of pairs
  // useful for interfacing with EulerTourTrees
  auto array = parlay::sequence<std::pair<int, int>>::uninitialized(se.size());

  parlay::parallel_for(0, se.size(), [&](int i) {
    array[i] = std::make_pair(se[i].first, se[i].second);
  });
  return array;
}

// We print the structure
void BatchDynamicConnectivity::PrintStructure() {
  // print the structure of the graph
  // for each level, print the edges in the level

  std::cout << "\n\n -- Printing the structure of the ETT -- " << std::endl;
  std::cout << " -- --          tree edges         -- -- " << std::endl;
  for (int i = 0; i < max_level_; i++) {
    PrintLevel(i);
  }
  std::cout << " -- --       non tree edges         -- -- " << std::endl;
  PrintNonTreeEdges();
  std::cout << " -- done -- " << std::endl;
  std::cout << "\n\n" << std::endl;
}

// We print the edges in each level of our structure
void BatchDynamicConnectivity::PrintLevel(int8_t level) {
  // print the edges in the level
  auto pLevelEulerTree = parallel_spanning_forests_[level];

  // if (pLevelEulerTree == nullptr || level >= max_level_ || level < 0 ||
  //     pLevelEulerTree->IsEmpty()) {
  //   return;
  // }

  std::cout << "Level " << (int)level << ": " << std::endl;
  pLevelEulerTree->Print();
}

// We print the non tree edges
void BatchDynamicConnectivity::PrintNonTreeEdges() {
  // print the non tree edges
  for (int8_t i = 0; i < max_level_; i++) {
    PrintNonTreeEdgesForLevel(i);
  }
}

// We print the non tree edges in each level of our structure
void BatchDynamicConnectivity::PrintNonTreeEdgesForLevel(int8_t level) {
  std::cout << "Level " << (int)level << ": " << std::endl;
  // contains all the non-tree edges in the level
  auto vtxLayer = non_tree_adjacency_lists_[level];

  // scan through and build a set of all the edges
  std::unordered_set<UndirectedEdge, UndirectedEdgeHash> edges;

  for (int i = 0; i < num_vertices_; i++) {
    auto vtxSet = vtxLayer[i];
    for (auto v : vtxSet.entries()) {
      // order the insertion to avoid duplicates
      if (v > i) {
        edges.insert(UndirectedEdge(i, v));
      }
    }
  }

  // print the edges
  for (auto e : edges) {
    std::cout << "(" << e.first << ", " << e.second << "), ";
  }
  std::cout << std::endl;
}

}  // namespace batchDynamicConnectivity

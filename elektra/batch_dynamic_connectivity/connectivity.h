#pragma once

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
#include "parallel_euler_tour_tree/euler_tour_tree.hpp"

namespace detail {

typedef int8_t Level;

enum class EdgeType {
  // Edge is in the spanning forest of the graph.
  kNonTree,
  // Edge is not in the spanning forest of the graph.
  kTree,
};

// A struct that contains information about a particular edge.
struct EdgeInfo {
  Level level;
  EdgeType type;
};

}  // namespace detail

typedef int64_t Vertex;

/** This class represents an undirected graph that can undergo efficient edge
 *  insertions, edge deletions, and connectivity queries. Multiple edges between
 *  a pair of vertices are supported.
 */

namespace batchDynamicConnectivity {

using UndirectedEdge = dynamicGraph::UndirectedEdge;
using UndirectedEdgeHash = dynamicGraph::UndirectedEdgeHash;
using BatchDynamicET = parallel_euler_tour_tree::EulerTourTree;

using treeSet = std::unordered_set<UndirectedEdge, UndirectedEdgeHash>;

class BatchDynamicConnectivity {
 public:
  /** Initializes an empty graph with a fixed number of vertices.
   *
   *  Efficiency: \f$ O(n \log n ) \f$ where \f$ n \f$ is the number of vertices
   *  in the graph.
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

  // TODO:make the API use sequences for everything.
  /** Returns true if vertices \p u and \p v are connected in the graph.
   *
   *  Efficiency: logarithmic in the size of the graph.
   *
   *  @param[in] u Vertex.
   *  @param[in] v Vertex.
   *  @returns True if \p u and \p v are connected, false if they are not.
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
  bool HasEdge(const UndirectedEdge &edge) const;

  /** Returns the number of vertices in `v`'s connected component.
   *
   * Efficiency: logarithmic in the size of the graph.
   *
   * @param[in] v Vertex.
   * @returns The number of vertices in \p v's connected component.
   */
  int64_t GetSizeOfConnectedComponent(Vertex v) const;

  /** Adds an edge to the graph.
   *
   *  The edge must not already be in the graph and must not be a self-loop
   * edge.
   *
   *  Efficiency: \f$ O\left( \log^2 n \right) \f$ amortized where \f$ n \f$ is
   *  the number of vertices in the graph.
   *
   *  @param[in] edge Edge to be added.
   */
  void BatchAddEdges(const parlay::sequence<UndirectedEdge> &se);

  /** Deletes an edge from the graph.
   *
   *  An exception will be thrown if the edge is not in the graph.
   *
   *  Efficiency: \f$ O\left( \log^2 n \right) \f$ amortized where \f$ n \f$ is
   *  the number of vertices in the graph.
   *
   *  @param[in] edge Edge to be deleted.
   */
  void BatchDeleteEdges(const parlay::sequence<UndirectedEdge> &se);

  parlay::sequence<Vertex> BatchFindRepr(const parlay::sequence<Vertex> &sv);

  void printStructure();
  void printLevel(int8_t level);
  void printNonTreeEdges();
  void printNonTreeEdgesForLevel();

 private:
  const int64_t num_vertices_;

  // TODO: convert this to int8_t
  const int8_t max_level_;

  // `spanning_forests_[i]` stores F_i, the spanning forest for the i-th
  // subgraph. In particular, `spanning_forests[0]` is a spanning forest for the
  // whole graph.
  parlay::sequence<BatchDynamicET *> parallel_spanning_forests_;

  // `adjacency_lists_by_level_[i][v]` contains the vertices connected to vertex
  // v by level-i non-tree edges.
  // `adjacency_lists_by_level_` is a vector of size `max_level_`. (i.e. the
  // first index represents the levels.)
  // `adjacency_lists_by_level_[i]` is a vector of size `num_vertices_`. (i.e.
  // the second index represents the vertices.)
  // TODO: make this concurrent map
  parlay::sequence<parlay::sequence<std::unordered_set<Vertex>>>
      non_tree_adjacency_lists_;

  // TODO: use a concurrent map here.
  // All edges in the graph.
  std::unordered_map<UndirectedEdge, detail::EdgeInfo, UndirectedEdgeHash>
      edges_;

  void AddNonTreeEdge(const UndirectedEdge &edge);

  void BatchAddNonTreeEdge(const parlay::sequence<UndirectedEdge> &se);

  void AddTreeEdge(const UndirectedEdge &edge);

  void BatchAddTreeEdge(const parlay::sequence<UndirectedEdge> &se);

  void AddEdgeToAdjacencyList(const UndirectedEdge &edge, detail::Level level);

  void BatchUpdateAdjacencyList(
      const parlay::sequence<std::pair<UndirectedEdge, detail::Level>> &sel);

  void DeleteEdgeFromAdjacencyList(const UndirectedEdge &edge,
                                   detail::Level level);

  void BatchDeleteEdgesInAdjacencyList(
      const parlay::sequence<std::pair<UndirectedEdge, detail::Level>> &sel);

  void ReplaceTreeEdge(const UndirectedEdge &edge, detail::Level level);

  UndirectedEdge componentSearch(int level, Vertex v);

  parlay::sequence<Vertex> parallelLevelSearch(
      const parlay::sequence<UndirectedEdge> &se,
      parlay::sequence<Vertex> &components,
      parlay::sequence<std::pair<int, int>> &promotedEdges, int level);

  treeSet getSpanningTree(const parlay::sequence<UndirectedEdge> &se);

  void replacementSearch(
      int level, parlay::sequence<int> components,
      parlay::sequence<pair<int, int>> &promotedEdges);

  template <typename Rank, typename Parent>
  treeSet constructTree(Rank &r, Parent &p,
                        const parlay::sequence<UndirectedEdge> &se);

  auto removeDuplicates(parlay::sequence<int> &seq);
};

parlay::sequence<std::unordered_set<Vertex>> generateInitialVertexLayer(
    int numVertices, int max_level_) {
  auto vtxLayer = parlay::sequence<std::unordered_set<Vertex>>(numVertices);

  parlay::parallel_for(0, numVertices, [&](int i) {
    auto vtxset = std::unordered_set<Vertex>();
    vtxLayer[i] = vtxset;
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
      parlay::sequence<parlay::sequence<std::unordered_set<Vertex>>>(
          max_level_);

  parlay::parallel_for(0, max_level_, [&](int i) {
    auto vtxLayer = generateInitialVertexLayer(numVertices, max_level_);
    non_tree_adjacency_lists_[i] = vtxLayer;
  });

  edges_ = std::unordered_map<UndirectedEdge, detail::EdgeInfo,
                              UndirectedEdgeHash>();
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
      parlay::sequence<parlay::sequence<std::unordered_set<Vertex>>>(
          max_level_);

  parlay::parallel_for(0, max_level_, [&](int i) {
    auto vtxLayer = generateInitialVertexLayer(numVertices, max_level_);
    non_tree_adjacency_lists_[i] = vtxLayer;
  });
  edges_ = std::unordered_map<UndirectedEdge, detail::EdgeInfo,
                              UndirectedEdgeHash>();

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

}  // namespace batchDynamicConnectivity

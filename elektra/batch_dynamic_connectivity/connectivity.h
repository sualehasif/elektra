#pragma once

#include "connectivity-helpers.h"

#define INITIAL_SIZE 50

/** This class represents an undirected graph that can undergo efficient edge
 *  insertions, edge deletions, and connectivity queries. Multiple edges between
 *  a pair of vertices are supported.
 *
 *  TODO(sualeh): What happens if we perform a delete on an edge
 *  (u,v) and multiple (u,v) edges are present? Do we delete one of
 *  them? All of them? It might be cleaner to assume/enforce that
 *  there's only a single edge between a pair of vertices.
 */

namespace bdcty {

class BatchDynamicConnectivity {
 public:
  /** Initializes an empty graph with a fixed number of vertices.
   *
   *  @param[in] num_vertices Number of vertices in the graph.
   *
   *  TODO(sualeh): What if the number of vertices is larger than
   *  INT_MAX? How about having a type (smt like uintV) which
   *  we can set to either uint32_t or uint64_t, depending on the
   *  need?
   *  Maybe make the "int" type here whatever we call "V" in
   *  graph.h?
   */
  explicit BatchDynamicConnectivity(int num_vertices);

  explicit BatchDynamicConnectivity(int num_vertices,
                                    const parlay::sequence<E> &se);

  /** Deallocates the data structure. */
  ~BatchDynamicConnectivity() = default;

  /** The default constructor is invalid because the number of vertices in the
   *  graph must be known. */
  BatchDynamicConnectivity() = delete;

  /** Copy constructor not implemented. */
  BatchDynamicConnectivity(const BatchDynamicConnectivity &other) = delete;

  /** Copy assignment not implemented. */
  auto operator=(const BatchDynamicConnectivity &other)
  -> BatchDynamicConnectivity & = delete;

  /** Move constructor. */
  BatchDynamicConnectivity(BatchDynamicConnectivity &&other) = delete;

  /** Move assignment not implemented. */
  auto operator=(BatchDynamicConnectivity &&other) noexcept
  -> BatchDynamicConnectivity & = delete;

  /** Returns true if vertices \p u and \p v are connected in the graph.
   *
   *  Efficiency:
   *
   *  @param[in] suv A sequence of pair's of vertices
   *
   *  TODO(sualeh): minor, but why is return type is char and not bool?
   */
  [[nodiscard]] auto BatchConnected(parlay::sequence<std::pair<V, V>> suv) const
  -> parlay::sequence<char>;

  /** Adds a batch of edges to the graph.
   *
   *  @param[in] se A sequence of edges
   */
  void BatchAddEdges(const parlay::sequence<E> &se);

  /** Deletes a batch of edges from the graph.
   *
   *  An exception will not be thrown if the edge is not in the graph.
   *
   *  Efficiency:
   *
   *  @param[in] se A sequence of edges
   */
  void BatchDeleteEdges(parlay::sequence<E> &se);

  //  parlay::sequence<V> BatchFindRepr(const parlay::sequence<V> &sv);

  [[maybe_unused]] void PrintStructure();
  void PrintLevel(int8_t level);

  void PrintNonTreeEdges();

  void PrintNonTreeEdgesForLevel(int8_t level);

 private:
  const int64_t num_vertices_;

  const Level max_level_;
  const bdcty::EInfo empty_info_ = {-1, bdcty::EType::K_NON_TREE};
  const std::tuple<std::pair<V, V>, bdcty::EInfo> empty_edge_ =
      std::make_tuple(std::make_pair(-1, -1), empty_info_);
  const std::tuple<std::pair<V, V>, bdcty::EInfo> tombstone_edge_ =
      std::make_tuple(std::make_pair(-1, -1), empty_info_);

  // `spanning_forests_[i]` stores F_i, the spanning forest for the i-th
  // subgraph. In particular, `spanning_forests[0]` is a spanning forest for the
  // whole graph.
  parlay::sequence<BatchDynamicEtt *> parallel_spanning_forests_;

  // `non_tree_adjacency_lists_[i][v]` contains the vertices connected to vertex
  // v by level-i non-tree edges.
  // `non_tree_adjacency_lists_` is a vector of size `max_level_`. (i.e. the
  // first index represents the levels.)
  // `non_tree_adjacency_lists_[i]` is a vector of size `num_vertices_`. (i.e.
  // the second index represents the vertices.)
  NonTreeAdjacencyList non_tree_adjacency_lists_;

  // TODO: use a concurrent map here.
  // All edges in the graph.
  // std::unordered_map<E, batchDynamicConnectivity::EInfo, EHash>
  //     edges_;
  elektra::resizable_table<std::pair<V, V>, bdcty::EInfo, HashIntPairStruct>
      edges_;

  auto GetSpanningTree(const parlay::sequence<E> &se) -> TreeSet;

  void ReplacementSearch(int level, parlay::sequence<int> components,
                         parlay::sequence<pair<int, int>> &promoted_edges);

  template<typename Rank, typename Parent>
  auto ConstructTree(Rank &r, Parent &p, const parlay::sequence<E> &se)
  -> TreeSet;

  static auto RemoveDuplicates(parlay::sequence<int> &seq)
  -> parlay::sequence<int>;
};

// auto generateInitialVertexLayer(int numVertices, int max_level_) {
//   auto vtxLayer =
//       parlay::sequence<vertexSet>::from_function(numVertices, [&](int n) {
//         return vertexSet(INITIAL_SIZE, parlay::hash_numeric<V>{});
//       });

//   return vtxLayer;
// }

BatchDynamicConnectivity::BatchDynamicConnectivity(int num_vertices)
    : num_vertices_(num_vertices), max_level_(parlay::log2_up(num_vertices)) {
  parallel_spanning_forests_ = parlay::sequence<BatchDynamicEtt *>(max_level_);

  parlay::parallel_for(0, max_level_, [&](int i) {
    auto *ett = new BatchDynamicEtt{num_vertices};
    parallel_spanning_forests_[i] = ett;
  });

  non_tree_adjacency_lists_ = NonTreeAdjacencyList(num_vertices, max_level_);

  // make the base size smaller and update n_elms if we ever update.
  edges_ =
      elektra::resizable_table<pair<V, V>, bdcty::EInfo, HashIntPairStruct>(
          num_vertices_ * num_vertices_, empty_edge_, tombstone_edge_,
          HashIntPairStruct());
}

BatchDynamicConnectivity::BatchDynamicConnectivity(
    int num_vertices, const parlay::sequence<E> &se)
    : num_vertices_(num_vertices), max_level_(parlay::log2_up(num_vertices)) {
  parallel_spanning_forests_ = parlay::sequence<BatchDynamicEtt *>(max_level_);

  parlay::parallel_for(0, max_level_, [&](int i) {
    auto *ett = new BatchDynamicEtt{num_vertices};
    parallel_spanning_forests_[i] = ett;
  });

  non_tree_adjacency_lists_ = NonTreeAdjacencyList(num_vertices, max_level_);

  edges_ =
      elektra::resizable_table<pair<V, V>, bdcty::EInfo, HashIntPairStruct>(
          num_vertices_ * num_vertices_, empty_edge_, tombstone_edge_,
          HashIntPairStruct());

  BatchAddEdges(se);
}

// TODO: fix the getRepresentative function
// parlay::sequence<V> BatchDynamicConnectivity::BatchFindRepr(
//     const parlay::sequence<V> &sv) {
//   auto pMaxLevelEulerTree = parallel_spanning_forests_[max_level_ - 1];
//   // get representatives from the highest level for ett
//   return parlay::map(sv, [&](V v) {
//     return (int64_t)pMaxLevelEulerTree->getRepresentative(v);
//   });
// }

// TODO: get the concurrent union find for this
template<typename Rank, typename Parent>
auto BatchDynamicConnectivity::ConstructTree(Rank &r, Parent &p,
                                             const parlay::sequence<E> &se)
-> TreeSet {
  // Given a sequence of edges, returns a set of a
  // spanning forest of the graph formed by them
  boost::disjoint_sets<Rank, Parent> dsu(r, p);
  TreeSet tree;

  // BUG POSSIBLE: check whether this causes a problem.
  for (auto v : se) {
    dsu.make_set(v.first);
    dsu.make_set(v.second);
  }

  for (auto v : se) {
    V first = v.first;
    V second = v.second;
    // TODO is there a race condition here if we parallelize this? How can we
    // resolve that

    if (dsu.find_set(first) != dsu.find_set(second)) {
      tree.insert(v);
      dsu.link(first, second);
    }
  }
  return tree;
}

// TODO: add parallel DSU structure to implement this
auto BatchDynamicConnectivity::GetSpanningTree(const parlay::sequence<E> &se)
-> TreeSet {
  // I am assuming the interface in
  // https://github.com/ParAlg/gbbs/blob/master/gbbs/union_find.h?fbclid=IwAR0U_Nbe1SpQF7mbmN0CEGLyF-5v362oy1q-9eQLvjQz916jhfTH69bMx9s
  // could be worth it to parallelize this

  using RankT = std::map<V, size_t>;
  using ParentT = std::map<V, V>;

  RankT rank_map;
  ParentT parent_map;

  boost::associative_property_map<RankT> rank_pmap(rank_map);
  boost::associative_property_map<ParentT> parent_pmap(parent_map);

  return ConstructTree(rank_pmap, parent_pmap, se);
}

auto EdgeBatchToPairArray(parlay::sequence<E> &se)
-> parlay::sequence<std::pair<int, int>> {
  // turns a sequence of edges to an array of pairs
  // useful for interfacing with EulerTourTrees
  auto array = parlay::sequence<std::pair<int, int>>::uninitialized(se.size());

  parlay::parallel_for(0, se.size(), [&](int i) {
    array[i] = std::make_pair(se[i].first, se[i].second);
  });
  return array;
}

// We print the structure
[[maybe_unused]] void BatchDynamicConnectivity::PrintStructure() {
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
  auto p_level_euler_tree = parallel_spanning_forests_[level];

  // if (p_level_euler_tree == nullptr || level >= max_level_ || level < 0 ||
  //     p_level_euler_tree->IsEmpty()) {
  //   return;
  // }

  std::cout << "Level " << (int) level << ": " << std::endl;
  p_level_euler_tree->Print();
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
  std::cout << "Level " << (int) level << ": " << std::endl;
  // contains all the non-tree edges in the level
  //
  // TODO(sualeh) Note that non_tree_adjacency_lists[level] is a
  // parlay::sequence and so the below will create a copy. Instead you
  // can do:
  // auto& vtx_layer = ...
  auto vtx_layer = non_tree_adjacency_lists_[level];

  // scan through and build a set of all the edges
  std::unordered_set<E, EHash> edges;

  for (int i = 0; i < num_vertices_; i++) {
    auto vtx_set = vtx_layer[i];
    for (auto &kv : vtx_set.entries()) {
      // get the vertex by getting the first element of the tuple
      auto v = std::get<0>(kv);
      // order the insertion to avoid duplicates
      if (v > i) {
        edges.insert(make_pair(i, v));
      }
    }
  }

  // print the edges
  for (auto e : edges) {
    std::cout << "(" << e.first << ", " << e.second << "), ";
  }
  std::cout << std::endl;
}

} // namespace bdcty

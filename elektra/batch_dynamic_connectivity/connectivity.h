#pragma once

#include "connectivity-helpers.h"

#define INITIAL_SIZE 50

namespace bdcty {
/** This class represents an undirected graph that can undergo efficient edge
 *  insertions, edge deletions, and connectivity queries. Multiple edges between
 *  a pair of vertices are not supported.
 */
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
  explicit BatchDynamicConnectivity(V num_vertices);

  explicit BatchDynamicConnectivity(V num_vertices,
                                    const parlay::sequence<E> &se);

  /** Deallocates the data structure.
   * ~BatchDynamicConnectivity() = default;

   * The default constructor is invalid because the number of vertices in the
   *  graph must be known.
  */
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
  void PrintLevel(Level level);

  void PrintNonTreeEdges();

  void PrintNonTreeEdgesForLevel(Level level);

 private:
  const uintE num_vertices_;

  const Level max_level_;
  const std::tuple<std::pair<V, V>, bdcty::EInfo> empty_edge_ =
      std::make_tuple(std::make_pair(-1, -1), kEmptyInfo);
  const std::tuple<std::pair<V, V>, bdcty::EInfo> tombstone_edge_ =
      std::make_tuple(std::make_pair(-1, -1), kEmptyInfo);
  const std::tuple<std::pair<V, V>, elektra::empty> edge_with_empty_struct_ =
      std::make_tuple(std::make_pair(-1, -1), elektra::empty{});
  const std::tuple<std::pair<V, V>, elektra::empty>
      tombstone_with_empty_struct_ =
          std::make_tuple(std::make_pair(-1, -1), elektra::empty{});

  // `spanning_forests_[i]` stores F_i, the spanning forest for the i-th
  // subgraph. In particular, `spanning_forests[0]` is a spanning forest for the
  // whole graph.
  parlay::sequence<unique_ptr<BatchDynamicEtt>> parallel_spanning_forests_;

  // `non_tree_adjacency_lists_[i][v]` contains the vertices connected to vertex
  // v by level-i non-tree edges.
  // `non_tree_adjacency_lists_` is a vector of size `max_level_`. (i.e. the
  // first index represents the levels.)
  // `non_tree_adjacency_lists_[i]` is a vector of size `num_vertices_`. (i.e.
  // the second index represents the vertices.)
  NonTreeAdjacencyList non_tree_adjacency_lists_;
  EdgeSet edges_;

  void ReplacementSearch(Level level, parlay::sequence<V> components,
                         parlay::sequence<pair<V, V>> &promoted_edges);

  void PushDownTreeEdgesFromComponents(Level l, parlay::sequence<V> &components);
  void PushDownNonTreeEdges(Level l, parlay::sequence<E> &non_tree_edges);

  static auto RemoveDuplicates(parlay::sequence<V> &seq) -> parlay::sequence<V>;
  static auto RemoveDuplicates(parlay::sequence<V> &&seq)
  -> parlay::sequence<V>;
  inline void InsertIntoEdgeTable(const pair<V, V> &e, EType e_type,
                                  Level level);
  inline auto GetSearchEdges(Level level, sequence<V> &components,
                             const unique_ptr<BatchDynamicEtt> &ett)
  -> vector<vector<E>>;
};

BatchDynamicConnectivity::BatchDynamicConnectivity(V num_vertices)
    : num_vertices_(num_vertices), max_level_(parlay::log2_up(num_vertices)) {
  parallel_spanning_forests_ =
      parlay::sequence<unique_ptr<BatchDynamicEtt>>::from_function(
          max_level_, [&](size_t i) {
            return make_unique<BatchDynamicEtt>(num_vertices_);
          });

  // parlay::parallel_for(0, max_level_, [&](auto i) {
  //   parallel_spanning_forests_[i] =
  //   make_unique<BatchDynamicEtt>(num_vertices_);
  // });

  non_tree_adjacency_lists_ = NonTreeAdjacencyList(num_vertices, max_level_);

  // make the base size smaller and update n_elms if we ever update.
  edges_ =
      elektra::resizable_table<pair<V, V>, bdcty::EInfo, HashIntPairStruct>(
          static_cast<size_t>(num_vertices_ * num_vertices_), empty_edge_,
          tombstone_edge_, HashIntPairStruct());
}

BatchDynamicConnectivity::BatchDynamicConnectivity(
    V num_vertices, const parlay::sequence<E> &se)
    : num_vertices_(num_vertices), max_level_(parlay::log2_up(num_vertices)) {
  parallel_spanning_forests_ =
      parlay::sequence<unique_ptr<BatchDynamicEtt>>::from_function(
          max_level_, [&](size_t i) {
            return make_unique<BatchDynamicEtt>(num_vertices_);
          });

  // parlay::parallel_for(0, max_level_, [&](auto i) {
  //   parallel_spanning_forests_[i] =
  //   make_unique<BatchDynamicEtt>(num_vertices_);
  // });

  non_tree_adjacency_lists_ = NonTreeAdjacencyList(num_vertices, max_level_);

  edges_ =
      elektra::resizable_table<pair<V, V>, bdcty::EInfo, HashIntPairStruct>(
          static_cast<size_t>(num_vertices_ * num_vertices_), empty_edge_,
          tombstone_edge_, HashIntPairStruct());

  BatchAddEdges(se);
}

template<typename E, typename N = E>
auto RepresentativeSpanningTree(const parlay::sequence<E> &se,
                                const unique_ptr<BatchDynamicEtt> &ett) {
  // Construct the auxiliary edges.
  sequence<pair<N, N>> aux_int_edges = parlay::map(se, [&](E e) {
    return make_pair(static_cast<N>(ett->GetRepresentative(e.first)),
                     static_cast<N>(ett->GetRepresentative(e.second)));
  });

  auto spanning_tree = sequence<pair<N, N>>::uninitialized(se.size());
  elektra::SpanningTree(aux_int_edges, spanning_tree);

  return spanning_tree;
}

template<typename T, class VType>
auto NewEdgeSequence(parlay::sequence<pair<VType, VType>> &se)
-> parlay::sequence<std::pair<T, T>> {
  // turns a sequence of edges to an array of pairs
  // useful for V with EulerTourTrees
  auto array = parlay::sequence<std::pair<T, T>>::uninitialized(se.size());

  parlay::parallel_for(0, se.size(), [&](auto i) {
    array[i] = std::make_pair(static_cast<T>(se[i].first),
                              static_cast<T>(se[i].second));
  });
  return array;
}

inline void BatchDynamicConnectivity::InsertIntoEdgeTable(const pair<V, V> &e,
                                                          EType e_type,
                                                          Level level) {
  EInfo ei = {level, e_type};
  edges_.insert(make_tuple(pair<V, V>(e.first, e.second), ei));

  // TODO(sualeh): Think about whether you can get away without
  // inserting the reverse edge add the reverse edge to the edges_ map
  EInfo ei_rev = {level, e_type};
  edges_.insert(make_tuple(pair<V, V>(e.second, e.first), ei_rev));
}

// TODO(tom): This needs to be supported as an augmentation.
inline auto BatchDynamicConnectivity::GetSearchEdges(
    Level level, sequence<V> &components, const unique_ptr<BatchDynamicEtt> &ett)
-> vector<vector<E>> {
  vector<vector<E>> non_tree_edges;
  non_tree_edges.reserve(components.size());
  for (auto _ : components) {
    non_tree_edges.emplace_back();
  }
  for (uintE i = 0; i < components.size(); i++) {
    auto comp_id = components[i];
    auto cc = ett->ComponentVertices(comp_id);
    // for each vertex in the component
    for (V u : cc) {
      // for each edge in the component
      for (auto &kw : non_tree_adjacency_lists_[level][u].entries()) {
        auto w = std::get<0>(kw);
        non_tree_edges[i].push_back(E(u, w));
      }
    }
  }
  return non_tree_edges;
}

auto BatchDynamicConnectivity::RemoveDuplicates(sequence<V> &seq)
-> sequence<V> {
  // TODO: possibly change this to use a semisort and not inplace sort
  // sort the sequence
  parlay::integer_sort_inplace(seq, [](V x) { return (unsigned) x; });

  return parlay::unique(seq);
}

auto BatchDynamicConnectivity::RemoveDuplicates(sequence<V> &&seq)
-> sequence<V> {
  parlay::integer_sort_inplace(seq, [](V x) { return (unsigned) x; });
  return parlay::unique(seq);
}

// We print the structure
[[maybe_unused]] void BatchDynamicConnectivity::PrintStructure() {
  // print the structure of the graph
  // for each level, print the edges in the level

  std::cout << "\n\n -- Printing the structure of the ETT -- " << std::endl;
  std::cout << " -- --          tree edges         -- -- " << std::endl;
  for (auto i = 0; i < max_level_; i++) {
    PrintLevel(i);
  }
  std::cout << " -- --       non tree edges         -- -- " << std::endl;
  PrintNonTreeEdges();
  std::cout << " -- done -- " << std::endl;
  std::cout << "\n\n" << std::endl;
}

// We print the edges in each level of our structure
void BatchDynamicConnectivity::PrintLevel(Level level) {
  // print the edges in the level
  auto &p_level_euler_tree = parallel_spanning_forests_[level];

  // if (p_level_euler_tree == nullptr || level >= max_level_ || level < 0 ||
  //     p_level_euler_tree->IsEmpty()) {
  //   return;
  // }

  std::cout << "Level " << level << ": " << std::endl;
  p_level_euler_tree->Print();
}

// We print the non tree edges
void BatchDynamicConnectivity::PrintNonTreeEdges() {
  // print the non tree edges
  for (auto i = 0; i < max_level_; i++) {
    PrintNonTreeEdgesForLevel(i);
  }
}

// We print the non tree edges in each level of our structure
void BatchDynamicConnectivity::PrintNonTreeEdgesForLevel(Level level) {
  std::cout << "Level " << level << ": " << std::endl;
  // contains all the non-tree edges in the level
  //
  // TODO(sualeh) Note that non_tree_adjacency_lists[level] is a
  // parlay::sequence and so the below will create a copy. Instead you
  // can do:
  // auto& vtx_layer = ...
  auto vtx_layer = non_tree_adjacency_lists_[level];

  // scan through and build a set of all the edges
  std::unordered_set<E, EHash> edges;

  for (int i = 0; i < (int) num_vertices_; i++) {
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

}  // namespace bdcty

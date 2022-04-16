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
   */
  explicit BatchDynamicConnectivity(V num_vertices);

  explicit BatchDynamicConnectivity(V num_vertices,
                                    const parlay::sequence<E> &se);

  /** Returns true if vertices \p u and \p v are connected in the graph.
   *
   *  Efficiency:
   *
   *  @param[in] suv A sequence of pair's of vertices
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

  /** Printing utilities. */
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

  //  void checkRep();

  parlay::sequence<E> ReplacementSearch(Level level,
                                        parlay::sequence<V> components);

  void PushDownTreeEdgesFromComponents(Level l,
                                       parlay::sequence<V> &components);
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

// void BatchDynamicConnectivity::checkRep() {
//   // Rep invariant:
//
//   // Basic checks
//   // 1. `edges_` is a set of edges.
//   // 2. `spanning_forests_[i].edges_` is a subset of `edges_`.
//   // 3. `(v, non_tree_adjacency_lists_[i][v])` is a subset of `edges_`
//   //    for all i and v.
//   // 4. `Set( ...E(v, non_tree_adjacency_lists_[i][v]),
//   //          ...spanning_forests_[j].edges_)`
//   //   = `edges_`.    for all i, v, j.
//   // 5. `spanning_forests_[i].edges_` is a subset of
//   //    `spanning_forests_[j].edges_`    for all i > j.
//   // 6. all edges in `edges_`, `spanning_forests_[i].edges_` are distinct.
//
//   // Component Size checks
//   // 1. All components of edges at level i are of size <= 2^i.
//   //    for all i in [0, max_level_].
//
//   // MST checks
//   // 1. M = {WeightedEdge( (e, w) | e in spanning_forests_[max_level_].edges_
//   //                                 and  w = edges_[e].level}.
//   //   is a minimum spanning forest of `edges_`.
//
//   auto edges_seq = edges_.entries();
//
//   for (Level level = 0; level <= max_level_; ++level) {
//     // Check that `spanning_forests_[i].edges_` is a subset of `edges_`.
//     auto spanning_forest_edges = parallel_spanning_forests_[level]->Edges_();
//     assert(spanning_forest_edges.size() <= edges_seq.size());
//     assert(
//         std::count_if(edges_seq.begin(), edges_seq.end(), [&](const auto &e)
//         {
//           auto [edge, value] = e;
//           auto [edge_level, e_type] = value;
//           return edge_level == level && e_type == EType::K_TREE;
//         }) == spanning_forest_edges.size());
//
//     // Check that `(v, non_tree_adjacency_lists_[i][v])` is a subset of
//     // `edges_`.
//     auto non_tree_level_edges = non_tree_adjacency_lists_[level];
//     for (V v = 0; v < num_vertices_; ++v) {
//       auto opposite_vertices = non_tree_level_edges[v].entries();
//       // construct a set of edges from the opposite vertices.
//       auto opposite_edges = vector<pair<V, V>>(opposite_vertices.size());
//       for (size_t i = 0; i < opposite_vertices.size(); ++i) {
//         auto [opposite_vertex, _] = opposite_vertices[i];
//         opposite_edges[i] = std::make_pair(v, opposite_vertex);
//       }
//       assert(opposite_edges.size() <= edges_seq.size());
//       assert(
//           std::count_if(edges_seq.begin(), edges_seq.end(), [&](const auto
//           &e) {
//             auto [edge, value] = e;
//             auto [edge_level, e_type] = value;
//             return edge_level == level && e_type == EType::K_NON_TREE;
//           }) == opposite_edges.size());
//     }
//
//     //    // Check that `Set( ...E(v, non_tree_adjacency_lists_[i][v]),
//     //    //          ...spanning_forests_[j].edges_)`
//     //    //   = `edges_`.
//     //    for (V v = 0; v < num_vertices_; ++v) {
//     //      for (const auto &e : non_tree_adjacency_lists_[level][v]) {
//     //        assert(spanning_forests_[level]->edges_.count(e));
//     //      }
//     //    }
//     //
//     //    // Check that `spanning_forests_[i].edges_` is a subset of
//     //    // `spanning_forests_[j].edges_`    for all i > j.
//     //    for (Level j = level + 1; j <= max_level_; ++j) {
//     //      assert(spanning_forests_[level]->edges_.size() <=
//     //             spanning_forests_[j]->edges_.size());
//     //      assert(spanning_forests_[level]->edges_.size() ==
//     //             std::count_if(spanning_forests_[j]->edges_.begin(),
//     //                           spanning_forests_[j]->edges_.end(),
//     //                           [&](const auto &e) {
//     //                             return
//     //                             spanning_forests_[level]->edges_.count(e);
//     //                           }));
//   }
//
//   // component size checks
//   for (Level level = 0; level <= max_level_; ++level) {
//     // construct the components from the tree edges
//     auto level_edges = parallel_spanning_forests_[level]->Edges_();
//
//     auto components = vector<set<V>>(num_vertices_);
//
//     for (const auto &e : level_edges) {
//       auto [u, v] = e;
//
//       for (const auto &component : components) {
//         if (component.count(u)) {
//           components[u].insert(v);
//         } else if (component.count(v)) {
//           components[v].insert(u);
//         } else {
//           auto new_component = set<V>();
//           new_component.insert(u);
//           new_component.insert(v);
//           components.push_back(new_component);
//         }
//       }
//     }
//
//     // now insert the non-tree edges
//     auto non_tree_edges = non_tree_adjacency_lists_[level];
//
//     for (V v = 0; v < num_vertices_; ++v) {
//       for (const auto &e : non_tree_edges[v].entries()) {
//         auto [u, _] = e;
//         for (const auto &component : components) {
//           if (component.count(u)) {
//             components[u].insert(v);
//           } else if (component.count(v)) {
//             components[v].insert(u);
//           } else {
//             auto new_component = set<V>();
//             new_component.insert(u);
//             new_component.insert(v);
//             components.push_back(new_component);
//           }
//         }
//       }
//     }
//
//     // assert that the size of each component is at most 2^(level)
//     for (const auto &component : components) {
//       assert(component.size() <= (1 << level));
//     }
//   }
//
//   // MST checks
//
// }

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

template <typename E, typename N = E>
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

template <typename T, class VType>
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
inline auto
BatchDynamicConnectivity::GetSearchEdges(Level level, sequence<V> &components,
                                         const unique_ptr<BatchDynamicEtt> &ett)
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
  parlay::integer_sort_inplace(seq, [](V x) { return (unsigned)x; });

  return parlay::unique(seq);
}

auto BatchDynamicConnectivity::RemoveDuplicates(sequence<V> &&seq)
    -> sequence<V> {
  parlay::integer_sort_inplace(seq, [](V x) { return (unsigned)x; });
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

  for (V i = 0; i < num_vertices_; i++) {
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

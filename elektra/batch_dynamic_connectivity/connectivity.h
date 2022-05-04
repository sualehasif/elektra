#pragma once

#include "connectivity-helpers.h"
#include "log.h"

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
  [[nodiscard]] auto BatchConnected(parlay::sequence<std::pair<V, V>> suv)
      -> parlay::sequence<char>;

  /** Adds a batch of edges to the graph.
   *
   *  @param[in] se A sequence of edges
   */
  void BatchAddEdges(const parlay::sequence<E> &se);
  // void BatchAddEdges(const parlay::slice<E, E> &se);

  /** Deletes a batch of edges from the graph.
   *
   *  An exception will not be thrown if the edge is not in the graph.
   *
   *  Efficiency:
   *
   *  @param[in] se A sequence of edges
   */
  void BatchDeleteEdges(const parlay::sequence<E> &se);

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
      std::make_tuple(std::make_pair(-2, -2), kEmptyInfo);
  const std::tuple<std::pair<V, V>, elektra::empty> edge_with_empty_struct_ =
      std::make_tuple(std::make_pair(-1, -1), elektra::empty{});
  const std::tuple<std::pair<V, V>, elektra::empty>
      tombstone_with_empty_struct_ =
          std::make_tuple(std::make_pair(-2, -2), elektra::empty{});

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

  void CheckRep();

  // Returns
  //   - a sequence E[] of level-i non-tree edges that should be promoted to be tree edges
  //   - a subset of E[] representing all elements of E[] that should simultaneously be
  //     pushed to level i - 1
  //
  // It is the caller's job to do the promotion of E[] and to push the subset,
  // though this function takes care of pushing down non-promoted edges.
  auto ReplacementSearch(Level level, parlay::sequence<V> components)
      -> std::pair<sequence<E>, elektra::resizable_table<E, elektra::empty, HashIntPairStruct>>;

  void PushDownTreeEdgesFromComponents(Level l,
                                       parlay::sequence<V> &components);
  void PushDownNonTreeEdges(Level l, const parlay::sequence<E> &non_tree_edges);

  static auto RemoveDuplicates(parlay::sequence<V> &seq) -> parlay::sequence<V>;
  static auto RemoveDuplicates(parlay::sequence<V> &&seq)
      -> parlay::sequence<V>;
  inline void InsertIntoEdgeTable(const pair<V, V> &e, EType e_type,
                                  Level level);
  inline void DeleteFromEdgeTable(const pair<V, V> &e);
};

void BatchDynamicConnectivity::CheckRep() {
  // Rep invariant:

  // Basic checks
  // 1. `edges_` is a set of edges.
  // 2. `spanning_forests_[i].edges_` is a subset of `edges_`.
  // 3. `(v, non_tree_adjacency_lists_[i][v])` is a subset of `edges_`
  //    for all i and v.
  // 4. `Set( ...E(v, non_tree_adjacency_lists_[i][v]) for all i and v.
  //          ...spanning_forests_[j].edges_ for all j.)`
  //   = `edges_`.
  // 5. `spanning_forests_[i].edges_` is a subset of
  //    `spanning_forests_[j].edges_`    for all i < j.
  //     i.e. as we go down the levels. edges are pushed down from the
  //     higher levels to the lower levels. They are inserted into the top.
  // 6. all edges in `edges_`, `spanning_forests_[i].edges_` are distinct.

  // Component Size checks
  // 1. All components of edges at level i are of size <= 2^i.
  //    for all i in [0, max_level_].

  // MST checks
  // 1. M = {WeightedEdge( (e, w) | e in spanning_forests_[max_level_].edges_
  //                                 and  w = edges_[e].level}.
  //   is a minimum spanning forest of `edges_`.

  cout << "Checking rep invariant..." << endl;
  auto edges_seq = edges_.entries();

  // ------------------------------------------------------

  // Basic checks

  cout << "Basic checks..." << endl;
  set<E> edges_set;

  for (Level level = 0; level < max_level_; ++level) {
    // Check that `spanning_forests_[i].edges_` is a subset of `edges_`.
    auto spanning_forest_edges = parallel_spanning_forests_[level]->EdgesBothDirs_();
    assert(spanning_forest_edges.size() <= edges_seq.size());
    auto sps =
        std::count_if(edges_seq.begin(), edges_seq.end(), [&](const auto &e) {
          auto [edge, value] = e;
          auto [edge_level, e_type] = value;
          return edge_level <= level && e_type == EType::K_TREE;
        });
    assert(sps == static_cast<uint32_t>(spanning_forest_edges.size()));

    // insert all edges into `edges_set`.
    for (const auto &e : spanning_forest_edges) {
      edges_set.insert(e);
    }

    // Check that `(v, non_tree_adjacency_lists_[i][v])` is a subset of
    // `edges_`.
    auto non_tree_level_edges = non_tree_adjacency_lists_[level];
    for (V v = 0; v < num_vertices_; ++v) {
      auto opposite_vertices = non_tree_level_edges[v].entries();
      // construct a set of edges from the opposite vertices.
      auto opposite_edges = vector<pair<V, V>>(opposite_vertices.size());
      for (size_t i = 0; i < opposite_vertices.size(); ++i) {
        auto [opposite_vertex, _] = opposite_vertices[i];
        opposite_edges[i] = std::make_pair(v, opposite_vertex);
        edges_set.insert(opposite_edges[i]);
      }
      assert(opposite_edges.size() <= edges_seq.size());
      assert(
          std::count_if(edges_seq.begin(), edges_seq.end(), [&](const auto &e) {
            auto [edge, value] = e;
            auto [edge_level, e_type] = value;
            return edge_level == level && e_type == EType::K_NON_TREE && edge.first == v;
          }) == static_cast<uint32_t>(opposite_edges.size()));
    }
  }

  //   Check that `Set( ...E(v, non_tree_adjacency_lists_[i][v]),
  //            ...spanning_forests_[j].edges_)`
  //     = `edges_`.
  assert(edges_set.size() == edges_seq.size());
  set<E> edges_set_from_table;
  for (const auto &e : edges_seq) {
    auto [edge, _] = e;
    edges_set_from_table.insert(edge);
  }
  assert(edges_set == edges_set_from_table);

  // all edges in `edges_`, `spanning_forests_[i].edges_` are distinct.
  assert(edges_set.size() == edges_seq.size());
  assert(edges_set_from_table.size() == edges_seq.size());

  cout << "Basic checks passed." << endl;

  // ------------------------------------------------------

  // MST setup:

  cout << "MST checks ..." << endl;

  // we utilize a union find structure to keep track of the merges.
  constexpr auto kFind{elektra::find_variants::find_compress};
  auto unite{elektra::unite_variants::Unite{kFind}};

  // initialize parents_ to the identity map
  auto parents = parlay::sequence<uintE>::from_function(
      num_vertices_, [](uintE i) { return i; });

  // MST checks:

  // The algorithm to ensure our invariant is as follows:
  // - start the lowest level iterating to the top level.
  // - union every tree edge. This should succeed!
  // - we then try to union every non-tree edge but that should fail!
  //    this is because of the MST invariant which expresses that there
  //    shouldn't be non-tree edges connecting components. If they do, they
  //    should have been in the MST.
  for (Level level = 0; level < max_level_; level++) {
    sequence<tuple<E, EInfo>> level_tree_edges =
        parlay::filter(edges_seq, [&](auto &edge) {
          auto [e, info] = edge;
          return info.level == level && info.type == EType::K_TREE;
        });

    for (const auto &edge : level_tree_edges) {
      auto [e, info] = edge;
      auto [u, v] = e;
      if (u < v) {
        auto did_unite = unite(u, v, parents);
        assert(did_unite != UINT_E_MAX && "union of tree edges failed");
      }
    }

    sequence<tuple<E, EInfo>> level_non_tree_edges =
        parlay::filter(edges_seq, [&](auto &edge) {
          auto [e, info] = edge;
          return info.level == level && info.type == EType::K_NON_TREE;
        });

    for (const auto &edge : level_non_tree_edges) {
      auto [e, info] = edge;
      auto [u, v] = e;
      auto did_unite = unite(u, v, parents);
      assert(did_unite == UINT_E_MAX && "union of non-tree edges succeeded. "
                                        "This should not happen!");
    }
  }

  cout << "MST checks passed." << endl;

  // ------------------------------------------------------

  // component size checks

  // for (Level level = 0; level < max_level_; ++level) {
  //   // construct the components from the tree edges
  //   auto level_edges = parallel_spanning_forests_[level]->EdgesBothDirs_();

  //   auto components = vector<set<V>>(num_vertices_);

  //   for (const auto &e : level_edges) {
  //     auto [u, v] = e;

  //     for (const auto &component : components) {
  //       if (component.count(u) != 0U) {
  //         components[u].insert(v);
  //       } else if (component.count(v) != 0U) {
  //         components[v].insert(u);
  //       } else {
  //         auto new_component = set<V>();
  //         new_component.insert(u);
  //         new_component.insert(v);
  //         components.push_back(new_component);
  //       }
  //     }
  //   }

  //   // now insert the non-tree edges
  //   auto non_tree_edges = non_tree_adjacency_lists_[level];

  //   for (V v = 0; v < num_vertices_; ++v) {
  //     for (const auto &e : non_tree_edges[v].entries()) {
  //       auto [u, _] = e;
  //       for (const auto &component : components) {
  //         if (component.count(u) != 0U) {
  //           components[u].insert(v);
  //         } else if (component.count(v) != 0U) {
  //           components[v].insert(u);
  //         } else {
  //           auto new_component = set<V>();
  //           new_component.insert(u);
  //           new_component.insert(v);
  //           components.push_back(new_component);
  //         }
  //       }
  //     }
  //   }

  // TODO(sualeh): fix this
  // // assert that all components are disjoint
  // for (size_t i = 0; i < components.size(); ++i) {
  //   for (size_t j = i + 1; j < components.size(); ++j) {
  //     set<V> intersection;
  //     std::set_intersection(
  //         components[i].begin(), components[i].end(),
  //         components[j].begin(), components[j].end(),
  //         std::inserter(intersection, intersection.begin()));
  //     assert(intersection.empty() && "components are not disjoint");
  //   }
  // }

  // // assert that the size of each component is at most 2^(level)
  // for (const auto &component : components) {
  //   auto error = "component size error: component size = " +
  //                std::to_string(component.size()) +
  //                " for level = " + std::to_string(level) + "\n" +
  //                " max component size = " + std::to_string(1U << level) +
  //                "\n";
  //   assert(component.size() <= (1U << level) && error.c_str());
  // }

  // sum up the size of each component

  // auto size_sum = 0U;
  // for (const auto &component : components) {
  //   size_sum += component.size();
  // }
  // assert(size_sum == num_vertices_ && "component size sum error");
  // }

  cout << "leaving invariant checks." << endl;
}

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
  // note(tom): probably? make EdgeSet into a class where insert/find/delete
  // always makes the edge have e.first < e.second
  EInfo ei_rev = {level, e_type};
  edges_.insert(make_tuple(pair<V, V>(e.second, e.first), ei_rev));
}

inline void BatchDynamicConnectivity::DeleteFromEdgeTable(const pair<V, V> &e) {
  edges_.deleteVal(e);
  edges_.deleteVal(pair<V, V>(e.second, e.first));
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

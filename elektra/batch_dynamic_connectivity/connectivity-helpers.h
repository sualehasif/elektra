#pragma once

#include <parlay/delayed_sequence.h>
#include <parlay/hash_table.h>
#include <parlay/primitives.h>
#include <parlay/sequence.h>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <ostream>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "elektra_assert.h"
#include "hash_pair.h"
#include "hash_table_utils.h"
#include "macros.h"
#include "parallel_euler_tour_tree/hdt_euler_tour_tree.h"
#include "resizable_table.h"
#include "spanning_tree.h"
#include "sparse_table.h"
#include "union_find.h"
#include "utilities/sequence_utils.h"

#define VERTEX_LAYER_SIZE 50

/** Represents a vertex in a graph. */

using V = uint32_t;
constexpr V kV_Max = INT32_MAX - 1;

namespace bdcty {
using std::make_pair;
using std::make_tuple;
using std::make_unique;
using std::min;
using std::pair;
using std::tuple;
using std::unique_ptr;
using std::vector;

using parlay::sequence;
using uintE = unsigned int;

using E = pair<V, V>;
using EHash = HashIntPairStruct;
using BatchDynamicEtt = parallel_euler_tour_tree::HdtEulerTourTree;

using TreeSet = std::unordered_set<E, EHash>;

using Level = int;
constexpr Level kLevelMax = std::numeric_limits<V>::digits;  // log2(max V)

enum class EType {
  // Edge is in the spanning forest of the graph.
  K_NON_TREE,
  // Edge is not in the spanning forest of the graph.
  K_TREE,
};

// A struct that contains information about a particular edge.
struct EInfo {
  Level level;
  EType type;

  // The equality operator
  auto operator==(const EInfo &other) const -> bool {
    return level == other.level && type == other.type;
  }

  // The inequality operator
  auto operator!=(const EInfo &other) const -> bool {
    return !(*this == other);
  }
};

const bdcty::EInfo kEmptyInfo = {-1, bdcty::EType::K_NON_TREE};

// make a hashtable with an empty value type
// hash32 is sufficient
struct hash_kv {
  auto operator()(const V &k) const -> uint64_t { return parlay::hash64(k); }
};

using VertexConcurrentSet =
    elektra::sparse_table<V, elektra::empty, elektra::HashUint32Empty>;

using EdgeSet =
    elektra::resizable_table<std::pair<V, V>, bdcty::EInfo, HashIntPairStruct>;

/**
 * @brief The structure that contains the Non-Tree edges.
 * @details The structure that contains the Non-Tree edges for every level.
 *          The Non-Tree edges are the edges that are not in the spanning
 *          forest of the graph.
 *          The Non-Tree edges are stored as a hashtable for each vertex and
 *          level.
 */
class NonTreeAdjacencyList {
public:
  /** Constructor.
   *  @param[in] numVertices Number of vertices in the graph.
   */
  NonTreeAdjacencyList() = default;
  NonTreeAdjacencyList(V num_vertices, V max_level);

  /** Adds a batch of edges to the graph.
   *  The insertion happens in parallel and ensures that the internal
   hashtables have enough space.
   *  @param[in] u One endpoint of the edge.
   *  @param[in] v The other endpoint of the edge.
   */
  void BatchAddEdgesToLevel(parlay::sequence<pair<V, V>> edges,
                            bdcty::Level level);

  /** Returns the sequence of hashtables that contain the level-i non-tree
   * edges.
   *  @param[in] v The vertex.
   *  @param[in] level The level.
   *  @returns The vertices connected to vertex v by level-i non-tree edges.
   */
  [[nodiscard]] auto operator[](bdcty::Level level)
      -> parlay::sequence<VertexConcurrentSet> & {
    return non_tree_adjacency_lists_[level];
  }

  /** Returns the vertices connected to vertex v by level-i non-tree edges.
   *  @param[in] v The vertex.
   *  @param[in] level The level.
   *  @returns The vertices connected to vertex v by level-i non-tree edges.
   */
  [[maybe_unused]] [[nodiscard]] auto GetEdges(V v, bdcty::Level level) const
      -> const VertexConcurrentSet & {
    return non_tree_adjacency_lists_[level][v];
  }

  /** Returns the number of edges in the graph that are not in the spanning
   *  forest.
   *  @returns Number of edges in the graph that are not in the spanning
   *  forest.
   */
  //  int NumNonTreeEdges() const;

  /** Returns the number of edges in the graph that are not in the spanning
   *  forest and have a particular level.
   *  @param[in] level The level of the edges to be counted.
   *  @returns Number of edges in the graph that are not in the spanning
   *  forest and have a particular level.
   */
  //  int NumNonTreeEdgesForLevel(bdcty::Level level) const;

private:
  V num_vertices_{0};
  int max_level_{0};

  // `non_tree_adjacency_lists_[i][v]` contains the vertices connected to vertex
  // v by level-i non-tree edges.
  // `non_tree_adjacency_lists_` is a vector of size `max_level_`. (i.e. the
  // first index represents the levels.)
  // `non_tree_adjacency_lists_[i]` is a vector of size `num_vertices_`.
  // (i.e. the second index represents the vertices.)
  parlay::sequence<parlay::sequence<VertexConcurrentSet>>
      non_tree_adjacency_lists_;

  // generate a vertex layer
  static auto GenerateVertexLayer(V num_vertices);
};

auto NonTreeAdjacencyList::GenerateVertexLayer(V num_vertices) {
  auto vtx_layer = parlay::sequence<VertexConcurrentSet>::from_function(
      num_vertices, [&](auto n) {
        return VertexConcurrentSet(VERTEX_LAYER_SIZE);
      });

  return vtx_layer;
}

NonTreeAdjacencyList::NonTreeAdjacencyList(V num_vertices, V max_level)
    : num_vertices_(num_vertices), max_level_(max_level) {
  // We set up the non-tree adjacency lists.
  non_tree_adjacency_lists_ =
      parlay::sequence<parlay::sequence<VertexConcurrentSet>>(max_level);

  parlay::parallel_for(0, max_level, [&](auto i) {
    auto vtx_layer = GenerateVertexLayer(num_vertices);
    non_tree_adjacency_lists_[i] = vtx_layer;
  });
}

void NonTreeAdjacencyList::BatchAddEdgesToLevel(
    parlay::sequence<pair<V, V>> edges, bdcty::Level level) {
  // resize the edges to contain the reverse edges
  auto ne = edges.size();
  if (edges.capacity() < ne * 2) {
    edges.resize(ne * 2);
  }

  // add all the reverse e's of the edges to the graph
  parlay::parallel_for(0, ne, [&](auto i) {
    auto e = edges[i];
    edges[ne + i] = std::make_pair(e.second, e.first);
  });

#ifdef DEBUG
  std::cout << "---- Batch add edges to level ----" << std::endl;
  std::cout << "Adding " << ne * 2 << " = " << edges.size()
            << " edges to level: " << level << std::endl;
#endif

  // first we sort the edges by their first vertex
  parlay::integer_sort_inplace(
      edges, [](pair<V, V> &a) -> elektra::uintV { return a.first; });

#ifdef DEBUG
  // print the sorted edges
  std::cout << "Sorted edges:" << std::endl;
  for (auto i = 0; i < (int)ne * 2; i++) {
    std::cout << "(" << edges[i].first << ", " << edges[i].second << ")"
              << std::endl;
  }
#endif

  // get the appropriate level to insert the edges
  auto &level_lists = non_tree_adjacency_lists_[level];

  // find all the unique vertices
  auto unique_starting_vertex_indices =
      elektra::get_offsets(edges, [](const pair<V, V> &a, const pair<V, V> &b) {
        return a.first == b.first;
      });

#ifdef DEBUG
  std::cout << "Unique starting vertices:" << std::endl;
  for (const auto pos : unique_starting_vertex_indices) {
    auto &vtx = edges[pos];
    std::cout << "(" << vtx.first << ", " << vtx.second << ")@" << pos << ", ";
  }
  std::cout << std::endl;
#endif

  // then we add the edges to the graph
  parlay::parallel_for(0, unique_starting_vertex_indices.size(), [&](size_t i) {
    auto start_pos = unique_starting_vertex_indices[i];
    auto end_pos =
      i == unique_starting_vertex_indices.size() - 1
      ? ne * 2
      : unique_starting_vertex_indices[i + 1];
    auto &vtx = edges[start_pos];
    auto &vertex_list = level_lists[vtx.first];
    // get the starting and ending positions of the edges for this vertex
    auto num_edges_for_vertex = (end_pos - start_pos);
    // make sure the vertex list has enough space for the extra edges
    vertex_list.maybe_resize(num_edges_for_vertex);

// print the vertex list
#ifdef DEBUG
    std::cout << "V insertion for vertex " << vtx.first << ": ("
              << num_edges_for_vertex << " edges to be inserted)" << std::endl;
    for (V j = 0; j < num_edges_for_vertex; j++) {
      auto &edge = edges[start_pos + j];
      std::cout << "(" << edge.first << ", " << edge.second << "), ";
    }
    std::cout << std::endl;
    std::cout << "--------------------------------------------" << std::endl;
#endif

    // add the edges to the graph
    for (V j = 0; j < num_edges_for_vertex; j++) {
      auto &edge = edges[start_pos + j];
      vertex_list.insert(make_tuple(edge.second, elektra::empty{}));
    }
  });
}

template <class Edge, class EdgeSet>
sequence<Edge> RemoveUnknownEdges(const sequence<Edge> &se, const EdgeSet& edges) {
  return parlay::filter(se, [&](E e) { return edges.contains(e); });
}

template <typename T> auto PrintSequence(T &seq, const std::string &&name) {
  // TODO: using tuple_size make this print if not a tuple
  std::cout << name << ": ";
  for (size_t i = 0; i < seq.size(); i++) {
    std::cout << std::get<0>(seq[i]) << ", ";
  }
  std::cout << std::endl;
}

auto PrintEdgeSequence(sequence<std::tuple<E, elektra::empty>> &seq,
                       const std::string &&name) {
  std::cout << name << ": ";

  auto sorted_seq = parlay::sort(seq, [](const auto &a, const auto &b) {
    if (std::get<0>(a).first < std::get<0>(b).first) {
      return true;
    } else if (std::get<0>(a).first > std::get<0>(b).first) {
      return false;
    } else {
      return std::get<0>(a).second < std::get<0>(b).second;
    }
  });

  for (auto &i : seq) {
    auto &e = std::get<0>(i);
    std::cout << "(" << e.first << ", " << e.second << "), \n";
  }
  std::cout << std::endl;
}

auto PrintEdgeSequence(sequence<std::tuple<E, bdcty::EInfo>> &seq,
                       const std::string &&name) {
  std::cout << name << ": ";
  auto sorted_seq = parlay::sort(seq, [](const auto &a, const auto &b) {
    if (std::get<0>(a).first < std::get<0>(b).first) {
      return true;
    } else if (std::get<0>(a).first > std::get<0>(b).first) {
      return false;
    } else {
      return std::get<0>(a).second < std::get<0>(b).second;
    }
  });

  for (auto &i : sorted_seq) {
    auto &e = std::get<0>(i);
    std::cout << "(" << e.first << ", " << e.second << "), \n";
  }
  std::cout << std::endl;
}

template <typename T>
auto PrintEdgeSequence(T &edges, const std::string &&name) {
  std::cout << name << ": ";
  auto sorted_seq = parlay::sort(edges, [](const E &a, const E &b) {
    if (a.first < b.first) {
      return true;
    } else if (a.first > b.first) {
      return false;
    } else {
      return a.second < b.second;
    }
  });

  for (const auto &e : sorted_seq) {
    std::cout << "(" << e.first << ", " << e.second << "), \n";
  }
  std::cout << std::endl;
}

template <typename T>
auto PrintDifferenceSequence(sequence<std::tuple<E, elektra::empty>> &seq1,
                             sequence<T> &seq2, const std::string &&name) {
  std::cout << name << ": ";
  for (auto &i : seq1) {
    auto &e = std::get<0>(i);
    if (!seq2.contains(e)) {
      std::cout << "(" << e.first << ", " << e.second << "), ";
    }
  }
  std::cout << std::endl;
}

} // namespace bdcty

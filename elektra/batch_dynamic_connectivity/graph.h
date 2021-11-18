#pragma once

#include <parlay/hash_table.h>
#include <parlay/sequence.h>

#include <cstdint>
#include <ostream>
#include <tuple>
#include <utility>

#include "macros.h"
#include "resizable_table.h"

using std::make_tuple;
using std::pair;

#define VERTEX_LAYER_SIZE 50

/** Represents a vertex in a graph. */
typedef int Vertex;

/** Represents an edge in a directed graph. */
typedef std::pair<Vertex, Vertex> DirectedEdge;

namespace dynamicGraph {
/** Represents an edge in an undirected graph. */
struct UndirectedEdge {
  /** Constructor.
   *  @param[in] u One endpoint of the edge.
   *  @param[in] v The other endpoint of the edge.
   */
  UndirectedEdge(Vertex u, Vertex v) : first(u), second(v){};
  // UndirectedEdge() = delete;

  /** One endpoint of the edge. */
  Vertex first;
  /** The other endpoint of the edge. */
  Vertex second;
};
/** Defines how to print UndirectedEdge in an output stream. */
std::ostream& operator<<(std::ostream& out, const UndirectedEdge& edge);
/** Equality operator for UndirectedEdge. */
bool operator==(const UndirectedEdge& e1, const UndirectedEdge& e2) {
  return e1.first == e2.first && e1.second == e2.second;
};

/** For storing undirected edges in hash containers. For instance:
 *  \code
 *    std::unordered_map<UndirectedEdge, std::string, UndirectedEdgeHash>
 *       hash_map;
 *  \endcode
 */
struct UndirectedEdgeHash {
  /** Returns the hash value of an edge.
   *  @param[in] edge Edge to be hashed.
   *  @returns Hash value of the edge.
   */
  std::size_t operator()(const UndirectedEdge& edge) const {
    return std::hash<Vertex>()((int64_t)edge.first) ^
           std::hash<Vertex>()((int64_t)edge.second);
  };
};

}  // namespace dynamicGraph

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

  // The equality operator
  bool operator==(const EdgeInfo& other) const {
    return level == other.level && type == other.type;
  }

  // The inequality operator
  bool operator!=(const EdgeInfo& other) const { return !(*this == other); }
};

}  // namespace detail

namespace batchDynamicConnectivity {
// make a hashtable with an empty value type
using vertexConcurrentSet =
    elektra::resizable_table<Vertex, elektra::empty,
                             parlay::hash_numeric<Vertex>>;

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
  NonTreeAdjacencyList(int numVertices, int max_level_);

  /** Destructor. */
  ~NonTreeAdjacencyList();

  /** Adds a batch of edges to the graph.
   *  The insertion happens in parallel and ensures that the internal hashtables
   * have enough space.
   *  @param[in] u One endpoint of the edge.
   *  @param[in] v The other endpoint of the edge.
   */
  void BatchAddEdgesToLevel(parlay::sequence<pair<Vertex, Vertex>> edges,
                            detail::Level level);

  /** Returns the sequence of hashtables that contain the level-i non-tree
   * edges.
   *  @param[in] v The vertex.
   *  @param[in] level The level.
   *  @returns The vertices connected to vertex v by level-i non-tree edges.
   */
  parlay::sequence<vertexConcurrentSet>& operator[](detail::Level level) {
    return non_tree_adjacency_lists_[level];
  }

  /** Returns the vertices connected to vertex v by level-i non-tree edges.
   *  @param[in] v The vertex.
   *  @param[in] level The level.
   *  @returns The vertices connected to vertex v by level-i non-tree edges.
   */
  const vertexConcurrentSet& getEdges(Vertex v, detail::Level level) const {
    return non_tree_adjacency_lists_[level][v];
  }

  /** Returns the number of edges in the graph that are not in the spanning
   *  forest.
   *  @returns Number of edges in the graph that are not in the spanning
   *  forest.
   */
  int NumNonTreeEdges() const;

  /** Returns the number of edges in the graph that are not in the spanning
   *  forest and have a particular level.
   *  @param[in] level The level of the edges to be counted.
   *  @returns Number of edges in the graph that are not in the spanning
   *  forest and have a particular level.
   */
  int NumNonTreeEdgesForLevel(detail::Level level) const;

 private:
  int numVertices_;
  int max_level_;

  // `non_tree_adjacency_lists_[i][v]` contains the vertices connected to vertex
  // v by level-i non-tree edges.
  // `non_tree_adjacency_lists_` is a vector of size `max_level_`. (i.e. the
  // first index represents the levels.)
  // `non_tree_adjacency_lists_[i]` is a vector of size `num_vertices_`.
  // (i.e. the second index represents the vertices.)
  parlay::sequence<parlay::sequence<vertexConcurrentSet>>
      non_tree_adjacency_lists_;

  // generate a vertex layer
  auto generateInitialVertexLayer(int numVertices, int max_level_);
};

auto NonTreeAdjacencyList::generateInitialVertexLayer(int numVertices,
                                                      int max_level_) {
  auto vtxLayer = parlay::sequence<vertexConcurrentSet>::from_function(
      numVertices, [&](int n) {
        return vertexConcurrentSet(VERTEX_LAYER_SIZE,
                                   std::make_tuple(INT_E_MAX, elektra::empty{}),
                                   parlay::hash_numeric<Vertex>());
      });

  return vtxLayer;
}

NonTreeAdjacencyList::NonTreeAdjacencyList(int numVertices, int max_level_)
    : numVertices_(numVertices), max_level_(max_level_) {
  // We set up the non-tree adjacency lists.
  non_tree_adjacency_lists_ =
      parlay::sequence<parlay::sequence<vertexConcurrentSet>>(max_level_);

  parlay::parallel_for(0, max_level_, [&](int i) {
    auto vtxLayer = generateInitialVertexLayer(numVertices, max_level_);
    non_tree_adjacency_lists_[i] = vtxLayer;
  });
};

const vertexConcurrentSet& NonTreeAdjacencyList::getEdges(
    Vertex v, detail::Level level) const {
  return non_tree_adjacency_lists_[level][v];
}

void NonTreeAdjacencyList::BatchAddEdgesToLevel(
    parlay::sequence<pair<Vertex, Vertex>> edges, detail::Level level) {
  // resize the edges to contain the reverse edges
  auto ne = edges.size();
  edges.resize(ne * 2);
  // add all the reverse e's of the edges to the graph
  parlay::parallel_for(0, ne, [&](int i) {
    auto e = edges[i];
    edges[ne + i] = std::make_pair(e.second, e.first);
  });

  // first we sort the edges by their first vertex
  parlay::integer_sort_inplace(
      edges, [](const pair<Vertex, Vertex>& a, const pair<Vertex, Vertex>& b) {
        return a.first < b.first;
      });
  // get the appropriate level to insert the edges
  auto& levelLists = non_tree_adjacency_lists_[level];

  // find all the unique vertices
  auto uniqueStartingVertices = parlay::unique(
      edges, [](const pair<Vertex, Vertex>& a, const pair<Vertex, Vertex>& b) {
        return a.first == b.first;
      });

  int epos[uniqueStartingVertices.size()];

  // find the position of each vertex in the unique vertices
  parlay::parallel_for(0, uniqueStartingVertices.size(), [&](int i) {
    auto pos = parlay::find(edges, uniqueStartingVertices[i]);
    epos[i] = pos - edges.begin();
  });

  // then we add the edges to the graph
  parlay::parallel_for(0, uniqueStartingVertices.size(), [&](int i) {
    auto& vtx = uniqueStartingVertices[i];
    auto& vertexList = levelLists[vtx.first];

    // get the starting and ending positions of the edges for this vertex
    auto startPos = epos[i];
    auto endPos = (i == uniqueStartingVertices.size() - 1)
                      ? uniqueStartingVertices.size()
                      : epos[i + 1];
    auto numEdgesForVertex = endPos - startPos;

    // make sure the vertex list has enough space
    vertexList.update_nelms();
    vertexList.maybe_resize(numEdgesForVertex);

    // add the edges to the graph
    for (int j = 0; j < numEdgesForVertex; j++) {
      auto& edge = edges[startPos + j];
      vertexList.insert(make_tuple(edge.second, elektra::empty{}));
    }
  });
};

}  // namespace batchDynamicConnectivity

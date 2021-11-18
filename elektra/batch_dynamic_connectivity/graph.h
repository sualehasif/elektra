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

// hash32 is sufficient
struct hash_kv {
  uint64_t operator()(const Vertex& k) { return parlay::hash64(k); }
};

using vertexConcurrentSet =
    elektra::resizable_table<Vertex, elektra::empty, hash_kv>;

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
  NonTreeAdjacencyList();
  NonTreeAdjacencyList(int numVertices, int max_level_);

  /** Destructor. */
  // ~NonTreeAdjacencyList();

  /** Adds a batch of edges to the graph.
   *  The insertion happens in parallel and ensures that the internal
   hashtables
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
  auto generateVertexLayer(int numVertices, int max_level_);
};

auto NonTreeAdjacencyList::generateVertexLayer(int numVertices,
                                               int max_level_) {
  auto vtxLayer = parlay::sequence<vertexConcurrentSet>::from_function(
      numVertices, [&](int n) {
        return vertexConcurrentSet(
            VERTEX_LAYER_SIZE, std::make_tuple(INT_E_MAX, elektra::empty{}),
            std::make_tuple(INT_E_MAX - 1, elektra::empty{}), hash_kv());
      });

  return vtxLayer;
}

// default constructor
NonTreeAdjacencyList::NonTreeAdjacencyList() : numVertices_(0), max_level_(0) {}

NonTreeAdjacencyList::NonTreeAdjacencyList(int numVertices, int max_level_)
    : numVertices_(numVertices), max_level_(max_level_) {
  // We set up the non-tree adjacency lists.
  non_tree_adjacency_lists_ =
      parlay::sequence<parlay::sequence<vertexConcurrentSet>>(max_level_);

  parlay::parallel_for(0, max_level_, [&](int i) {
    auto vtxLayer = generateVertexLayer(numVertices, max_level_);
    non_tree_adjacency_lists_[i] = vtxLayer;
  });
};

void NonTreeAdjacencyList::BatchAddEdgesToLevel(
    parlay::sequence<pair<Vertex, Vertex>> edges, detail::Level level) {
  // resize the edges to contain the reverse edges
  auto ne = edges.size();
  if (edges.capacity() < ne * 2) {
    edges.resize(ne * 2);
  }

  // add all the reverse e's of the edges to the graph
  parlay::parallel_for(0, ne, [&](int i) {
    auto e = edges[i];
    edges[ne + i] = std::make_pair(e.second, e.first);
  });

#ifdef DEBUG
  std::cout << "Adding " << ne * 2 << " = " << edges.size() << " edges to level"
            << level << std::endl;
#endif

  // first we sort the edges by their first vertex
  parlay::integer_sort_inplace(
      edges, [](pair<Vertex, Vertex>& a) -> elektra::uintV { return a.first; });

#ifdef DEBUG
  // print the sorted edges
  std::cout << "Sorted edges:" << std::endl;
  for (int i = 0; i < ne * 2; i++) {
    std::cout << "(" << edges[i].first << ", " << edges[i].second << ")"
              << std::endl;
  }
#endif

  // get the appropriate level to insert the edges
  auto& levelLists = non_tree_adjacency_lists_[level];

  // find all the unique vertices
  auto uniqueStartingVertices = parlay::unique(
      edges, [](const pair<Vertex, Vertex>& a, const pair<Vertex, Vertex>& b) {
        return a.first == b.first;
      });

#ifdef DEBUG
  std::cout << "Unique starting vertices:" << std::endl;
  for (int i = 0; i < uniqueStartingVertices.size(); i++) {
    std::cout << "(" << uniqueStartingVertices[i].first << ", "
              << uniqueStartingVertices[i].second << "), ";
  }
  std::cout << std::endl;
#endif

  int epos[uniqueStartingVertices.size() + 1];

  // find the position of each vertex in the unique vertices
  parlay::parallel_for(0, uniqueStartingVertices.size(), [&](int i) {
    auto pos = parlay::find(edges, uniqueStartingVertices[i]);
    epos[i] = pos - edges.begin();
  });

  epos[uniqueStartingVertices.size()] = ne * 2;

  // print the positions of the vertices in DEBUG
#ifdef DEBUG
  std::cout << "Positions of unique starting vertices:" << std::endl;
  for (int i = 0; i < uniqueStartingVertices.size(); i++) {
    std::cout << epos[i] << ", ";
  }
  std::cout << std::endl;
#endif

  // then we add the edges to the graph
  parlay::parallel_for(0, uniqueStartingVertices.size(), [&](int i) {
    auto& vtx = uniqueStartingVertices[i];
    auto& vertexList = levelLists[vtx.first];

    // get the starting and ending positions of the edges for this vertex
    auto startPos = epos[i];
    auto endPos = epos[i + 1];
    auto numEdgesForVertex = (int)(endPos - startPos);

    // make sure the vertex list has enough space for the extra edges
    vertexList.maybe_resize(numEdgesForVertex);

// print the vertex list
#ifdef DEBUG
    std::cout << "Vertex insertion for vertex " << vtx.first << ": ("
              << numEdgesForVertex << " edges to be inserted)" << std::endl;
    for (int i = 0; i < numEdgesForVertex; i++) {
      auto& edge = edges[startPos + i];
      std::cout << "(" << edge.first << ", " << edge.second << "), ";
    }
    std::cout << std::endl;
#endif

    // add the edges to the graph
    for (int j = 0; j < numEdgesForVertex; j++) {
      auto& edge = edges[startPos + j];
      vertexList.insert(make_tuple(edge.second, elektra::empty{}));
    }
  });
};

}  // namespace batchDynamicConnectivity

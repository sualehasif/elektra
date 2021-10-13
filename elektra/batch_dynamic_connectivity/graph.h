#pragma once

#include <cstdint>
#include <ostream>
#include <utility>

#include "macros.h"

/** Represents a vertex in a graph. */
typedef int64_t Vertex;

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
  const Vertex first;
  /** The other endpoint of the edge. */
  const Vertex second;
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

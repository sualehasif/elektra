/** @file graph.hpp
 *  Some useful graph class declarations.
 *
 *  @author Tom Tseng (tomtseng)
 */
#pragma once

#include <cstdint>
#include <ostream>
#include <utility>

/** Represents a vertex in a graph. */
typedef int64_t Vertex;

/** Represents an edge in a directed graph. */
typedef std::pair<Vertex, Vertex> DirectedEdge;

/** Represents an edge in an undirected graph. */
struct UndirectedEdgeR {
  /** Constructor.
   *  @param[in] u One endpoint of the edge.
   *  @param[in] v The other endpoint of the edge.
   */
  UndirectedEdgeR(Vertex u, Vertex v);
  UndirectedEdgeR() = delete;

  /** One endpoint of the edge. */
  const Vertex first;
  /** The other endpoint of the edge. */
  const Vertex second;
};
/** Defines how to print UndirectedEdge in an output stream. */
std::ostream& operator<<(std::ostream& out, const UndirectedEdgeR& edge);
/** Equality operator for UndirectedEdge. */
bool operator==(const UndirectedEdgeR& e1, const UndirectedEdgeR& e2);

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
  std::size_t operator()(const UndirectedEdgeR& edge) const;
};

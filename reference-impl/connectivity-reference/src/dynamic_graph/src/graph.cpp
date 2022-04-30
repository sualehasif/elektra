#include <algorithm>
#include <dynamic_graph/graph.hpp>

#include "utils_hash.hpp"

UndirectedEdge::UndirectedEdge(Vertex u, Vertex v)
  : first{std::min(u, v)}
  , second{std::max(u, v)} {}

std::ostream& operator<<(std::ostream& out, const UndirectedEdge& edge) {
  out << '{' << edge.first << ", " << edge.second << '}';
  return out;
}

bool operator==(const UndirectedEdge& e1, const UndirectedEdge& e2) {
  return e1.first == e2.first && e1.second == e2.second;
}

std::size_t UndirectedEdgeHash::operator()(const UndirectedEdge& edge) const {
  return CombineHashes(Hash(edge.first), Hash(edge.second));
}

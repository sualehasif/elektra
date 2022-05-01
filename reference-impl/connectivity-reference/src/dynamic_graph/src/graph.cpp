#include <algorithm>
#include <dynamic_graph/graph.hpp>

#include "utils_hash.hpp"

UndirectedEdgeR::UndirectedEdgeR(Vertex u, Vertex v)
  : first{std::min(u, v)}
  , second{std::max(u, v)} {}

std::ostream& operator<<(std::ostream& out, const UndirectedEdgeR& edge) {
  out << '{' << edge.first << ", " << edge.second << '}';
  return out;
}

bool operator==(const UndirectedEdgeR& e1, const UndirectedEdgeR& e2) {
  return e1.first == e2.first && e1.second == e2.second;
}

std::size_t UndirectedEdgeHash::operator()(const UndirectedEdgeR& edge) const {
  return CombineHashes(Hash(edge.first), Hash(edge.second));
}

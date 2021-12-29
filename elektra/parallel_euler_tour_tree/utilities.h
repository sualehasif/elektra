#pragma once

#include <stdexcept>
#include <sstream>
#include <utility>

#include <parlay/alloc.h>

namespace parallel_euler_tour_tree {

// Integer type able to hold the number of vertices in the graph.
using v_int = uint32_t;

namespace _internal {

// For an undirected edge {u,v}, we allocate elements representing (u,v) and
// (v,u) contiguously.
template <typename Elem>
using EdgesAllocator = parlay::type_allocator<Elem[2]>;

// Given an element representing edge (u, v), returns (v, u) (assuming that that
// (u, v) and (v, u) are allocated contiguously).
template <typename Elem>
Elem* OppositeEdge(const Elem* uv) {
  const auto [u, v] = uv->id_;
  if (u < v) {
    return const_cast<Elem*>(uv + 1);
  } else if (u > v) {
    return const_cast<Elem*>(uv - 1);
  } else {
    std::ostringstream errorMessage;
    errorMessage << "expected edge, got ("
      << uv->id_.first << ' '
      << uv->id_.second;
    throw std::invalid_argument(errorMessage.str());
  }
}

// For an edge {u,v}, returns elements representing (u,v) and (v,u).
template <typename Elem>
std::pair<Elem*, Elem*> AllocEdges(v_int u, v_int v) {
  Elem* edges = *EdgesAllocator<Elem>::alloc();
  if (u < v) {
    return {edges, edges + 1};
  } else {
    return {edges + 1, edges};
  }
}

// For an edge {u,v}, destructs and frees the elements representing (u,v) and (v,u).
template <typename Elem>
void DestroyEdges(Elem* uv, Elem* vu) {
  uv->~Elem();
  vu->~Elem();
  if (uv < vu) {
    EdgesAllocator<Elem>::free(reinterpret_cast<Elem(*)[2]>(uv));
  } else {
    EdgesAllocator<Elem>::free(reinterpret_cast<Elem(*)[2]>(vu));
  }
}
// Given edge (u, v), destroys and frees both (u, v) and (v, u).
template <typename Elem>
void DestroyEdges(Elem* uv) {
  DestroyEdges(uv, OppositeEdge(uv));
}

}  // namespace _internal
}  // namespace parallel_euler_tour_tree

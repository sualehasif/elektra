#pragma once

#include <parlay/parallel.h>
#include <parlay/primitives.h>

#include "concurrentMap.h"
#include "element.h"
#include "hash_pair.h"
#include "utilities.h"

namespace parallel_euler_tour_tree {

namespace _internal {

// Used in Euler tour tree for mapping directed edges (pairs of ints) to the
// sequence element in the Euler tour representing the edge.
//
// Only one of (u, v) and (v, u) should be added to the map; we can find the
// other edge using OppositeEdge().
template <typename Elem>
class EdgeMap {
 public:
  EdgeMap() = delete;
  explicit EdgeMap(int num_vertices);
  ~EdgeMap();

  bool Insert(int u, int v, Elem* edge);
  bool Delete(int u, int v);
  Elem* Find(int u, int v);

  int IsEmpty() const;

  // Deallocate all elements held in the map. This assumes that all elements
  // in the map were allocated through `allocator`.
  void FreeElements();

  // Returns an estimate of the amount of memory in bytes this element occupies
  // beyond what's captured by sizeof(EdgeMap). This includes the amount of
  // memory occupied by the `Elem`s allocated and inserted as values into the
  // map.
  size_t AllocatedMemorySize() const;

 private:
  concurrent_map::concurrentHT<std::pair<int, int>, Elem*, HashIntPairStruct>
      map_;
};

template <typename Elem>
EdgeMap<Elem>::EdgeMap(int num_vertices)
    : map_{nullptr, static_cast<size_t>(num_vertices - 1),
           std::make_pair(-1, -1), std::make_pair(-2, -2)} {}

template <typename Elem>
EdgeMap<Elem>::~EdgeMap() { map_.del(); }

// Check whether the hash table is empty.
template <typename Elem>
int EdgeMap<Elem>::IsEmpty() const { return map_.n_elms == 0; }

template <typename Elem>
bool EdgeMap<Elem>::Insert(int u, int v, Elem* edge) {
  if (u > v) {
    std::swap(u, v);
    edge = OppositeEdge(edge);
  }
  return map_.insert(make_pair(u, v), edge);
}

template <typename Elem>
bool EdgeMap<Elem>::Delete(int u, int v) {
  if (u > v) {
    std::swap(u, v);
  }
  return map_.deleteVal(make_pair(u, v));
}

template <typename Elem>
Elem* EdgeMap<Elem>::Find(int u, int v) {
  if (u > v) {
    Elem* vu{*map_.find(make_pair(v, u))};
    return vu == nullptr ? nullptr : OppositeEdge(vu);
  } else {
    return *map_.find(make_pair(u, v));
  }
}

template <typename Elem>
void EdgeMap<Elem>::FreeElements() {
  parlay::parallel_for(0, map_.capacity, [&](size_t i) {
    auto kv{map_.table[i]};
    auto key{get<0>(kv)};
    if (key != map_.empty_key && key != map_.tombstone) {
      Elem* edge{get<1>(kv)};
      DestroyEdges(edge);
    }
  });
}

template <typename Elem>
size_t EdgeMap<Elem>::AllocatedMemorySize() const {
  return sizeof(map_.table[0]) * map_.capacity +
    parlay::reduce(
      parlay::map(map_.entries(), [&](std::tuple<std::pair<int, int>, Elem*> kv) {
        const Elem* elem{std::get<1>(kv)};
        return 2 * sizeof(*elem)
          + elem->AllocatedMemorySize() + OppositeEdge(elem)->AllocatedMemorySize();
      })
    );
}

}  // namespace _internal
}  // namespace parallel_euler_tour_tree

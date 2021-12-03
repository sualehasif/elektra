#pragma once

#include <parlay/alloc.h>
#include <parlay/parallel.h>
#include <parlay/primitives.h>

#include "concurrentMap.h"
#include "element.h"
#include "hash_pair.h"

namespace parallel_euler_tour_tree {

namespace _internal {

// The element allocator is used to allocate and deallocate elements.
// It is used to allocate the elements in the skip list.
// This is borrowed from parlay's type allocator.
using ElementAllocator = parlay::type_allocator<Element>;

// Used in Euler tour tree for mapping directed edges (pairs of ints) to the
// sequence element in the Euler tour representing the edge.
//
// Only one of (u, v) and (v, u) should be added to the map; we can find the
// other edge using the `twin_` pointer in `Element`.
class EdgeMap {
 public:
  EdgeMap() = delete;
  explicit EdgeMap(int num_vertices);
  ~EdgeMap();

  bool Insert(int u, int v, Element* edge);
  bool Delete(int u, int v);
  Element* Find(int u, int v);

  int IsEmpty() const;

  // Deallocate all elements held in the map. This assumes that all elements
  // in the map were allocated through `allocator`.
  void FreeElements();

  // Returns an estimate of the amount of memory in bytes this element occupies
  // beyond what's captured by sizeof(EdgeMap). This includes the amount of
  // memory occupied by the `Element`s allocated and inserted as values into the
  // map.
  size_t AllocatedMemorySize() const;

 private:
  concurrent_map::concurrentHT<std::pair<int, int>, Element*, HashIntPairStruct>
      map_;
};

EdgeMap::EdgeMap(int num_vertices)
    : map_{nullptr, static_cast<size_t>(num_vertices - 1),
           std::make_pair(-1, -1), std::make_pair(-2, -2)} {}

EdgeMap::~EdgeMap() { map_.del(); }

// Check whether the hash table is empty.
int EdgeMap::IsEmpty() const { return map_.n_elms == 0; }

bool EdgeMap::Insert(int u, int v, Element* edge) {
  if (u > v) {
    std::swap(u, v);
    edge = edge->twin_;
  }
  return map_.insert(make_pair(u, v), edge);
}

bool EdgeMap::Delete(int u, int v) {
  if (u > v) {
    std::swap(u, v);
  }
  return map_.deleteVal(make_pair(u, v));
}

Element* EdgeMap::Find(int u, int v) {
  if (u > v) {
    Element* vu{*map_.find(make_pair(v, u))};
    return vu == nullptr ? nullptr : vu->twin_;
  } else {
    return *map_.find(make_pair(u, v));
  }
}

// Assumes that all elements were allocator using ElementAllocator.

void EdgeMap::FreeElements() {
  parlay::parallel_for(0, map_.capacity, [&](size_t i) {
    auto kv{map_.table[i]};
    auto key{get<0>(kv)};
    if (key != map_.empty_key && key != map_.tombstone) {
      Element* element{get<1>(kv)};
      element->twin_->~Element();
      ElementAllocator::free(element->twin_);
      element->~Element();
      ElementAllocator::free(element);
    }
  });
}

size_t EdgeMap::AllocatedMemorySize() const {
  return sizeof(map_.table[0]) * map_.capacity +
    parlay::reduce(
      parlay::map(map_.entries(), [&](std::tuple<std::pair<int, int>, Element*> kv) {
        const Element* elem{std::get<1>(kv)};
        // Count both elem and elem->twin_ since only one appears in map_.
        return 2 * sizeof(*elem) + elem->AllocatedMemorySize() + elem->twin_->AllocatedMemorySize();
      })
    );
}

}  // namespace _internal
}  // namespace parallel_euler_tour_tree

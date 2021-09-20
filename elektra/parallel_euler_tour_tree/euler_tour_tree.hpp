#include <parlay/alloc.h>

#include <utility>

#include "elektra/concurrentMap.h"
#include "elektra/hash_pair.h"

// TODO: port over the skip list from the previous implementation
// TODO: port over the list allocator from the previous implementation
// TOOD:

namespace parallel_euler_tour_tree {

namespace _internal {

class Element : public parallel_skip_list::ElementBase<Element> {
 public:
  Element() : parallel_skip_list::ElementBase<Element>{} {}
  explicit Element(size_t random_int)
      : parallel_skip_list::ElementBase<Element>{random_int} {}

  // If this element represents edge (u, v), `twin` should point towards (v, u).
  Element *twin_{nullptr};
  // When batch splitting, we mark this as `true` for an edge that we will
  // splice out in the current round of recursion.
  bool split_mark_{false};

 private:
  friend class parallel_skip_list::ElementBase<Element>;
  static void DerivedInitialize() {}
  static void DerivedFinish() {}
};

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

  bool Insert(int u, int v, Element *edge);
  bool Delete(int u, int v);
  Element *Find(int u, int v);

  // Deallocate all elements held in the map. This assumes that all elements
  // in the map were allocated through `allocator`.
  void FreeElements(Element *el);

 private:
  concurrent_map::concurrentHT<std::pair<int, int>, Element *,
                               HashIntPairStruct>
      map_;
};

EdgeMap::EdgeMap(int num_vertices)
    : map_{nullptr, static_cast<size_t>(num_vertices - 1),
           std::make_pair(-1, -1), std::make_pair(-2, -2)} {}

EdgeMap::~EdgeMap() { map_.del(); }

bool EdgeMap::Insert(int u, int v, Element *edge) {
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

Element *EdgeMap::Find(int u, int v) {
  if (u > v) {
    Element *vu{*map_.find(make_pair(v, u))};
    return vu == nullptr ? nullptr : vu->twin_;
  } else {
    return *map_.find(make_pair(u, v));
  }
}

void EdgeMap::FreeElements(Element *el) { ElementAllocator::free(el); }

}  // namespace _internal

// Euler tour trees represent forests. We may add an edge using `Link`, remove
// an edge using `Cut`, and query whether two vertices are in the same tree
// using `IsConnected`. This implementation can also exploit parallelism when
// many edges are added at once through `BatchLink` or many edges are deleted at
// once through `BatchCut`.
class EulerTourTree {
 public:
  EulerTourTree() = delete;
  // Initializes n-vertex forest with no edges.
  explicit EulerTourTree(int num_vertices);
  ~EulerTourTree();
  EulerTourTree(const EulerTourTree &) = delete;
  EulerTourTree(EulerTourTree &&) = delete;
  EulerTourTree &operator=(const EulerTourTree &) = delete;
  EulerTourTree &operator=(EulerTourTree &&) = delete;

  // Returns true if `u` and `v` are in the same tree in the represented forest.
  bool IsConnected(int u, int v) const;
  // Adds edge {`u`, `v`} to forest. The addition of this edge must not create a
  // cycle in the graph.
  void Link(int u, int v);
  // Removes edge {`u`, `v`} from forest. The edge must be present in the
  // forest.
  void Cut(int u, int v);

  // Adds all edges in the `len`-length array `links` to the forest. Adding
  // these edges must not create cycles in the graph.
  void BatchLink(std::pair<int, int> *links, int len);
  // Removes all edges in the `len`-length array `cuts` from the forest. These
  // edges must be present in the forest and must be distinct.
  void BatchCut(std::pair<int, int> *cuts, int len);

 private:
  void BatchCutRecurse(std::pair<int, int> *cuts, int len, bool *ignored,
                       _internal::Element **join_targets,
                       _internal::Element **edge_elements);

  int num_vertices_;
  _internal::Element *vertices_;
  _internal::EdgeMap edges_;
  pbbs::random randomness_;
};

}  // namespace parallel_euler_tour_tree
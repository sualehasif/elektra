#pragma once

#include <parlay/alloc.h>
#include <parlay/delayed_sequence.h>
#include <parlay/parallel.h>
#include <parlay/sequence.h>

#include <unordered_set>
#include <utility>

#include "hash_pair.h"
#include "edge_map.h"
#include "element.h"

namespace parallel_euler_tour_tree {

using std::pair;

// Euler tour trees represent forests. We may add an edge using `Link`, remove
// an edge using `Cut`, and query whether two vertices are in the same tree
// using `IsConnected`. This implementation can also exploit parallelism when
// many edges are added at once through `BatchLink` or many edges are deleted at
// once through `BatchCut`.
//
// This class template takes a template argument `Elem` that should be a
// subclass of parallel_euler_tour_tree::ElementBase. The argument specifies the
// underlying skip list elements storing ETT data. To augment the ETT (e.g., to
// be able to quickly answer queries like "how many elements are in this
// tree?"), `Elem` itself needs to be augmented appropriately.
//
// note(tomtseng): Ideally I would make this easier to use by making the
// template argument an augmentation function rather than a full augmented skip
// list implementation, but that's too much work to be worth the effort for now.
// To make this work, I'll need to write a library of skip list operations that
// make it easy to do useful things with augmented values so that users don't
// need to know much about the skip list implementation.
template <typename Elem>
class EulerTourTreeBase {
 public:
  EulerTourTreeBase() = delete;
  // Initializes n-vertex forest with no edges.
  explicit EulerTourTreeBase(int num_vertices);
  ~EulerTourTreeBase();
  EulerTourTreeBase(const EulerTourTreeBase&) = delete;
  EulerTourTreeBase(EulerTourTreeBase&&) = delete;
  EulerTourTreeBase& operator=(const EulerTourTreeBase&) = delete;
  EulerTourTreeBase& operator=(EulerTourTreeBase&&) = delete;

  // Returns true if `u` and `v` are in the same tree in the represented forest.
  bool IsConnected(int u, int v) const;

  // Adds edge {`u`, `v`} to forest. The addition of this edge must not create a
  // cycle in the graph.
  void Link(int u, int v);
  // Removes edge {`u`, `v`} from forest. The edge must be present in the
  // forest.
  void Cut(int u, int v);

  // Returns the representative of the vertex with id `u`.
  int GetRepresentative(int u) const;

  // Adds all edges in the `len`-length array `links` to the forest. Adding
  // these edges must not create cycles in the graph.
  void BatchLink(const parlay::sequence<std::pair<int, int>>& links);
  // Removes all edges in the `len`-length array `cuts` from the forest. These
  // edges must be present in the forest and must be distinct.
  void BatchCut(const parlay::sequence<std::pair<int, int>>& cuts);

  // Returns if the forest is empty.
  bool IsEmpty() const;

  // Prints the tree to stdout.
  void Print();
  // Returns an estimate of the amount of system memory in bytes this data
  // structure occupies.
  size_t MemorySize() const;

 protected:
  using ElementAllocator = parlay::type_allocator<Elem>;
  int num_vertices_;
  parlay::sequence<Elem> vertices_;
  _internal::EdgeMap edges_;
  parlay::random randomness_;

  // This version of `Link` can be used if `Elem`s representing edges should be
  // constructed with the Elem(size_t random_int, std::pair<int, int>&& id) constructor.
  // The caller must instead perform construct `Elems` representing edges (u, v)
  // and (v, u).
  void Link(int u, int v, Elem* uv, Elem* vu);
  // This version of `BatchLink` can be used if `Elem`s representing edges should be
  // constructed with the Elem(size_t random_int, std::pair<int, int>&& id) constructor.
  // The caller must instead provide a function that will perform the
  // construction.
  //
  // `construct(const parlay::random& randomness, size_t i, Elem* uv, Elem* vu) -> void`
  // should use placement new to construct uv and vu. If random numbers are
  // needed, use `randomness.ith_rand(2 * i)` and `randomness.ith_rand(2 * i + 1)`.
  // `i` corresponds to an index into `links[]`.
  template <typename ElemConstructor>
  void BatchLink(const parlay::sequence<std::pair<int, int>>& links, ElemConstructor construct);

 private:
  template <typename ElemConstructor>
  void BatchLinkSequential(
      const parlay::sequence<std::pair<int, int>>& links,
      ElemConstructor construct_elements);
  void BatchCutSequential(const parlay::sequence<std::pair<int, int>>& cuts);
  void BatchCutRecurse(
      const parlay::sequence<std::pair<int, int>>& cuts,
      parlay::sequence<bool>& ignored,
      parlay::sequence<Elem*>& join_targets,
      parlay::sequence<Elem*>& edge_elements);
};

// Euler tour tree augmented to be able to fetch component sizes.
class EulerTourTree : public EulerTourTreeBase<parallel_euler_tour_tree::Element> {
 public:
  explicit EulerTourTree(int num_vertices);

  // Returns the number of vertices in the connected component of vertex `v`.
  size_t ComponentSize(int v) const;
  // Returns all vertices in the connected component of vertex `v`.
  parlay::sequence<int> ComponentVertices(int v) const;
  // Returns all edges in the connected component of vertex `v`.
  parlay::sequence<std::pair<int, int>> ComponentEdges(int v) const;

 private:
  using Base = EulerTourTreeBase<Element>;
};

// Some structural details:
// To represent a tree in a forest, for each vertex v, add a loop edge (v, v),
// and for each edge {u, v}, replace it with two directed edges (u, v) and (v,
// u). Because the transformed tree has equal indegree and outdegree for each
// vertex, it admits an Euler tour.
//
// Take any Euler tour of the tree and place it in a circular sequence data
// structure like a skip list. Euler tours behave nicely under edge additions
// and edge deletions, so links and cuts reduce to a few splits and joins on
// sequences.

namespace {

// On BatchCut, randomly ignore 1/`kBatchCutRecursiveFactor` cuts and recurse
// on them later.
constexpr int kBatchCutRecursiveFactor{100};

}  // namespace

template <typename E>
EulerTourTreeBase<E>::EulerTourTreeBase(int num_vertices)
    : num_vertices_{num_vertices}, edges_{num_vertices_}, randomness_{} {
  // TODO: figure out the appropriate new array function for this.
  // vertices_ = pbbs::new_array_no_init<Element>(num_vertices_);

  vertices_ =
      parlay::sequence<E,
                       parlay::internal::sequence_default_allocator<E>,
                       false>::uninitialized(num_vertices_);

  parlay::parallel_for(0, num_vertices_, [&](size_t i) {
    new (&vertices_[i]) E{randomness_.ith_rand(i), make_pair(i, i)};
    // The Euler tour on a vertex v (a singleton tree) is simply (v, v).
    E::JoinWithoutUpdate(&vertices_[i], &vertices_[i]);
  });

  randomness_ = randomness_.next();
}

template <typename E>
EulerTourTreeBase<E>::~EulerTourTreeBase() {
  edges_.FreeElements();
}

template <typename E>
bool EulerTourTreeBase<E>::IsConnected(int u, int v) const {
  return vertices_[u].FindRepresentative() == vertices_[v].FindRepresentative();
}

template <typename E>
int EulerTourTreeBase<E>::GetRepresentative(int u) const {
  return vertices_[u].FindRepresentativeVertex();
}

// Prints the tree to stdout.
template <typename E>
void EulerTourTreeBase<E>::Print() {
  // make a set of all edges
  std::unordered_set<std::pair<int, int>, HashIntPairStruct> edges;

  for (int i = 0; i < num_vertices_; i++) {
    for (const auto& edge : vertices_[i].GetEdges()) {
      edges.insert(edge);
    }
  }

  // print all the edges
  for (const auto& edge : edges) {
    std::cout << "(" << edge.first << ", " << edge.second << "), ";
  }
  std::cout << std::endl;
}

template <typename E>
size_t EulerTourTreeBase<E>::MemorySize() const {
  // Estimate size of vertices_.
  const auto vertexSizes{parlay::delayed_seq<size_t>(vertices_.size(), [&](const size_t i) {
    return vertices_[i].AllocatedMemorySize() + sizeof(vertices_[i]);
  })};

  return sizeof(*this) + parlay::reduce(vertexSizes) + edges_.AllocatedMemorySize();
}

// Checks if the tree is empty.
template <typename E>
bool EulerTourTreeBase<E>::IsEmpty() const { return edges_.IsEmpty(); }

template <typename E>
void EulerTourTreeBase<E>::Link(int u, int v, E* uv, E* vu) {
  uv->twin_ = vu;
  vu->twin_ = uv;
  edges_.Insert(u, v, uv);
  E* u_left{&vertices_[u]};
  E* v_left{&vertices_[v]};
  E* u_right{u_left->SplitWithoutUpdate()};
  E* v_right{v_left->SplitWithoutUpdate()};
  E::JoinWithoutUpdate(u_left, uv);
  E::JoinWithoutUpdate(uv, v_right);
  E::JoinWithoutUpdate(v_left, vu);
  E::JoinWithoutUpdate(vu, u_right);
  E::BatchUpdate(parlay::sequence<E*>{{u_left, uv, v_left, vu}});
}

template <typename E>
void EulerTourTreeBase<E>::Link(int u, int v) {
  E* uv = ElementAllocator::alloc();
  E* vu = ElementAllocator::alloc();
  new (uv) E{randomness_.ith_rand(0), make_pair(u, v)};
  new (vu) E{randomness_.ith_rand(1), make_pair(v, u)};
  randomness_ = randomness_.next();
  Link(u, v, uv, vu);
}

template <typename E>
template <typename ElemConstructor>
void EulerTourTreeBase<E>::BatchLinkSequential(
    const parlay::sequence<std::pair<int, int>>& links,
    ElemConstructor construct_elements) {
  // TODO(tomtseng): We should do all links without doing any augmented value
  // updates, then do all augmented value updates at the end in a single
  // BatchUpdate.
  for (size_t i = 0; i < links.size(); i++) {
    E* uv = ElementAllocator::alloc();
    E* vu = ElementAllocator::alloc();
    construct_elements(randomness_, i, uv, vu);
    Link(links[i].first, links[i].second, uv, vu);
  }
  randomness_ = randomness_.next();
}

template <typename E>
template <typename ElemConstructor>
void EulerTourTreeBase<E>::BatchLink(
    const parlay::sequence<std::pair<int, int>>& links,
    ElemConstructor construct_elements) {
  const size_t len = links.size();
  if (len <= 75) {
    return BatchLinkSequential(links, construct_elements);
  }
  // For each added edge {x, y}, allocate elements (x, y) and (y, x).
  // For each vertex x that shows up in an added edge, split on (x, x). Let
  // succ(x) denote the successor of (x, x) prior to splitting.
  // For each vertex x, identify which y_1, y_2, ... y_k that x will be newly
  // connected to by performing a semisort on {(x, y), (y, x) : {x, y} is an
  // added edge}.
  // If x has new neighbors y_1, y_2, ..., y_k, join (x, x) to (x, y_1). Join
  // (y_i,x) to (x, y_{i+1}) for each i < k. Join (y_k, x) to succ(x).

  // pair<int, int>* links_both_dirs{pbbs::new_array_no_init<pair<int, int>>
  //                                                    (2 *len)};

  auto links_both_dirs =
      parlay::sequence<pair<int, int>>::uninitialized(2 * len);

  // parallel_for(int i = 0; i < len; i++) {
  //   links_both_dirs[2 * i] = links[i];
  //   links_both_dirs[2 * i + 1] = make_pair(links[i].second, links[i].first);
  // }

  parlay::parallel_for(0, len, [&](size_t i) {
    const int u{links[i].first};
    const int v{links[i].second};
    links_both_dirs[2 * i] = links[i];
    links_both_dirs[2 * i + 1] = {v, u};

    E* uv{ElementAllocator::alloc()};
    E* vu{ElementAllocator::alloc()};
    construct_elements(randomness_, i, uv, vu);
    uv->twin_ = vu;
    vu->twin_ = uv;
    edges_.Insert(u, v, uv);
  });
  randomness_ = randomness_.next();

  // intSort::iSort(links_both_dirs, 2 * len, num_vertices_ + 1,
  //                firstF<int, int>());
  // parlay::integer_sort(links_both_dirs, firstF<int, int>());
  auto getFirst = [](std::pair<int, int> a) { return (uint)a.first; };
  parlay::integer_sort_inplace(links_both_dirs, getFirst);

  // Element** split_successors{pbbs::new_array_no_init<Element*>(2 * len)};
  auto split_successors = parlay::sequence<E*>::uninitialized(2 * len);

  parlay::parallel_for(0, 2 * len, [&](size_t i) {
    int u, v;
    std::tie(u, v) = links_both_dirs[i];

    // split on each vertex that appears in the input
    if (i == 2 * len - 1 || u != links_both_dirs[i + 1].first) {
      split_successors[i] = vertices_[u].SplitWithoutUpdate();
    }
  });

  auto update_targets = parlay::sequence<E*>::uninitialized(2 * len);
  parlay::parallel_for(0, 2 * len, [&](size_t i) {
    int u, v;
    std::tie(u, v) = links_both_dirs[i];
    E* uv{edges_.Find(u, v)};
    E* vu{uv->twin_};
    if (i == 0 || u != links_both_dirs[i - 1].first) {
      update_targets[i] = &vertices_[u];
      E::JoinWithoutUpdate(&vertices_[u], uv);
    }
    if (i == 2 * len - 1 || u != links_both_dirs[i + 1].first) {
      update_targets[i] = vu;
      E::JoinWithoutUpdate(vu, split_successors[i]);
    } else {
      int u2, v2;
      std::tie(u2, v2) = links_both_dirs[i + 1];
      update_targets[i] = vu;
      E::JoinWithoutUpdate(vu, edges_.Find(u2, v2));
    }
  });
  // We only need to update the elements that participated in joins since the elements that
  // participated in splits participated in joins as well.
  E::BatchUpdate(update_targets);
}

template <typename E>
void EulerTourTreeBase<E>::BatchLink(const parlay::sequence<std::pair<int, int>>& links) {
  const auto construct{[&](const parlay::random& randomness, size_t i, E* uv, E* vu) {
    const int u{links[i].first};
    const int v{links[i].second};
    new (uv) E{randomness.ith_rand(2 * i), make_pair(u, v)};
    new (vu) E{randomness.ith_rand(2 * i + 1), make_pair(v, u)};
  }};
  BatchLink(links, construct);
}

template <typename E>
void EulerTourTreeBase<E>::Cut(int u, int v) {
  E* uv{edges_.Find(u, v)};
  E* vu{uv->twin_};
  edges_.Delete(u, v);
  E* u_left{uv->GetPreviousElement()};
  E* v_left{vu->GetPreviousElement()};
  E* v_right{uv->SplitWithoutUpdate()};
  E* u_right{vu->SplitWithoutUpdate()};
  u_left->SplitWithoutUpdate();
  v_left->SplitWithoutUpdate();

  uv->~E();
  ElementAllocator::free(uv);
  vu->~E();
  ElementAllocator::free(vu);

  E::JoinWithoutUpdate(u_left, u_right);
  E::JoinWithoutUpdate(v_left, v_right);
  E::BatchUpdate(parlay::sequence<E*>{{u_left, v_left}});
}

template <typename E>
void EulerTourTreeBase<E>::BatchCutSequential(const parlay::sequence<std::pair<int, int>>& cuts) {
  // TODO(tomtseng): We should do all cuts without doing any augmented value
  // updates, then do all augmented value updates at the end in a single
  // BatchUpdate. Or is this too difficult to do correctly due to edge elements
  // being deleted and returned to the memory allocator?
  for (const auto& cut : cuts) {
    Cut(cut.first, cut.second);
  }
}

// `ignored`, `join_targets`, and `edge_elements` are scratch space.
// `ignored[i]` will be set to true if `cuts[i]` will not be executed in this
// round of recursion.
// `join_targets` will store sequence elements that need to be joined to each other.
// `edge_elements[i]` will store a pointer to the sequence element corresponding to
// edge `cuts[i]`.
template <typename E>
void EulerTourTreeBase<E>::BatchCutRecurse(const parlay::sequence<std::pair<int, int>>& cuts,
                                    parlay::sequence<bool>& ignored,
                                    parlay::sequence<E*>& join_targets,
                                    parlay::sequence<E*>& edge_elements) {
  const size_t len = cuts.size();
  if (len <= 75) {
    return BatchCutSequential(cuts);
  }

  // Notation: "(x, y).next" is the next element in the tour (x, y) is in. "(x,
  // y).prev" is the previous element. "(x, y).twin" is (y, x).
  // For each edge {x, y} to cut:
  // Sequentially, we'd want to join (y, x).prev to (x, y).next and (x, y).prev
  // to (y, x).next. We can't correctly do this if any of those four elements
  // are to be cut and removed as well. Instead, for dealing with connecting (y,
  // x).prev to (x, y).next (dealing with connecting (x, y).prev to (y,x).next
  // is symmetric), we do the following:
  // - If (y, x).prev is to be cut, then do nothing --- some other thread will
  // deal with this.
  // - Otherwise, start with element e = (x, y).next. So long as e is to be cut,
  // traverse to the next possible join location at e := e.next.twin. Join
  // (y,x).prev to e.
  // This strategy doesn't have good depth since we may have to traverse on e
  // for a long time. To fix this, we randomly ignore some cuts so that all
  // traversal lengths are O(log n) with high probability. We perform all
  // unignored cuts as described above, and recurse on the ignored cuts
  // afterwards.

  parlay::parallel_for(0, len, [&](size_t i) {
    ignored[i] = randomness_.ith_rand(i) % kBatchCutRecursiveFactor == 0;

    if (!ignored[i]) {
      int u, v;
      std::tie(u, v) = cuts[i];
      E* uv{edges_.Find(u, v)};
      edge_elements[i] = uv;
      E* vu{uv->twin_};
      uv->split_mark_ = vu->split_mark_ = true;
    }
  });

  randomness_ = randomness_.next();

  parlay::parallel_for(0, len, [&](size_t i) {
    if (ignored[i]) {
      join_targets[4 * i] = join_targets[4 * i + 2] = nullptr;
      return;
    }
    E* uv{edge_elements[i]};
    E* vu{uv->twin_};

    E* left_target{uv->GetPreviousElement()};
    if (left_target->split_mark_) {
      join_targets[4 * i] = nullptr;
    } else {
      E* right_target{vu->GetNextElement()};
      while (right_target->split_mark_) {
        right_target = right_target->twin_->GetNextElement();
      }
      join_targets[4 * i] = left_target;
      join_targets[4 * i + 1] = right_target;
    }

    left_target = vu->GetPreviousElement();
    if (left_target->split_mark_) {
      join_targets[4 * i + 2] = nullptr;
    } else {
      E* right_target{uv->GetNextElement()};
      while (right_target->split_mark_) {
        right_target = right_target->twin_->GetNextElement();
      }
      join_targets[4 * i + 2] = left_target;
      join_targets[4 * i + 3] = right_target;
    }
  });

  parlay::parallel_for(0, len, [&](size_t i) {
    if (!ignored[i]) {
      E* uv{edge_elements[i]};
      E* vu{uv->twin_};
      uv->SplitWithoutUpdate();
      vu->SplitWithoutUpdate();
      E* predecessor{uv->GetPreviousElement()};
      if (predecessor != nullptr) {
        predecessor->SplitWithoutUpdate();
      }
      predecessor = vu->GetPreviousElement();
      if (predecessor != nullptr) {
        predecessor->SplitWithoutUpdate();
      }
    }
  });

  parlay::parallel_for(0, len, [&](size_t i) {
    if (!ignored[i]) {
      // Here we must use `edge_elements[i]` instead of `edges_.Find(u, v)`
      // because the concurrent hash table cannot handle simultaneous lookups
      // and deletions.
      E* uv{edge_elements[i]};
      E* vu{uv->twin_};

      uv->~E();
      ElementAllocator::free(uv);
      vu->~E();
      ElementAllocator::free(vu);

      int u, v;
      std::tie(u, v) = cuts[i];
      edges_.Delete(u, v);

      if (join_targets[4 * i] != nullptr) {
        E::JoinWithoutUpdate(join_targets[4 * i], join_targets[4 * i + 1]);
      }
      if (join_targets[4 * i + 2] != nullptr) {
        E::JoinWithoutUpdate(join_targets[4 * i + 2], join_targets[4 * i + 3]);
      }
    }
  });

  // We only need to update the elements that participated in joins since the elements that
  // participated in splits were either deleted or participated in joins.
  const auto update_targets{parlay::delayed_seq<E*>(2 * len, [&join_targets](const size_t i) {
    return join_targets[2 * i];
  })};
  E::BatchUpdate(update_targets);

  // parlay::sequence<std::pair<int, int>> full_cuts = *cuts;

  // seq::sequence<pair<int, int>> cuts_seq{seq::sequence<pair<int, int>>(cuts,
  // len)};
  // auto cuts_seq = parlay::sequence<pair<int, int>>(cuts);

  // seq::sequence<bool> ignored_seq{seq::sequence<bool>(ignored, len)};
  // auto ignored_seq = parlay::sequence<bool>(ignored);

  // seq::sequence<pair<int, int>> next_cuts_seq{
  //     pbbs::pack(cuts_seq, ignored_seq)};
  auto next_cuts_seq =
      parlay::sequence<std::pair<int, int>>(parlay::pack(cuts, ignored));

  // parlay::sequence<std::pair<int, int>>* next_cuts_seq_ptr = &next_cuts_seq;

  BatchCutRecurse(next_cuts_seq, ignored, join_targets, edge_elements);
  // pbbs::delete_array(next_cuts_seq.as_array(), next_cuts_seq.size());
}

template <typename E>
void EulerTourTreeBase<E>::BatchCut(const parlay::sequence<std::pair<int, int>>& cuts) {
  const size_t len = cuts.size();
  if (len <= 75) {
    return BatchCutSequential(cuts);
  }

  auto ignored = parlay::sequence<bool>::uninitialized(len);
  auto join_targets = parlay::sequence<E*>::uninitialized(4 * len);
  auto edge_elements = parlay::sequence<E*>::uninitialized(len);
  BatchCutRecurse(cuts, ignored, join_targets, edge_elements);
}

EulerTourTree::EulerTourTree(int num_vertices) : Base{num_vertices} {}

parlay::sequence<int> EulerTourTree::ComponentVertices(int v) const {
  // TODO(tomtseng) this is a bad, sequential way to get vertices of a
  // component, i'm just not optimizing this right now because we're not using
  // this anywhere important
  const auto edges{ComponentEdges(v)};
  auto vertices{parlay::sequence<int>::uninitialized(edges.size() + 1)};
  if (edges.empty()) {
    vertices[0] = v;
  } else {
    size_t curr_idx{0};
    std::unordered_set<int> seen_vertices;
    for (const auto& edge : edges) {
      for (const auto u : {edge.first, edge.second}) {
        if (seen_vertices.count(u) > 0) {
          continue;
        }
        seen_vertices.insert(u);
        vertices[curr_idx++] = u;
      }
    }
  }
  return vertices;
}

parlay::sequence<std::pair<int, int>> EulerTourTree::ComponentEdges(int v) const {
  return vertices_[v].GetEdges();
}

size_t EulerTourTree::ComponentSize(int v) const {
  return vertices_[v].GetComponentSize();
}

}  // namespace parallel_euler_tour_tree

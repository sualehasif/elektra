#pragma once

#include <parlay/alloc.h>
#include <parlay/delayed_sequence.h>
#include <parlay/parallel.h>
#include <parlay/sequence.h>

#include <unordered_set>
#include <utility>

#include "hash_pair.h"
#include "edge_map.hpp"
#include "element.hpp"

namespace parallel_euler_tour_tree {

using Element = _internal::Element;
using ElementAllocator = parlay::type_allocator<Element>;

using std::pair;

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
  EulerTourTree(const EulerTourTree&) = delete;
  EulerTourTree(EulerTourTree&&) = delete;
  EulerTourTree& operator=(const EulerTourTree&) = delete;
  EulerTourTree& operator=(EulerTourTree&&) = delete;

  // Returns true if `u` and `v` are in the same tree in the represented forest.
  bool IsConnected(int u, int v) const;
  // Returns all vertices in the connected component of vertex `v`.
  parlay::sequence<int> ComponentVertices(int v);
  // Returns all edges in the connected component of vertex `v`.
  parlay::sequence<std::pair<int, int>> ComponentEdges(int v);
  // Returns the number of vertices in the connected component of vertex `v`.
  size_t ComponentSize(int v);

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
  void BatchLink(parlay::sequence<std::pair<int, int>>& links);
  // Removes all edges in the `len`-length array `cuts` from the forest. These
  // edges must be present in the forest and must be distinct.
  void BatchCut(parlay::sequence<std::pair<int, int>>& cuts);

  // Returns if the forest is empty.
  bool IsEmpty() const;

  // Prints the tree to stdout.
  void Print();
  // Returns an estimate of the amount of system memory in bytes this data
  // structure occupies.
  size_t MemorySize() const;

 private:
  void BatchCutRecurse(
      const parlay::sequence<std::pair<int, int>>& cuts,
      parlay::sequence<bool>& ignored,
      parlay::sequence<Element*>& join_targets,
      parlay::sequence<Element*>& edge_elements);

  int num_vertices_;
  parlay::sequence<parallel_euler_tour_tree::Element> vertices_;
  _internal::EdgeMap edges_;
  parlay::random randomness_;
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

void BatchCutSequential(EulerTourTree* ett,
                        const parlay::sequence<std::pair<int, int>>& cuts) {
  // TODO(tomtseng): We should do all cuts without doing any augmented value
  // updates, then do all augmented value updates at the end in a single
  // BatchUpdate. Or is this too difficult to do correctly due to edge elements
  // being deleted and returned to the memory allocator?
  for (const auto& cut : cuts) {
    ett->Cut(cut.first, cut.second);
  }
}

void BatchLinkSequential(EulerTourTree* ett,
                         const parlay::sequence<std::pair<int, int>>& links) {
  // TODO(tomtseng): We should do all links without doing any augmented value
  // updates, then do all augmented value updates at the end in a single
  // BatchUpdate.
  for (const auto& link : links) {
    ett->Link(link.first, link.second);
  }
}

}  // namespace

EulerTourTree::EulerTourTree(int num_vertices)
    : num_vertices_{num_vertices}, edges_{num_vertices_}, randomness_{} {
  // TODO: figure out the appropriate new array function for this.
  // vertices_ = pbbs::new_array_no_init<Element>(num_vertices_);

  vertices_ =
      parlay::sequence<Element,
                       parlay::internal::sequence_default_allocator<Element>,
                       false>::uninitialized(num_vertices_);

  parlay::parallel_for(0, num_vertices_, [&](size_t i) {
    new (&vertices_[i]) Element{randomness_.ith_rand(i), make_pair(i, i)};
    // The Euler tour on a vertex v (a singleton tree) is simply (v, v).
    Element::JoinWithoutUpdate(&vertices_[i], &vertices_[i]);
  });

  randomness_ = randomness_.next();
}

EulerTourTree::~EulerTourTree() {
  edges_.FreeElements();
}

bool EulerTourTree::IsConnected(int u, int v) const {
  return vertices_[u].FindRepresentative() == vertices_[v].FindRepresentative();
}

int EulerTourTree::GetRepresentative(int u) const {
  return vertices_[u].FindRepresentativeVertex();
}

parlay::sequence<int> EulerTourTree::ComponentVertices(int v) {
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

parlay::sequence<std::pair<int, int>> EulerTourTree::ComponentEdges(int v) {
  return vertices_[v].GetEdges();
}

size_t EulerTourTree::ComponentSize(int v) {
  return vertices_[v].GetComponentSize();
}

// Prints the tree to stdout.
void EulerTourTree::Print() {
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

size_t EulerTourTree::MemorySize() const {
  // Estimate size of vertices_.
  const auto vertexSizes{parlay::delayed_seq<size_t>(vertices_.size(), [&](const size_t i) {
    return vertices_[i].AllocatedMemorySize() + sizeof(vertices_[i]);
  })};

  return sizeof(*this) + parlay::reduce(vertexSizes) + edges_.AllocatedMemorySize();
}

// Checks if the tree is empty.
bool EulerTourTree::IsEmpty() const { return edges_.IsEmpty(); }

void EulerTourTree::Link(int u, int v) {
  Element* uv = ElementAllocator::alloc();
  new (uv) Element{randomness_.ith_rand(0), make_pair(u, v)};
  Element* vu = ElementAllocator::alloc();
  new (vu) Element{randomness_.ith_rand(1), make_pair(v, u)};
  randomness_ = randomness_.next();
  uv->twin_ = vu;
  vu->twin_ = uv;
  edges_.Insert(u, v, uv);
  Element* u_left{&vertices_[u]};
  Element* v_left{&vertices_[v]};
  Element* u_right{u_left->SplitWithoutUpdate()};
  Element* v_right{v_left->SplitWithoutUpdate()};
  Element::JoinWithoutUpdate(u_left, uv);
  Element::JoinWithoutUpdate(uv, v_right);
  Element::JoinWithoutUpdate(v_left, vu);
  Element::JoinWithoutUpdate(vu, u_right);
  Element::BatchUpdate(parlay::sequence<Element*>{{u_left, uv, v_left, vu}});
}

template <class E1, class E2>
struct firstF {
  E1 operator()(std::pair<E1, E2> a) { return a.first; }
};

template <class E1, class E2>
struct secondF {
  E2 operator()(std::pair<E1, E2> a) { return a.second; }
};

void EulerTourTree::BatchLink(parlay::sequence<std::pair<int, int>>& links) {
  const size_t len = links.size();
  if (len <= 75) {
    BatchLinkSequential(this, links);
    return;
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
    links_both_dirs[2 * i] = links[i];
    links_both_dirs[2 * i + 1] = {links[i].second, links[i].first};
  });

  // intSort::iSort(links_both_dirs, 2 * len, num_vertices_ + 1,
  //                firstF<int, int>());
  // parlay::integer_sort(links_both_dirs, firstF<int, int>());
  auto getFirst = [](std::pair<int, int> a) { return (uint)a.first; };
  parlay::integer_sort_inplace(links_both_dirs, getFirst);

  // Element** split_successors{pbbs::new_array_no_init<Element*>(2 * len)};
  auto split_successors = parlay::sequence<Element*>::uninitialized(2 * len);

  parlay::parallel_for(0, 2 * len, [&](size_t i) {
    int u, v;
    std::tie(u, v) = links_both_dirs[i];

    // split on each vertex that appears in the input
    if (i == 2 * len - 1 || u != links_both_dirs[i + 1].first) {
      split_successors[i] = vertices_[u].SplitWithoutUpdate();
    }

    // allocate edge element
    if (u < v) {
      Element* uv = ElementAllocator::alloc();
      // new (uv) Element{randomness_.ith_rand(2 * i)};
      new (uv) Element{randomness_.ith_rand(2 * i), make_pair(u, v)};
      Element* vu = ElementAllocator::alloc();
      // new (vu) Element{randomness_.ith_rand(2 * i + 1)};
      new (vu) Element{randomness_.ith_rand(2 * i + 1), make_pair(v, u)};
      uv->twin_ = vu;
      vu->twin_ = uv;
      edges_.Insert(u, v, uv);
    }
  });

  randomness_ = randomness_.next();

  auto update_targets = parlay::sequence<Element*>::uninitialized(2 * len);
  parlay::parallel_for(0, 2 * len, [&](size_t i) {
    int u, v;
    std::tie(u, v) = links_both_dirs[i];
    Element* uv{edges_.Find(u, v)};
    Element* vu{uv->twin_};
    if (i == 0 || u != links_both_dirs[i - 1].first) {
      update_targets[i] = &vertices_[u];
      Element::JoinWithoutUpdate(&vertices_[u], uv);
    }
    if (i == 2 * len - 1 || u != links_both_dirs[i + 1].first) {
      update_targets[i] = vu;
      Element::JoinWithoutUpdate(vu, split_successors[i]);
    } else {
      int u2, v2;
      std::tie(u2, v2) = links_both_dirs[i + 1];
      update_targets[i] = vu;
      Element::JoinWithoutUpdate(vu, edges_.Find(u2, v2));
    }
  });
  // We only need to update the elements that participated in joins since the elements that
  // participated in splits participated in joins as well.
  Element::BatchUpdate(update_targets);
}

void EulerTourTree::Cut(int u, int v) {
  Element* uv{edges_.Find(u, v)};
  Element* vu{uv->twin_};
  edges_.Delete(u, v);
  Element* u_left{uv->GetPreviousElement()};
  Element* v_left{vu->GetPreviousElement()};
  Element* v_right{uv->SplitWithoutUpdate()};
  Element* u_right{vu->SplitWithoutUpdate()};
  u_left->SplitWithoutUpdate();
  v_left->SplitWithoutUpdate();

  uv->~Element();
  ElementAllocator::free(uv);
  vu->~Element();
  ElementAllocator::free(vu);

  Element::JoinWithoutUpdate(u_left, u_right);
  Element::JoinWithoutUpdate(v_left, v_right);
  Element::BatchUpdate(parlay::sequence<Element*>{{u_left, v_left}});
}

// `ignored`, `join_targets`, and `edge_elements` are scratch space.
// `ignored[i]` will be set to true if `cuts[i]` will not be executed in this
// round of recursion.
// `join_targets` will store sequence elements that need to be joined to each other.
// `edge_elements[i]` will store a pointer to the sequence element corresponding to
// edge `cuts[i]`.
void EulerTourTree::BatchCutRecurse(const parlay::sequence<std::pair<int, int>>& cuts,
                                    parlay::sequence<bool>& ignored,
                                    parlay::sequence<Element*>& join_targets,
                                    parlay::sequence<Element*>& edge_elements) {
  const size_t len = cuts.size();
  if (len <= 75) {
    BatchCutSequential(this, cuts);
    return;
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
      Element* uv{edges_.Find(u, v)};
      edge_elements[i] = uv;
      Element* vu{uv->twin_};
      uv->split_mark_ = vu->split_mark_ = true;
    }
  });

  randomness_ = randomness_.next();

  parlay::parallel_for(0, len, [&](size_t i) {
    if (ignored[i]) {
      join_targets[4 * i] = join_targets[4 * i + 2] = nullptr;
      return;
    }
    Element* uv{edge_elements[i]};
    Element* vu{uv->twin_};

    Element* left_target{uv->GetPreviousElement()};
    if (left_target->split_mark_) {
      join_targets[4 * i] = nullptr;
    } else {
      Element* right_target{vu->GetNextElement()};
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
      Element* right_target{uv->GetNextElement()};
      while (right_target->split_mark_) {
        right_target = right_target->twin_->GetNextElement();
      }
      join_targets[4 * i + 2] = left_target;
      join_targets[4 * i + 3] = right_target;
    }
  });

  parlay::parallel_for(0, len, [&](size_t i) {
    if (!ignored[i]) {
      Element* uv{edge_elements[i]};
      Element* vu{uv->twin_};
      uv->SplitWithoutUpdate();
      vu->SplitWithoutUpdate();
      Element* predecessor{uv->GetPreviousElement()};
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
      Element* uv{edge_elements[i]};
      Element* vu{uv->twin_};

      uv->~Element();
      ElementAllocator::free(uv);
      vu->~Element();
      ElementAllocator::free(vu);

      int u, v;
      std::tie(u, v) = cuts[i];
      edges_.Delete(u, v);

      if (join_targets[4 * i] != nullptr) {
        Element::JoinWithoutUpdate(join_targets[4 * i], join_targets[4 * i + 1]);
      }
      if (join_targets[4 * i + 2] != nullptr) {
        Element::JoinWithoutUpdate(join_targets[4 * i + 2], join_targets[4 * i + 3]);
      }
    }
  });

  // We only need to update the elements that participated in joins since the elements that
  // participated in splits were either deleted or participated in joins.
  const auto update_targets{parlay::delayed_seq<Element*>(2 * len, [&join_targets](const size_t i) {
    return join_targets[2 * i];
  })};
  Element::BatchUpdate(update_targets);

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

void EulerTourTree::BatchCut(parlay::sequence<std::pair<int, int>>& cuts) {
  const size_t len = cuts.size();
  if (len <= 75) {
    BatchCutSequential(this, cuts);
    return;
  }

  auto ignored = parlay::sequence<bool>::uninitialized(len);
  auto join_targets = parlay::sequence<Element*>::uninitialized(4 * len);
  auto edge_elements = parlay::sequence<Element*>::uninitialized(len);
  BatchCutRecurse(cuts, ignored, join_targets, edge_elements);
}

}  // namespace parallel_euler_tour_tree

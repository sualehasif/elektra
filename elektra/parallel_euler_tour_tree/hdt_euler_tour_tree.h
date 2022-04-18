#pragma once

#include <parlay/internal/group_by.h>
#include <parlay/primitives.h>

#include "euler_tour_tree.h"
#include "hdt_element.h"
#include "utilities.h"

namespace parallel_euler_tour_tree {

class HdtEulerTourTree : public EulerTourTreeBase<HdtElement> {
 public:
  HdtEulerTourTree(v_int num_vertices);

  // Adds an edge to the forest.
  //
  // `is_level_i_edge` refers to whether this element represents a level-i edge,
  // where "level-i" refers to the level (in the HDT algorithm) of this ETT.
  void Link(v_int u, v_int v, bool is_level_i_edge);
  // Adds edges to the forest, where `is_level_i_edge[i]` corresponds to whether
  // links[i] is a level-i edge.
  template <typename BoolSeq>
  void BatchLink(const parlay::sequence<std::pair<v_int, v_int>>& links, const BoolSeq& is_level_i_edge);
  // Adds edges to the forest, where `is_level_i_edge` corresponds to whether
  // links[i] is a level-i edge.
  void BatchLink(const parlay::sequence<std::pair<v_int, v_int>>& links, bool is_level_i_edge);

  // Get the number of vertices in v's connected component.
  size_t ComponentSize(v_int v) const;
  // Returns all vertices in v's connected component.
  parlay::sequence<v_int> ComponentVertices(v_int v) const;
  // Returns all edges in `v`'s connected component..
  parlay::sequence<std::pair<v_int, v_int>> ComponentEdges(v_int v) const;

  // In v's connected component, return all level-i tree edges and mark them as
  // no longer being level-i tree edges (because they're going to be pushed down
  // to the next level).
  parlay::sequence<std::pair<v_int, v_int>> GetAndClearLevelIEdges(v_int v);

  // Increments vertices[i]'s count of "number of level-i non-tree edges incident to this
  // vertex" to be new_values[i].
  template <typename ElemSeq, typename IntSeq>
  void IncrementNontreeEdgeCounts(const ElemSeq& vertices, const IntSeq& increments);

  // Creates a NontreeEdgeFinder for v's connected component. The
  // NontreeEdgeFinder can be used to find level-i non-tree edges incident to
  // the component.
  NontreeEdgeFinder CreateNontreeEdgeFinder(v_int v) const;

 private:
  using Base = EulerTourTreeBase<HdtElement>;
  using Elem = HdtElement;
};

// A NontreeEdgeFinder finds level-i non-tree edges incident to a particular
// connected component.
//
// The NontreeEdgeFinder is no longer valid if the component is modified by a
// link or a cut, though IncrementNontreeEdgeCounts() modifications are OK.
//
// In the comments in this class, when talking about the number of level-i
// non-tree edges incident to the component, an edge is counted twice if both of
// its endpoints are in the component.
class NontreeEdgeFinder {
 public:
  // The returns the number of level-i non-tree edges incident to the connected
  // component.
  uint64_t NumEdges() const;

  // Among the level-i non-tree edges incident on this component, this function
  // finds the `edges_begin`-th (inclusive) to `edges_end`-th (exclusive) such
  // edges. (It will find fewer edges if there are fewer than `ends_end` such
  // edges.) The function then applies `f` to the vertices in the component that
  // these edges are incident on to allow the caller to iterate over the edges.
  //
  // f should be a function capable of running in concurrently with other copies
  // of itself and should have the interface `void f(v_int vertex_id, v_int
  // begin, v_int end)` where the arguments signify that f should operate upon
  // the `begin`-th (inclusive) to `end`-th (exclusive) level-i non-tree edges
  // of vertex `vertex_id`.
  template <typename Func>
  void ForEachIncidentVertex(uint64_t edges_begin, uint64_t edges_end, Func f) const;

 private:
  friend HdtEulerTourTree;

  NontreeEdgeFinder(const HdtElement* top_element);

  template <typename F>
  void ForEachIncidentVertexImpl(
      F f,
      uint64_t num_desired_edges,
      uint64_t search_offset,
      int level,
      const HdtElement* element) const;

  // A fixed representative element of the list representing the connected component.
  const HdtElement* top_element_;
};

// Given a batch of non-tree edges, update the ETT level-i non-tree edge counts (by
// incrementing the non-tree edge counts for each edge for insertions
// (`is_insertion == true`) or decrementing the counts for each edge for
// deletions (`is_insertion == false`)).
template <typename ElemSeq>
void UpdateNontreeEdges(
    HdtEulerTourTree* ett,
    bool is_insertion,
    const ElemSeq& edges
) {
  auto vertices = parlay::delayed_seq<v_int>(
      2 * edges.size(),
      [&](const size_t i) {
        return i % 2 ? edges[i / 2].first : edges[i / 2].second;
      }
  );
  const auto vertex_counts = parlay::histogram_by_key(vertices);
  const auto vertices_to_increment = parlay::delayed_seq<v_int>(
      vertex_counts.size(),
      [&](const size_t i) {
        return vertex_counts[i].first;
      }
  );
  const int multiplier = is_insertion ? 1 : -1;
  const auto increments = parlay::delayed_seq<int32_t>(
      vertex_counts.size(),
      [&](const size_t i) {
        return multiplier * vertex_counts[i].second;
      }
  );
  ett->IncrementNontreeEdgeCounts(vertices_to_increment, increments);
}

///////////////////////////////////////////////////////////////////////////////
//                  HdtEulerTourTree implementation below.                   //
///////////////////////////////////////////////////////////////////////////////

HdtEulerTourTree::HdtEulerTourTree(v_int num_vertices) : Base{num_vertices} {}

void HdtEulerTourTree::Link(v_int u, v_int v, bool is_level_i_edge) {
  auto [uv, vu] = _internal::AllocEdges<Elem>(u, v);
  new (uv) Elem{randomness_.ith_rand(0), make_pair(u, v), is_level_i_edge};
  new (vu) Elem{randomness_.ith_rand(1), make_pair(v, u), is_level_i_edge};
  randomness_ = randomness_.next();
  Base::Link(u, v, uv, vu);
}

template <typename BoolSeq>
void HdtEulerTourTree::BatchLink(
    const parlay::sequence<std::pair<v_int, v_int>>& links,
    const BoolSeq& is_level_i_edge) {
  const auto construct{[&](const parlay::random& randomness, size_t i, Elem* uv, Elem* vu) {
    const v_int u{links[i].first};
    const v_int v{links[i].second};
    const bool is_level_i{is_level_i_edge[i]};
    new (uv) Elem{randomness.ith_rand(2 * i), make_pair(u, v), is_level_i};
    new (vu) Elem{randomness.ith_rand(2 * i + 1), make_pair(v, u), is_level_i};
  }};
  Base::BatchLink(links, construct);
}

void HdtEulerTourTree::BatchLink(
    const parlay::sequence<std::pair<v_int, v_int>>& links,
    bool is_level_i_edge) {
  BatchLink(links, parlay::delayed_seq<bool>(links.size(), [&](size_t) { return is_level_i_edge; }));
}

size_t HdtEulerTourTree::ComponentSize(v_int v) const {
  return vertices_[v].GetComponentSize();
}

parlay::sequence<v_int> HdtEulerTourTree::ComponentVertices(v_int v) const {
  // TODO(tomtseng) implement this more efficiently if we actually need this
  // function
  parlay::sequence<v_int> vertices;
  vertices.reserve(ComponentSize(v));
  const Elem* const start_element{&vertices_[v]};
  const Elem* curr{start_element};
  do {
    const auto [u, v] = curr->id_;
    if (u == v) {
      vertices.emplace_back(u);
    }
    curr = curr->GetNextElement();
  } while (curr != start_element);
  return vertices;
}

// Returns all edges in `v`'s connected component..
parlay::sequence<std::pair<v_int, v_int>> HdtEulerTourTree::ComponentEdges(v_int v) const {
  // TODO(tomtseng) implement this more efficiently if we actually need this
  // function
  parlay::sequence<std::pair<v_int, v_int>> edges;
  edges.reserve(ComponentSize(v) - 1);
  const Elem* const start_element{&vertices_[v]};
  const Elem* curr{start_element};
  do {
    const auto [u, v] = curr->id_;
    if (u < v) {
      edges.emplace_back(curr->id_);
    }
    curr = curr->GetNextElement();
  } while (curr != start_element);
  return edges;
}

parlay::sequence<std::pair<v_int, v_int>> HdtEulerTourTree::GetAndClearLevelIEdges(v_int v) {
  return vertices_[v].GetAndClearLevelIEdges();
}

template <typename ElemSeq, typename IntSeq>
void HdtEulerTourTree::IncrementNontreeEdgeCounts(
    const ElemSeq& vertices,
    const IntSeq& increments) {
  const auto elements{parlay::delayed_seq<Elem*>(
      vertices.size(),
      [&](size_t i) { return &vertices_[vertices[i]]; })};
  Elem::IncrementNontreeEdgeCounts(elements, increments);
}

NontreeEdgeFinder HdtEulerTourTree::CreateNontreeEdgeFinder(v_int v) const {
  return NontreeEdgeFinder{vertices_[v].FindRepresentative()};
}

///////////////////////////////////////////////////////////////////////////////
//                  NontreeEdgeFinder implementation below.                  //
///////////////////////////////////////////////////////////////////////////////

// Implementation notes:
// - We implement this as separate class from HdtElement just to avoid recalculating
//   top_element_ every time we call ForEachIncidentVertex. This optimization is
//   probably completely insignificant.

NontreeEdgeFinder::NontreeEdgeFinder(const HdtElement* top_element)
  : top_element_{top_element} {
}

uint64_t NontreeEdgeFinder::NumEdges() const {
  uint64_t num_edges{0};
  const HdtElement* curr{top_element_};
  const int level{top_element_->height_ - 1};
  do {
    num_edges += std::get<2>(curr->values_[level]);
    curr = curr->neighbors_[level].next;
  } while (curr != top_element_);
  return num_edges;
}

// Helper function for ForEachIncidentVertex. Starting at `element`, the
// function traverses along `level` until reaching its right parent or returning
// back to `top_element_`. While traversing, the function searches for the first
// `num_desired_edges` level-i non-tree edges after skipping the first
// `search_offset` edges.
template <typename F>
void NontreeEdgeFinder::ForEachIncidentVertexImpl(
    F f,
    uint64_t num_desired_edges,
    uint64_t search_offset,
    int level,
    const HdtElement* element) const {
  if (level <= 6) {
    // run sequentially if we're near the bottom of the list and not doing as
    // much work per thread
    do {
      const uint64_t num_edges{std::get<2>(element->values_[level])};
      if (search_offset >= num_edges) {
        search_offset -= num_edges;
      } else {
        if (level == 0) {
          f(element->id_.first, search_offset, std::min<uint64_t>(search_offset + num_desired_edges, num_edges));
        } else {
          ForEachIncidentVertexImpl(f, num_desired_edges, search_offset, level - 1, element);
        }
        num_desired_edges -= std::min<uint64_t>(num_edges - search_offset, num_desired_edges);
        search_offset = 0;
      }
      element = element->neighbors_[level].next;
    } while (element->height_ <= level + 1 && num_desired_edges > 0 && element != top_element_);
    return;
  }

  // run in parallel
  struct LoopState {
    const HdtElement* element;
    uint64_t num_desired_edges;
    uint64_t search_offset;
  } loop_state = { element, num_desired_edges, search_offset };
  const auto loop_condition{[&](const LoopState& state) {
    return state.element->height_ <= level + 1 &&
      num_desired_edges > 0 &&
      element != top_element_;
  }};
  const auto loop_action{[&](const LoopState& state) {
    const uint64_t num_edges{std::get<2>(state.element->values_[level])};
    if (state.search_offset >= num_edges) {
      return;
    }
    ForEachIncidentVertexImpl(f, state.num_desired_edges, state.search_offset, level - 1, state.element);
  }};
  const auto loop_update{[&](LoopState state) -> LoopState {
    uint64_t num_edges{std::get<2>(state.element->values_[level])};
    if (state.search_offset >= num_edges) {
      state.search_offset -= num_edges;
    } else {
      state.num_desired_edges -= std::min<uint64_t>(num_edges - state.search_offset, state.num_desired_edges);
      state.search_offset = 0;
    }
    state.element = state.element->neighbors_[level].next;
    return state;
  }};
  constexpr bool is_do_while{true};
  elektra::ParallelWhile(loop_condition, loop_action, loop_update, std::move(loop_state), is_do_while);
}

template <typename F>
void NontreeEdgeFinder::ForEachIncidentVertex(uint64_t l, uint64_t r, F f) const {
  if (l >= r) {
    return;
  }
  const int level{top_element_->height_ - 1};
  ForEachIncidentVertexImpl(f, r - l, l, level, top_element_);
}

}  // namespace parallel_euler_tour_tree

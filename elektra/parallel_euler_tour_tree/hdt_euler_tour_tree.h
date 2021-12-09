#pragma once

#include "euler_tour_tree.h"
#include "hdt_element.h"

namespace parallel_euler_tour_tree {

class HdtEulerTourTree : public EulerTourTreeBase<HdtElement> {
 public:
  HdtEulerTourTree(int num_vertices);

  // Adds an edge to the forest.
  //
  // `is_level_i_edge` refers to whether this element represents a level-i edge,
  // where "level-i" refers to the level (in the HDT algorithm) of this ETT.
  void Link(int u, int v, bool is_level_i_edge);
  // Adds edges to the forest, where `is_level_i_edge[i]` corresponds to whether
  // links[i] is a level-i edge.
  template <typename BoolSeq>
  void BatchLink(
      const parlay::sequence<std::pair<int, int>>& links,
      const BoolSeq& is_level_i_edge);

  // Get the number of vertices in v's connected component.
  size_t ComponentSize(int v) const;

  // In v's connected component, return all level-i tree edges and mark them as
  // no longer being level-i tree edges (because they're going to be pushed down
  // to the next level).
  parlay::sequence<std::pair<int, int>> ClearLevelIEdges(int v);

  // Updates vertices[i]'s count of "number of level-i non-tree edges incident to this
  // vertex" to be new_values[i].
  template <typename IntSeq>
  void UpdateNontreeEdgeCounts(
      const parlay::sequence<int>& vertices,
      IntSeq&& new_values);

  // Creates a NontreeEdgeFinder for v's connected component. The
  // NontreeEdgeFinder can be used to find level-i non-tree edges incident to
  // the component
  NontreeEdgeFinder CreateNontreeEdgeFinder(int v) const;

 private:
  using Base = EulerTourTreeBase<HdtElement>;
  using Elem = HdtElement;
};

// A NontreeEdgeFinder finds level-i non-tree edges incident to a particular
// connected component.
//
// The NontreeEdgeFinder is no longer valid if the component is modified by a
// link, a cut, or a UpdateNontreeEdgeCounts().
class NontreeEdgeFinder {
 public:
  // The returns the number of level-i non-tree edges incident to the connected
  // component.
  uint64_t NumEdges() const;

  // This function finds vertices in the component that have incident level-i
  // non-tree edges and applies `f` to the vertices. The edges found are the
  // first `num_edges` edges after skipping the first `search_offset` edges.
  //
  // Params:
  // - f: A parallel function to apply to vertices that have incident edges. The
  //   interface of f should be `void f(uint32_t vertex_id, uint32_t begin,
  //   uint32_t end)` where the arguments signify that f should operate upon the
  //   `begin`-th (inclusive) to `end`-th (exclusive) level-i non-tree edges of
  //   vertex `vertex_id`.
  // - num_edges: Specifies how many level-i non-tree edges to search for.
  // - search_offset: Specifies prefix of non-tree edges to skip in the search.
  //   For instance, suppose you've already found 50 edges using
  //   FindNontreeIncidentVertices(f, 50, 0). Then you can find the next 20 edges
  //   using FindNontreeIncidentVertices(f, 20, 50).
  template <typename Func>
  void ForEachIncidentVertex(Func f, uint64_t num_desired_edges, uint64_t search_offset) const;

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
  // The number of level-i non-tree edges incident to the component.
  uint64_t num_incident_edges_{0};
};

///////////////////////////////////////////////////////////////////////////////
//                  HdtEulerTourTree Implementation below.                   //
///////////////////////////////////////////////////////////////////////////////

HdtEulerTourTree::HdtEulerTourTree(int num_vertices) : Base{num_vertices} {}

void HdtEulerTourTree::Link(int u, int v, bool is_level_i_edge) {
  Elem* uv = ElementAllocator::alloc();
  Elem* vu = ElementAllocator::alloc();
  new (uv) Elem{randomness_.ith_rand(0), make_pair(u, v), is_level_i_edge};
  new (vu) Elem{randomness_.ith_rand(1), make_pair(v, u), is_level_i_edge};
  randomness_ = randomness_.next();
  Base::Link(u, v, uv, vu);
}

template <typename BoolSeq>
void HdtEulerTourTree::BatchLink(
    const parlay::sequence<std::pair<int, int>>& links,
    const BoolSeq& is_level_i_edge) {
  const auto construct{[&](parlay::random* randomness, size_t i, Elem* uv, Elem* vu) {
    const int u{links[i].first};
    const int v{links[i].second};
    const bool is_level_i{is_level_i_edge[i]};
    new (uv) Elem{randomness->ith_rand(2 * i), make_pair(u, v), is_level_i};
    new (vu) Elem{randomness->ith_rand(2 * i + 1), make_pair(v, u), is_level_i};
  }};
  Base::BatchLink(links, construct);
}

size_t HdtEulerTourTree::ComponentSize(int v) const {
  return vertices_[v].GetComponentSize();
}

parlay::sequence<std::pair<int, int>> HdtEulerTourTree::ClearLevelIEdges(int v) {
  return vertices_[v].ClearLevelIEdges();
}

template <typename IntSeq>
void HdtEulerTourTree::UpdateNontreeEdgeCounts(
    const parlay::sequence<int>& vertices,
    IntSeq&& new_values) {
  const auto elements{parlay::delayed_seq<Elem*>(
      vertices.size(),
      [&](size_t i) { return &vertices_[vertices[i]]; })};
  Elem::UpdateNontreeEdgeCounts(elements, std::move(new_values));
}

NontreeEdgeFinder HdtEulerTourTree::CreateNontreeEdgeFinder(int v) const {
  return NontreeEdgeFinder{vertices_[v].FindRepresentative()};
}

///////////////////////////////////////////////////////////////////////////////
//                  NontreeEdgeFinder Implementation below.                  //
///////////////////////////////////////////////////////////////////////////////

// Implementation notes:
// - We implement this as separate class from HdtElement just to avoid recalculating
//   top_element_ every time we call ForEachIncidentVertex. This optimization is
//   probably completely insignificant.

NontreeEdgeFinder::NontreeEdgeFinder(const HdtElement* top_element)
  : top_element_{top_element} {
  const HdtElement* curr{top_element};
  const int level{top_element_->height_ - 1};
  do {
    num_incident_edges_ += std::get<2>(curr->values_[level]);
    curr = curr->neighbors_[level].next;
  } while (curr != top_element);
}

uint64_t NontreeEdgeFinder::NumEdges() const {
  return num_incident_edges_;
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
void NontreeEdgeFinder::ForEachIncidentVertex(
    F f,
    uint64_t num_desired_edges,
    uint64_t search_offset) const {
  const int level{top_element_->height_ - 1};
  ForEachIncidentVertexImpl(f, num_desired_edges, search_offset, level, top_element_);
}

}  // namespace parallel_euler_tour_tree
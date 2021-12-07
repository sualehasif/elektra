#pragma once

#include <parlay/delayed_sequence.h>
#include <parlay/primitives.h>
#include <parlay/sequence.h>

#include "element.h"
#include "utilities/while_loop.h"

namespace parallel_euler_tour_tree {

namespace _internal {

struct HdtAugmentation {
  // 0. Augment each vertex element to get "number of vertices in this list".
  // 1. Augment each edge element to get "number of level-i edges in this list".
  // 2. Augment each vertex element to get "number of level-i non-tree edges
  //    incident to vertices in this list".
  using T = std::tuple<uint32_t, uint32_t, uint64_t>;
  static T f(const T& a, const T& b) {
    return {
      std::get<0>(a) + std::get<0>(b),
      std::get<1>(a) + std::get<1>(b),
      std::get<2>(a) + std::get<2>(b),
    };
  }
};

template <typename Tuple, int Index>
struct IntTupleGetter {
  using T = std::tuple_element_t<Index, Tuple>;
  static inline T& Get(Tuple& elem) {
    return std::get<Index>(elem);
  }
  static T f(T a, T b) {
    return a + b;
  }
};

// Structs for use in `HdtElement::BatchUpdate()` that apply `HdtAugmentation`
// at a single index of the tuple `HdtAugmentation::T`.
using IsLevelIEdgeGetter = IntTupleGetter<HdtAugmentation::T, 1>;
using NontreeEdgeCountGetter = IntTupleGetter<HdtAugmentation::T, 2>;

}  // namespace _internal

class NontreeEdgeFinder;

// Skip list element augmented to be useful for HDT-like dynamic connectivity
// algorithm. (HDT = Holm + De Lichtenberg + Thorup, referring to authors of a
// dynamic connectivity algorithm.)
//
// Most of the functions in this class assume that the list is circular.
class HdtElement : public ElementBase<HdtElement, _internal::HdtAugmentation> {
 public:
  using Base = ElementBase<HdtElement, _internal::HdtAugmentation>;

  // `is_level_i_edge` refers to whether this element represents a level-i tree
  // edge, where "level-i" refers to be the level (in the HDT algorithm) of
  // whatever ETT this skip list lives in.
  explicit HdtElement(size_t random_int, std::pair<int, int>&& id, bool is_level_i_edge);

  // Get the number of vertices in the sequence that contains this element.
  size_t GetComponentSize() const;

  // TODO(tomtseng) comment
  NontreeEdgeFinder CreateNontreeEdgeFinder() const;

  // In the list that contains this element, return all level-i tree edges
  // and mark them as no longer being level-i tree edges (because they're going
  // to be pushed down to the next level).
  parlay::sequence<std::pair<int, int>> ClearLevelIEdges();

  // Updates the counts of "number of level-i non-tree edges incident to this
  // vertex" for each element (which should represent a vertex) in `elements`.
  //
  // `IntSeq` should behave like `parlay::sequence<int>`.
  template <typename IntSeq>
  static void UpdateNontreeEdgeCounts(
      const parlay::sequence<HdtElement*>& elements,
      IntSeq&& new_values);

 private:
  friend Base::Base::Base;  // make parallel_skip_list::ElementBase a friend
  friend NontreeEdgeFinder;

  // Helper function for ClearLevelIEdges. Gets edges held in descendants of
  // this element at writes them into the sequence starting at the offset.
  void GetLevelIEdgesBelow(parlay::sequence<HdtElement*>* s, int level, uint32_t offset) const;
};

// TODO(tomtseng): comment and use this
//
// no longer valid if u join/link//update the underlying component
//
// impl note: we just design this this way so we don't fetch the top_element
// repeatedly lol
class NontreeEdgeFinder {
 public:
  // TODO(tomtseng): comment
  //
  //  f should take a vertex_id, start_index, end_index, and should be able to
  //  run in parallel
  //
  //  should return how many edges remaining
  //
  // In the skip list that contains this element, return a list of vertices adjacent
  // to `num_edges` level-i non-tree edges.
  //
  // Params:
  // - num_edges: Specifies how many level-i non-tree edges to search for
  // - search_offset: Specifies prefix of non-tree edges to skip in the search.
  //     For instance, suppose you've already found 50 edges using
  //     FindNontreeIncidentVertices(50, 0). Then you can find the next 100
  //     edges using FindNontreeIncidentVertices(100, 50).
  //
  template <typename Func>
  void ForEachIncidentVertex(Func f, size_t num_non_tree_edges, size_t search_offset) const;

  // TODO(tomtseng): comment
  uint64_t NumEdges() const;

 private:
  friend HdtElement;

  NontreeEdgeFinder(const HdtElement* top_element);

  const HdtElement* top_element_;
  // TODO(tomtseng): populate and comment
  uint64_t num_incident_edges_{0};
};

///////////////////////////////////////////////////////////////////////////////
//                           Implementation below.                           //
///////////////////////////////////////////////////////////////////////////////

HdtElement::HdtElement(size_t random_int, std::pair<int, int>&& id, bool is_level_i_edge)
  // Using std::move(id) and also accessing id's fields is sketchy but should
  // be OK since std::move's immediate effect is really just a cast.
  : Base{random_int, std::move(id), {id.first == id.second, id.first < id.second && is_level_i_edge, 0}} {}

size_t HdtElement::GetComponentSize() const {
  HdtElement* const top_element{FindRepresentative()};
  const int level{top_element->height_ - 1};

  size_t num_vertices{0};
  HdtElement* curr{top_element};
  do {
    num_vertices += std::get<0>(curr->values_[level]);
    curr = curr->neighbors_[level].next;
  } while (curr != top_element);
  return num_vertices;
}

NontreeEdgeFinder HdtElement::CreateNontreeEdgeFinder() const {
  return NontreeEdgeFinder{FindRepresentative()};
}

void HdtElement::GetLevelIEdgesBelow(parlay::sequence<HdtElement*>* s, int level, uint32_t offset) const {
  if (level == 0) {
    if (std::get<1>(values_[0])) {
      (*s)[offset] = const_cast<HdtElement*>(this);
    }
    return;
  }
  if (std::get<1>(values_[level]) == 0) {
    return;
  }

  const HdtElement* curr{this};
  if (level <= 6) {
    // run sequentially once we're near the bottom of the list and not doing as
    // much work per thread
    do {
      curr->GetLevelIEdgesBelow(s, level - 1, offset);
      offset += std::get<1>(curr->values_[level - 1]);
      curr = curr->neighbors_[level - 1].next;
    } while (curr->height_ < level + 1);
  } else {  // run in parallel
    struct LoopState {
      const HdtElement* curr;
      uint32_t offset;
    } loop_state = { this, offset };
    const auto loop_condition{[&](const LoopState& state) { return state.curr->height_ < level + 1; }};
    const auto loop_action{[&](const LoopState& state) {
      state.curr->GetLevelIEdgesBelow(s, level - 1, state.offset);
    }};
    const auto loop_update{[&](const LoopState& state) -> LoopState {
      return { state.curr->neighbors_[level - 1].next, state.offset + std::get<1>(state.curr->values_[level - 1]) };
    }};
    constexpr bool is_do_while{true};
    elektra::ParallelWhile(loop_condition, loop_action, loop_update, std::move(loop_state), is_do_while);
  }
}

parlay::sequence<std::pair<int, int>> HdtElement::ClearLevelIEdges() {
  HdtElement* const top_element{FindRepresentative()};
  const int level{top_element->height_ - 1};

  size_t num_edges{0};
  {
    HdtElement* curr = top_element;
    do {
      num_edges += std::get<1>(curr->values_[level]);
      curr = curr->neighbors_[level].next;
    } while (curr != top_element);
  }

  parlay::sequence<HdtElement*> edge_elements(num_edges);
  struct LoopState {
    const HdtElement* curr;
    uint32_t offset;
  } loop_state = { top_element, 0 };
  const auto loop_condition{[&](const LoopState& state) { return state.curr != top_element; }};
  const auto loop_action{[&](const LoopState& state) {
    state.curr->GetLevelIEdgesBelow(&edge_elements, level, state.offset);
  }};
  const auto loop_update{[&](const LoopState& state) -> LoopState {
    return { state.curr->neighbors_[level].next, state.offset + std::get<1>(state.curr->values_[level]) };
  }};
  constexpr bool is_do_while{true};
  elektra::ParallelWhile(loop_condition, loop_action, loop_update, std::move(loop_state), is_do_while);

  BatchUpdate<_internal::IsLevelIEdgeGetter>(
      edge_elements,
      parlay::delayed_seq<int>(num_edges, [](size_t) { return 0; }));
  parlay::sequence<std::pair<int, int>> edges{
    parlay::map(edge_elements, [](const HdtElement* elem) { return elem->id_; })
  };
  return edges;
}

template <typename IntSeq>
void HdtElement::UpdateNontreeEdgeCounts(
    const parlay::sequence<HdtElement*>& elements,
    IntSeq&& new_values) {
  BatchUpdate<_internal::NontreeEdgeCountGetter>(elements, std::move(new_values));
}

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

template <typename F>
void NontreeEdgeFinder::ForEachIncidentVertex(
    F f,
    size_t num_non_tree_edges,
    size_t search_offset) const {
  const int level{top_element_->height_ - 1};
  exit(-1);
  // TODO(tomtseng) impl
}

}  // namespace parallel_euler_tour_tree

#pragma once

#include <parlay/delayed_sequence.h>
#include <parlay/primitives.h>
#include <parlay/sequence.h>

#include "element.h"
#include "utilities/while_loop.h"
#include "utilities.h"

namespace parallel_euler_tour_tree {

namespace _internal {

struct HdtAugmentation {
  // 0. Augment each vertex element to get "number of vertices in this list".
  // 1. Augment each edge element to get "number of level-i edges in this list".
  // 2. Augment each vertex element to get "number of level-i non-tree edges
  //    incident to vertices in this list".
  using T = std::tuple<v_int, v_int, uint64_t>;
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
  // edge, where "level-i" refers to the level (in the HDT algorithm) of
  // whatever ETT this skip list lives in.
  explicit HdtElement(size_t random_int, std::pair<v_int, v_int>&& id, bool is_level_i_edge = false);

  // Gets the number of vertices in the list that contains this element.
  size_t GetComponentSize() const;

  // Creates a NontreeEdgeFinder based on the list that contains this element.
  // The NontreeEdgeFinder can be used to find level-i non-tree edges incident
  // to vertices in this list.
  NontreeEdgeFinder CreateNontreeEdgeFinder() const;

  // In the list that contains this element, returns all level-i tree edges
  // and mark them as no longer being level-i tree edges (because they're going
  // to be pushed down to the next level).
  parlay::sequence<std::pair<v_int, v_int>> GetAndClearLevelIEdges();

  // For each HdtElement elements[i] (which should represent a vertex), increments
  // its count of "number of level-i non-tree edges incident to this element" to
  // be new_values[i].
  template <typename ElemSeq, typename IntSeq>
  static void IncrementNontreeEdgeCounts(const ElemSeq& elements, const IntSeq& increment);

 private:
  friend Base::Base::Base;  // make parallel_skip_list::ElementBase a friend
  friend NontreeEdgeFinder;

  void GetLevelIEdgesBelow(parlay::sequence<HdtElement*>* s, int level, v_int offset) const;
};

///////////////////////////////////////////////////////////////////////////////
//                           Implementation below.                           //
///////////////////////////////////////////////////////////////////////////////

HdtElement::HdtElement(size_t random_int, std::pair<v_int, v_int>&& id, bool is_level_i_edge)
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

// Helper function for GetAndClearLevelIEdges. Gets edges held in descendants of
// this element and writes them into the sequence starting at the offset.
void HdtElement::GetLevelIEdgesBelow(parlay::sequence<HdtElement*>* s, int level, v_int offset) const {
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
      v_int offset;
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

parlay::sequence<std::pair<v_int, v_int>> HdtElement::GetAndClearLevelIEdges() {
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
    v_int offset;
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
      parlay::delayed_seq<v_int>(num_edges, [](size_t) { return 0; }));

  // edge_elements = parlay::filter(
  //     edge_elements, [](const HdtElement* e) { return e != nullptr; });

  parlay::sequence<std::pair<v_int, v_int>> edges{
    parlay::map(edge_elements, [](const HdtElement* elem) { return elem->id_; })
  };
  return edges;
}

template <typename ElemSeq, typename IntSeq>
void HdtElement::IncrementNontreeEdgeCounts(const ElemSeq& elements, const IntSeq& increments) {
  BatchIncrement<_internal::NontreeEdgeCountGetter>(elements, increments);
}

}  // namespace parallel_euler_tour_tree

#pragma once

#include <parlay/monoid.h>
#include <parlay/sequence.h>

#include <limits>
#include <utility>

#include "parallel_skip_list/augmented_skip_list.h"
#include "utilities/while_loop.h"
#include "utilities.h"

namespace parallel_euler_tour_tree {

// Skip list element for ETT with augmentation function specified by Func.
template <typename Derived, typename Func>
class ElementBase : public parallel_skip_list::AugmentedElementBase<Derived, Func> {
 public:
  using Base = parallel_skip_list::AugmentedElementBase<Derived, Func>;
  using Value = typename Func::T;

  ElementBase(size_t random_int, std::pair<v_int, v_int>&& id, const Value& value);

  // Returns a representative vertex from the sequence the element lives in.
  // Whereas `FindRepresentative()` returns a representative element that might
  // be an edge, this function returns a vertex ID. For implementation
  // simplicity, this function also assumes that the sequence is circular (as is
  // the case for a sequence representing an ETT component that is not currently
  // performing a join or split).
  v_int FindRepresentativeVertex() const;

  // Returns an estimate of the amount of memory in bytes this element occupies
  // beyond what's captured by sizeof(Element).
  size_t AllocatedMemorySize() const;

  // If this element represents a vertex v, then id == (v, v). Otherwise if
  // this element represents a directed edge (u, v), then id == (u,v).
  std::pair<v_int, v_int> id_;
  // When batch splitting, we mark this as `true` for an edge that we will
  // splice out in the current round of recursion.
  bool split_mark_{false};

 protected:
  friend Base;
  using Base::height_;
  using Base::neighbors_;
  using Base::values_;
};

// Example of an ETT skip list element augmented with the ability to
// - fetch component sizes
// - fetch all elements in a list representing (u, v) with u < v.
class Element : public ElementBase<Element, parlay::addm<v_int>> {
  using Base = ElementBase<Element, parlay::addm<v_int>>;

 public:
  Element(size_t random_int, std::pair<v_int, v_int>&& id)
    // Augment with the function "count the number of edges (u, v) such that u <
    // v in this list".
    //
    // Using std::move(id) and also accessing id's fields is sketchy but should
    // be OK since std::move's immediate effect is really just a cast.
      : Base{random_int, std::move(id), id.first < id.second} {}

  // Get all edges {u, v} in the sequence that contains this element, assuming
  // that the sequence represents an ETT component.
  parlay::sequence<std::pair<v_int, v_int>> GetEdges() const;
  // Get the number of vertices in the sequence that contains this element.
  size_t GetComponentSize() const;

 private:
  friend Base::Base::Base;  // make parallel_skip_list::ElementBase a friend

  // Gets edges held in descendants of this element and writes them into
  // the sequence starting at the offset. `values_` needs to be up to date.
  void GetEdgesBelow(parlay::sequence<std::pair<v_int, v_int>>* s, int level, v_int offset) const;
};

template <typename D, typename F>
ElementBase<D, F>::ElementBase(
    size_t random_int,
    std::pair<v_int, v_int>&& id,
    const Value& value)
    : Base{random_int, value}
    , id_{std::move(id)} {}

template <typename D, typename F>
v_int ElementBase<D, F>::FindRepresentativeVertex() const {
  const D* current_element{static_cast<const D*>(this)};
  const D* seen_element{nullptr};
  int current_level{current_element->height_ - 1};

  // walk up while moving forward to find top level
  while (seen_element != current_element) {
    if (seen_element == nullptr) {
      seen_element = current_element;
    }
    current_element = current_element->neighbors_[current_level].next;
    const int top_level{current_element->height_ - 1};
    if (current_level < top_level) {
      current_level = top_level;
      seen_element = nullptr;
    }
  }

  // look for minimum ID vertex in top level, or try again at lower levels if no
  // vertex is found
  v_int min_vertex{std::numeric_limits<v_int>::max()};
  while (current_level >= 0 && min_vertex == std::numeric_limits<v_int>::max()) {
    do {
      const std::pair<v_int, v_int>& id{current_element->id_};
      if (id.first == id.second && id.first < min_vertex) {
        min_vertex = id.first;
      }
      current_element = current_element->neighbors_[current_level].next;
    } while (current_element != seen_element);
    current_level--;
  }
  return min_vertex;
}

template <typename D, typename F>
size_t ElementBase<D, F>::AllocatedMemorySize() const {
  // Estimate size of neighbors_ and values_, knowing that they're allocated
  // using concurrent_array_allocator.
  return (sizeof(neighbors_[0]) + sizeof(values_[0])) * (1 << parlay::log2_up(height_));
}

void Element::GetEdgesBelow(parlay::sequence<std::pair<v_int, v_int>>* s, int level, v_int offset) const {
  if (level == 0) {
    if (values_[0]) {
      (*s)[offset] = id_;
    }
    return;
  }
  if (values_[level] == 0) {
    return;
  }

  const Element* curr{this};
  if (level <= 6) {
    // run sequentially once we're near the bottom of the list and not doing as
    // much work per thread
    do {
      curr->GetEdgesBelow(s, level - 1, offset);
      offset += curr->values_[level - 1];
      curr = curr->neighbors_[level - 1].next;
    } while (curr->height_ < level + 1);
  } else {  // run in parallel
    struct LoopState {
      const Element* curr;
      v_int offset;
    } loop_state = { this, offset };
    const auto loop_condition{[&](const LoopState& state) { return state.curr->height_ < level + 1; }};
    const auto loop_action{[&](const LoopState& state) {
      state.curr->GetEdgesBelow(s, level - 1, state.offset);
    }};
    const auto loop_update{[&](const LoopState& state) -> LoopState {
      return { state.curr->neighbors_[level - 1].next, state.offset + state.curr->values_[level - 1] };
    }};
    constexpr bool is_do_while{true};
    elektra::ParallelWhile(loop_condition, loop_action, loop_update, std::move(loop_state), is_do_while);
  }
}

parlay::sequence<std::pair<v_int, v_int>> Element::GetEdges() const {
  Element* const top_element{FindRepresentative()};
  const int level = top_element->height_ - 1;

  size_t num_edges{0};
  {
    Element* curr = top_element;
    do {
      num_edges += curr->values_[level];
      curr = curr->neighbors_[level].next;
    } while (curr != top_element);
  }

  parlay::sequence<std::pair<v_int, v_int>> edges(num_edges);

  struct LoopState {
    const Element* curr;
    v_int offset;
  } loop_state = { top_element, 0 };
  const auto loop_condition{[&](const LoopState& state) { return state.curr != top_element; }};
  const auto loop_action{[&](const LoopState& state) {
    state.curr->GetEdgesBelow(&edges, level, state.offset);
  }};
  const auto loop_update{[&](const LoopState& state) -> LoopState {
    return { state.curr->neighbors_[level].next, state.offset + state.curr->values_[level] };
  }};
  constexpr bool is_do_while{true};
  elektra::ParallelWhile(loop_condition, loop_action, loop_update, std::move(loop_state), is_do_while);

  return edges;
}

size_t Element::GetComponentSize() const {
  Element* const top_element{FindRepresentative()};
  const int level = top_element->height_ - 1;

  size_t num_edges{0};
  Element* curr = top_element;
  do {
    num_edges += curr->values_[level];
    curr = curr->neighbors_[level].next;
  } while (curr != top_element);

  const size_t num_vertices = num_edges + 1;
  return num_vertices;
}

}  // namespace parallel_euler_tour_tree

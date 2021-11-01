#pragma once

#include "parallel_skip_list/skip_list.h"

#include <algorithm>

#include <parlay/parallel.h>
#include <parlay/sequence.h>
#include <parlay/utilities.h>

namespace parallel_skip_list {

// Batch-parallel augmented skip list. Currently, the augmentation is
// hardcoded to the sum function with the value 1 assigned to each element. As
// such, `GetSum()` returns the size of the list.
//
// TODO(tomtseng): Allow user to pass in their own arbitrary associative
// augmentation functions. The contract for `GetSum` on a cyclic list should be
// that the function will be applied starting from `this`, because where we
// begin applying the function matters for non-commutative functions.
class AugmentedElement : private ElementBase<AugmentedElement> {
  friend class ElementBase<AugmentedElement>;

 public:
  // See comments on `ElementBase<>`.
  AugmentedElement();
  explicit AugmentedElement(size_t random_int);
  ~AugmentedElement();

  // Can run concurrently with other `JoinWithoutUpdate` calls, but augmented
  // values must be updated separately afterwards.
  static void JoinWithoutUpdate(AugmentedElement* left, AugmentedElement* right);
  // Can run concurrently with other `SplitWithoutUpdate` calls, but augmented
  // values must be updated separately afterwards.
  AugmentedElement* SplitWithoutUpdate();

  // For each `{left, right}` in the `len`-length array `joins`, concatenate the
  // list that `left` lives in to the list that `right` lives in.
  //
  // `left` must be the last node in its list, and `right` must be the first
  // node of in its list. Each `left` must be unique, and each `right` must be
  // unique.
  static void BatchJoin(
      const parlay::sequence<std::pair<AugmentedElement*, AugmentedElement*>>& joins);

  // For each `v` in the `len`-length array `splits`, split `v`'s list right
  // after `v`.
  static void BatchSplit(const parlay::sequence<AugmentedElement*>& splits);

  // For each `i`=0,1,...,`len`-1, assign value `new_values[i]` to element
  // `elements[i]`.
  static void BatchUpdate(
      const parlay::sequence<AugmentedElement*>& elements,
      const parlay::sequence<int>& new_values);

  // Get the result of applying the augmentation function over the subsequence
  // between `left` and `right` inclusive.
  //
  // `left` and `right` must live in the same list, and `left` must precede
  // `right` in the list.
  //
  // This function does not modify the data structure, so it may run
  // concurrently with other `GetSubsequenceSum` calls and const function calls.
  static int GetSubsequenceSum(const AugmentedElement* left, const AugmentedElement* right);

  // Get result of applying the augmentation function over the whole list that
  // the element lives in.
  int GetSum() const;

  using ElementBase<AugmentedElement>::FindRepresentative;
  using ElementBase<AugmentedElement>::GetPreviousElement;
  using ElementBase<AugmentedElement>::GetNextElement;

 private:
  static int* AllocateValues(int height);

  static void UpdateTopDownImpl(int level, AugmentedElement* curr, bool is_loop_start = true);
  // Update aggregate value of node and clear `join_update_level` after joins.
  void UpdateTopDown(int level);
  void UpdateTopDownSequential(int level);

  static concurrent_array_allocator::Allocator<int> values_allocator_;

  int* values_;
  // When updating augmented values, this marks the lowest index at which the
  // `values_` needs to be updated.
  int update_level_;
};

///////////////////////////////////////////////////////////////////////////////
//                           Implementation below.                           //
///////////////////////////////////////////////////////////////////////////////

namespace _internal {

constexpr int NA{-1};

template <typename T>
inline bool write_min(T* a, T b) {
  T c;
  bool r = 0;
  do {
    c = *a;
  } while (b < c && !(r = CAS(a, c, b)));
  return r;
}

}  // namespace _internal

concurrent_array_allocator::Allocator<int> AugmentedElement::values_allocator_{};

int* AugmentedElement::AllocateValues(int height) {
  int* values{values_allocator_.Allocate(height)};
  for (int i = 0; i < height; i++) {
    values[i] = 1;
  }
  return values;
}

AugmentedElement::AugmentedElement() :
  ElementBase<AugmentedElement>{}, update_level_{_internal::NA} {
  values_ = AllocateValues(height_);
}

AugmentedElement::AugmentedElement(size_t random_int) :
  ElementBase<AugmentedElement>{random_int}, update_level_{_internal::NA} {
  values_ = AllocateValues(height_);
}

AugmentedElement::~AugmentedElement() {
  values_allocator_.Free(values_, height_);
}

void AugmentedElement::JoinWithoutUpdate(AugmentedElement* left, AugmentedElement* right) {
  Join(left, right);
}

AugmentedElement* AugmentedElement::SplitWithoutUpdate() {
  return Split();
}

void AugmentedElement::UpdateTopDownSequential(int level) {
  if (level == 0) {
    if (height_ == 1) {
      update_level_ = _internal::NA;
    }
    return;
  }

  if (update_level_ < level) {
    UpdateTopDownSequential(level - 1);
  }
  int sum{values_[level - 1]};
  AugmentedElement* curr{neighbors_[level - 1].next};
  while (curr != nullptr && curr->height_ < level + 1) {
    if (curr->update_level_ != _internal::NA && curr->update_level_ < level) {
      curr->UpdateTopDownSequential(level - 1);
    }
    sum += curr->values_[level - 1];
    curr = curr->neighbors_[level - 1].next;
  }
  values_[level] = sum;

  if (height_ == level + 1) {
    update_level_ = _internal::NA;
  }
}

// Helper function for AugmentedElement::UpdateTopDown that updates the
// augmented values of the descendants of `curr`.
void AugmentedElement::UpdateTopDownImpl(
    int level,
    AugmentedElement* curr,
    bool is_loop_start)
{
  while (is_loop_start || (curr != nullptr && curr->height_ < level + 1)) {
    if (curr->update_level_ != _internal::NA && curr->update_level_ < level) {
      parlay::par_do(
          [&]() { curr->UpdateTopDown(level - 1); },
          [&]() { UpdateTopDownImpl(level, curr->neighbors_[level - 1].next, false); }
      );
      return;
    }
    curr = curr->neighbors_[level - 1].next;
  }
}


// `v.UpdateTopDown(level)` updates the augmented values of descendants of `v`'s
// `level`-th node. `update_level_` is used to determine what nodes need
// updating. `update_level_` is reset to `NA` for all traversed nodes at end of
// this function.
void AugmentedElement::UpdateTopDown(int level) {
  if (level <= 6) {
    UpdateTopDownSequential(level);
    return;
  }

  UpdateTopDownImpl(level, this);

  // Now that children have correct augmented valeus, update self's augmented
  // value.
  int sum{values_[level - 1]};
  AugmentedElement* curr = neighbors_[level - 1].next;
  while (curr != nullptr && curr->height_ < level + 1) {
    sum += curr->values_[level - 1];
    curr = curr->neighbors_[level - 1].next;
  }
  values_[level] = sum;

  if (height_ == level + 1) {
    update_level_ = _internal::NA;
  }
}

// If `new_values` is non-empty, for each `i`=0,1,...,`len`-1, assign value
// `new_vals[i]` to element `elements[i]`.
//
// If `new_values` is empty, update the augmented values of the ancestors of
// `elements`, where the "ancestors" of element `v` refer to `v`,
// `v->FindLeftParent(0)`, `v->FindLeftParent(0)->FindLeftParent(1)`,
// `v->FindLeftParent(0)->FindLeftParent(2)`, and so on. This functionality is
// used privately to keep the augmented values correct when the list has
// structurally changed.
void AugmentedElement::BatchUpdate(
    const parlay::sequence<AugmentedElement*>& elements,
    const parlay::sequence<int>& new_values) {
  const size_t len{elements.size()};
  if (!new_values.empty()) {
    parlay::parallel_for(0, len, [&](size_t i) {
      elements[i]->values_[0] = new_values[i];
    });
  }

  // The nodes whose augmented values need updating are the ancestors of
  // `elements`. Some nodes may share ancestors. `top_nodes` will contain,
  // without duplicates, the set of all ancestors of `elements` with no left
  // parents. From there we can walk down from those ancestors to update all
  // required augmented values.
  auto top_nodes{parlay::sequence<AugmentedElement*>::uninitialized(len)};

  parlay::parallel_for(0, len, [&](const size_t i) {
    int level{0};
    AugmentedElement* curr{elements[i]};
    while (true) {
      int curr_update_level{curr->update_level_};
      if (curr_update_level == _internal::NA && CAS(&curr->update_level_, _internal::NA, level)) {
        level = curr->height_ - 1;
        AugmentedElement* parent{curr->FindLeftParent(level)};
        if (parent == nullptr) {
          top_nodes[i] = curr;
          break;
        } else {
          curr = parent;
          level++;
        }
      } else {
        // Someone other execution is shares this ancestor and has already
        // claimed it, so there's no need to walk further up.
        if (curr_update_level > level) {
          _internal::write_min(&curr->update_level_, level);

        }
        top_nodes[i] = nullptr;
        break;
      }
    }
  });

  parlay::parallel_for(0, len, [&](const size_t i) {
    if (top_nodes[i] != nullptr) {
      top_nodes[i]->UpdateTopDown(top_nodes[i]->height_ - 1);
    }
  });
}

void AugmentedElement::BatchJoin(
    const parlay::sequence<std::pair<AugmentedElement*, AugmentedElement*>>& joins)
{
  const size_t len{joins.size()};
  auto join_lefts{parlay::sequence<AugmentedElement*>::uninitialized(len)};
  parlay::parallel_for(0, len, [&](const size_t i) {
    Join(joins[i].first, joins[i].second);
    join_lefts[i] = joins[i].first;
  });
  BatchUpdate(join_lefts, {});
}

void AugmentedElement::BatchSplit(const parlay::sequence<AugmentedElement*>& splits) {
  const size_t len{splits.size()};
  parlay::parallel_for(0, len, [&](const size_t i) {
    splits[i]->Split();
  });
  parlay::parallel_for(0, len, [&](const size_t i) {
    AugmentedElement* curr{splits[i]};
    // `can_proceed` breaks ties when there are duplicate splits. When two
    // splits occur at the same place, only one of them should walk up and
    // update.
    bool can_proceed{
        curr->update_level_ == _internal::NA && CAS(&curr->update_level_, _internal::NA, 0)};
    if (can_proceed) {
      // Update values of `curr`'s ancestors.
      int sum{curr->values_[0]};
      int level{0};
      while (true) {
        if (level < curr->height_ - 1) {
          level++;
          curr->values_[level] = sum;
        } else {
          curr = curr->neighbors_[level].prev;
          if (curr == nullptr) {
            break;
          } else {
            sum += curr->values_[level];
          }
        }
      }
    }
  });
  parlay::parallel_for(0, len, [&](const size_t i) {
    splits[i]->update_level_ = _internal::NA;
  });
}

int AugmentedElement::GetSubsequenceSum(
    const AugmentedElement* left, const AugmentedElement* right) {
  int level{0};
  int sum{right->values_[level]};
  while (left != right) {
    level = std::min(left->height_, right->height_) - 1;
    if (level == left->height_ - 1) {
      sum += left->values_[level];
      left = left->neighbors_[level].next;
    } else {
      right = right->neighbors_[level].prev;
      sum += right->values_[level];
    }
  }
  return sum;
}

int AugmentedElement::GetSum() const {
  // Here we use knowledge of the implementation of `FindRepresentative()`.
  // `FindRepresentative()` gives some element that reaches the top level of the
  // list. For acyclic lists, the element is the leftmost one.
  AugmentedElement* root{FindRepresentative()};
  // Sum the values across the top level of the list.
  int level{root->height_ - 1};
  int sum{root->values_[level]};
  AugmentedElement* curr{root->neighbors_[level].next};
  while (curr != nullptr && curr != root) {
    sum += curr->values_[level];
    curr = curr->neighbors_[level].next;
  }
  if (curr == nullptr) {
    // The list is not circular, so we need to traverse backwards to beginning
    // of list and sum values to the left of `root`.
    curr = root;
    while (true) {
      while (level >= 0 && curr->neighbors_[level].prev == nullptr) {
        level--;
      }
      if (level < 0) {
        break;
      }
      while (curr->neighbors_[level].prev != nullptr) {
        curr = curr->neighbors_[level].prev;
        sum += curr->values_[level];
      }
    }
  }
  return sum;
}

}  // namespace parallel_skip_list

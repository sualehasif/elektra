#pragma once

#include "parallel_skip_list/skip_list.h"

#include <algorithm>

#include <parlay/monoid.h>
#include <parlay/parallel.h>
#include <parlay/sequence.h>
#include <parlay/utilities.h>

namespace parallel_skip_list {

// Batch-parallel augmented skip list. Each element is assigned a value, and
// partial sums over connected elements can be computed.
//
// A minimal example of using this class is given in the class AugmentedElement.
//
// Template arguments:
// - Derived: When using this class or creating a subclass of it, `Derived`
//   should be the (sub)class itself.
// - Func: This is the augmentation function and should be a class containing a
//   type `T` and an associative, commutative function `f` with type signature
//   `static T f(T a, T b)`.
template <typename Derived, typename Func>
class AugmentedElementBase : public ElementBase<Derived> {
 public:
  using Element = AugmentedElementBase<Derived, Func>;
  using Base = ElementBase<Derived>;
  using Value = typename Func::T;

  // Same as constructors of `ElementBase<>` except that there's an argument for
  // an initial `value` for the element.
  AugmentedElementBase(Value value);
  AugmentedElementBase(size_t random_int, Value value);
  ~AugmentedElementBase();

  // Can run concurrently with other `JoinWithoutUpdate` calls, but the augmented
  // value of `left` must be updated separately afterwards.
  static void JoinWithoutUpdate(Derived* left, Derived* right);
  // Can run concurrently with other `SplitWithoutUpdate` calls, but the augmented
  // value of the element must be updated separately afterwards.
  Derived* SplitWithoutUpdate();

  // For each `{left, right}` in the `len`-length array `joins`, concatenate the
  // list that `left` lives in to the list that `right` lives in.
  //
  // `left` must be the last node in its list, and `right` must be the first
  // node of in its list. Each `left` must be unique, and each `right` must be
  // unique.
  static void BatchJoin(const parlay::sequence<std::pair<Derived*, Derived*>>& joins);

  // For each `v` in the `len`-length array `splits`, split `v`'s list right
  // after `v`.
  static void BatchSplit(const parlay::sequence<Derived*>& splits);

  // For each `i`=0,1,...,`len`-1, assign value `new_values[i]` to element
  // `elements[i]`.
  static void BatchUpdate(
      const parlay::sequence<Derived*>& elements,
      const parlay::sequence<Value>& new_values);
  // Updates augmented values of the elements' ancestors according to whatever
  // values already exist at elements[]->values_[0]. This is useful after calls
  // to JoinWithoutUpdate() and SplitWithoutUpdate().
  //
  // The type `Seq` should behave like parlay::sequence<Derived*>.
  //
  // (The "ancestors" of an element e refers e->FindLeftParent(0),
  // e->FindLeftParent(0)->FindLeftParent(1),
  // e->FindLeftParent(0)->FindLeftParent(1)->FindLeftParent(2)`, and so on.)
  template <typename Seq>
  static void BatchUpdate(const Seq& elements);

  // Get the result of applying the augmentation function over the subsequence
  // between `left` and `right` inclusive.
  //
  // `left` and `right` must live in the same list, and `left` must precede
  // `right` in the list.
  //
  // This function does not modify the data structure, so it may run
  // concurrently with other `GetSubsequenceSum` calls and const function calls.
  static Value GetSubsequenceSum(const Derived* left, const Derived* right);

  // Get result of applying the augmentation function over the whole list that
  // the element lives in.
  Value GetSum() const;

  using Base::FindRepresentative;
  using Base::GetPreviousElement;
  using Base::GetNextElement;

 protected:
  static Value* AllocateValues(int height, Value default_value);

  static void UpdateTopDownImpl(int level, Derived* curr, bool is_loop_start = true);
  // Update aggregate value of node and clear `join_update_level` after joins.
  void UpdateTopDown(int level);
  void UpdateTopDownSequential(int level);

  static concurrent_array_allocator::Allocator<Value> values_allocator_;

  Value* values_;
  // When updating augmented values, this marks the lowest index at which the
  // `values_` needs to be updated.
  int update_level_;

  friend Base;
  using Base::Join;
  using Base::Split;
  using Base::height_;
  using Base::neighbors_;
};

// Basic augmented skip list augmented such that calling `elem.GetSum()` on
// element `elem` returns the size of the list containing `elem`. See interface of
// `AugmentedElementBase`.
class AugmentedElement : public AugmentedElementBase<AugmentedElement, parlay::addm<int>> {
  using Base = AugmentedElementBase<AugmentedElement, parlay::addm<int>>;

 public:
  AugmentedElement() : Base(1) {}
  AugmentedElement(size_t random_int) : Base(random_int, 1) {}

 private:
  friend Base;
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

template <typename D, typename F>
concurrent_array_allocator::Allocator<typename F::T> AugmentedElementBase<D, F>::values_allocator_{};

template <typename D, typename F>
typename F::T* AugmentedElementBase<D, F>::AllocateValues(int height, typename F::T default_value) {
  typename F::T* values{values_allocator_.Allocate(height)};
  for (int i = 0; i < height; i++) {
    values[i] = default_value;
  }
  return values;
}

template <typename D, typename F>
AugmentedElementBase<D, F>::AugmentedElementBase(typename F::T value) :
    ElementBase<D>{}, update_level_{_internal::NA} {
  values_ = AllocateValues(height_, value);
}

template <typename D, typename F>
AugmentedElementBase<D, F>::AugmentedElementBase(size_t random_int, typename F::T value) :
    ElementBase<D>{random_int}, update_level_{_internal::NA} {
  values_ = AllocateValues(height_, value);
}

template <typename D, typename F>
AugmentedElementBase<D, F>::~AugmentedElementBase() {
  values_allocator_.Free(values_, height_);
}

template <typename D, typename F>
void AugmentedElementBase<D, F>::JoinWithoutUpdate(D* left, D* right) {
  Join(left, right);
}

template <typename D, typename F>
D* AugmentedElementBase<D, F>::SplitWithoutUpdate() {
  return Split();
}

template <typename D, typename F>
void AugmentedElementBase<D, F>::UpdateTopDownSequential(int level) {
  if (level == 0) {
    if (height_ == 1) {
      update_level_ = _internal::NA;
    }
    return;
  }

  if (update_level_ < level) {
    UpdateTopDownSequential(level - 1);
  }
  typename F::T sum{values_[level - 1]};
  AugmentedElementBase* curr{neighbors_[level - 1].next};
  while (curr != nullptr && curr->height_ < level + 1) {
    if (curr->update_level_ != _internal::NA && curr->update_level_ < level) {
      curr->UpdateTopDownSequential(level - 1);
    }
    sum = F::f(sum, curr->values_[level - 1]);
    curr = curr->neighbors_[level - 1].next;
  }
  values_[level] = sum;

  if (height_ == level + 1) {
    update_level_ = _internal::NA;
  }
}

// Helper function for UpdateTopDown() that updates the augmented values of the
// descendants of `curr`.
template <typename D, typename F>
void AugmentedElementBase<D, F>::UpdateTopDownImpl(int level, D* curr, bool is_loop_start) {
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
template <typename D, typename F>
void AugmentedElementBase<D, F>::UpdateTopDown(int level) {
  if (level <= 6) {
    UpdateTopDownSequential(level);
    return;
  }

  UpdateTopDownImpl(level, static_cast<D*>(this));

  // Now that children have correct augmented valeus, update self's augmented
  // value.
  typename F::T sum{values_[level - 1]};
  AugmentedElementBase<D, F>* curr = neighbors_[level - 1].next;
  while (curr != nullptr && curr->height_ < level + 1) {
    sum = F::f(sum, curr->values_[level - 1]);
    curr = curr->neighbors_[level - 1].next;
  }
  values_[level] = sum;

  if (height_ == level + 1) {
    update_level_ = _internal::NA;
  }
}

template <typename D, typename F>
void AugmentedElementBase<D, F>::BatchUpdate(
    const parlay::sequence<D*>& elements,
    const parlay::sequence<typename F::T>& new_values) {
  parlay::parallel_for(0, elements.size(), [&](size_t i) {
    if (elements[i] != nullptr) {
      elements[i]->values_[0] = new_values[i];
    }
  });
  AugmentedElementBase<D, F>::BatchUpdate(elements);
}

template <typename D, typename F>
template <typename Seq>
void AugmentedElementBase<D, F>::BatchUpdate(const Seq& elements) {
  const size_t len{elements.size()};
  // The nodes whose augmented values need updating are the ancestors of
  // `elements`. Some nodes may share ancestors. `top_nodes` will contain,
  // without duplicates, the set of all ancestors of `elements` with no left
  // parents. From there we can walk down from those ancestors to update all
  // required augmented values.
  auto top_nodes{parlay::sequence<AugmentedElementBase<D, F>*>::uninitialized(len)};

  parlay::parallel_for(0, len, [&](const size_t i) {
    if (elements[i] == nullptr) {
      top_nodes[i] = nullptr;
      return;
    }

    int level{0};
    AugmentedElementBase<D, F>* curr{elements[i]};
    while (true) {
      int curr_update_level{curr->update_level_};
      if (curr_update_level == _internal::NA && CAS(&curr->update_level_, _internal::NA, level)) {
        level = curr->height_ - 1;
        AugmentedElementBase* parent{curr->FindLeftParent(level)};
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

template <typename D, typename F>
void AugmentedElementBase<D, F>::BatchJoin(const parlay::sequence<std::pair<D*, D*>>& joins)
{
  const size_t len{joins.size()};
  auto join_lefts{parlay::sequence<D*>::uninitialized(len)};
  parlay::parallel_for(0, len, [&](const size_t i) {
    Join(joins[i].first, joins[i].second);
    join_lefts[i] = joins[i].first;
  });
  BatchUpdate(join_lefts);
}

template <typename D, typename F>
void AugmentedElementBase<D, F>::BatchSplit(const parlay::sequence<D*>& splits) {
  const size_t len{splits.size()};
  parlay::parallel_for(0, len, [&](const size_t i) {
    splits[i]->Split();
  });
  parlay::parallel_for(0, len, [&](const size_t i) {
    AugmentedElementBase<D, F>* curr{splits[i]};
    // `can_proceed` breaks ties when there are duplicate splits. When two
    // splits occur at the same place, only one of them should walk up and
    // update.
    bool can_proceed{
        curr->update_level_ == _internal::NA && CAS(&curr->update_level_, _internal::NA, 0)};
    if (can_proceed) {
      // Update values of `curr`'s ancestors.
      typename F::T sum{curr->values_[0]};
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
            sum = F::f(sum, curr->values_[level]);
          }
        }
      }
    }
  });
  parlay::parallel_for(0, len, [&](const size_t i) {
    splits[i]->update_level_ = _internal::NA;
  });
}

template <typename D, typename F>
typename F::T AugmentedElementBase<D, F>::GetSubsequenceSum(const D* left, const D* right) {
  int level{0};
  typename F::T sum{right->values_[level]};
  while (left != right) {
    level = std::min(left->height_, right->height_) - 1;
    if (level == left->height_ - 1) {
      sum = F::f(sum, left->values_[level]);
      left = left->neighbors_[level].next;
    } else {
      right = right->neighbors_[level].prev;
      sum = F::f(sum, right->values_[level]);
    }
  }
  return sum;
}

template <typename D, typename F>
typename F::T AugmentedElementBase<D, F>::GetSum() const {
  // Here we use knowledge of the implementation of `FindRepresentative()`.
  // `FindRepresentative()` gives some element that reaches the top level of the
  // list. For acyclic lists, the element is the leftmost one.
  AugmentedElementBase<D, F>* root{FindRepresentative()};
  // Sum the values across the top level of the list.
  int level{root->height_ - 1};
  typename F::T sum{root->values_[level]};
  AugmentedElementBase<D, F>* curr{root->neighbors_[level].next};
  while (curr != nullptr && curr != root) {
    sum = F::f(sum, curr->values_[level]);
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
        sum = F::f(sum, curr->values_[level]);
      }
    }
  }
  return sum;
}

}  // namespace parallel_skip_list

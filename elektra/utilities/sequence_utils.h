#pragma once

#include <parlay/primitives.h>

namespace elektra {

// Given a sorted sequence, returns the indices of the first
// appearance of each unique element in order.
template <typename IndexType = size_t, PARLAY_RANGE_TYPE R, typename BinaryPredicate>
parlay::sequence<IndexType> get_offsets(const R& range, BinaryPredicate eq) {
  const auto unique_flags = parlay::delayed_seq<bool>(
      range.size(),
      [&](const size_t i) {
        return i == 0 || !eq(range[i], range[i - 1]);
      }
  );
  return parlay::pack_index<IndexType>(unique_flags);
}

template <typename IndexType = size_t, PARLAY_RANGE_TYPE R>
parlay::sequence<IndexType> get_offsets(const R& range) {
  return get_offsets<IndexType>(range, std::equal_to<parlay::range_value_type_t<R>>());
}

}  // namespace elektra

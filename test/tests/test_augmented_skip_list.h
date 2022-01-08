#pragma once

#include <parlay/sequence.h>

#include <algorithm>
#include <vector>

#include "../../elektra/parallel_skip_list/augmented_skip_list.h"
#include "gtest/gtest.h"

namespace elektra::testing {
namespace augmented_skip_list {

using Element = parallel_skip_list::AugmentedElement;

// Create `n` skip list elements.
parlay::sequence<Element> CreateElements(size_t n, size_t random_seed = 0) {
  auto elements_seq = parlay::sequence<Element>::uninitialized(n);
  Element* elements = elements_seq.data();
  parlay::random r{random_seed};
  parlay::parallel_for(
      0, n, [&](const size_t i) { new (&elements[i]) Element(r.ith_rand(i)); });
  return elements_seq;
}

// Convenience function for generating a sequence of joins given indices into
// `elements`.
parlay::sequence<std::pair<Element*, Element*>> CreateJoins(
    parlay::sequence<Element>& elements,
    const std::vector<std::vector<size_t>>& join_indices) {
  auto joins = parlay::sequence<std::pair<Element*, Element*>>::uninitialized(
      join_indices.size());
  parlay::parallel_for(0, join_indices.size(), [&](const size_t i) {
    joins[i] = std::make_pair(&elements[join_indices[i][0]],
                              &elements[join_indices[i][1]]);
  });
  return joins;
}

// Convenience function for generating a sequence of splits given indices into
// `elements`.
parlay::sequence<Element*> CreateSplits(
    parlay::sequence<Element>& elements,
    const std::vector<size_t>& split_indices) {
  auto splits = parlay::sequence<Element*>::uninitialized(split_indices.size());
  parlay::parallel_for(0, split_indices.size(), [&](const size_t i) {
    splits[i] = &elements[split_indices[i]];
  });
  return splits;
}

// Check that lists of size 1 have the correct sum.
TEST(AugmentedParallelSkipListTest, SumLengthOneLists) {
  auto elements = CreateElements(3);
  EXPECT_EQ(elements[0].GetSum(), 1);
  EXPECT_EQ(elements[1].GetSum(), 1);
  EXPECT_EQ(elements[2].GetSum(), 1);
}

// Check that lists have the correct sum after batch joins and splits to form
// linear lists.
TEST(AugmentedParallelSkipListTest, SumLinearList) {
  auto elements = CreateElements(8);
  parlay::sequence<std::pair<Element*, Element*>> joins =
      CreateJoins(elements, {
                                {0, 1},
                                {1, 2},
                                {3, 4},
                                {6, 7},
                            });
  Element::BatchJoin(joins);
  EXPECT_EQ(elements[0].GetSum(), 3);
  EXPECT_EQ(elements[2].GetSum(), 3);
  EXPECT_EQ(elements[3].GetSum(), 2);
  EXPECT_EQ(elements[4].GetSum(), 2);
  EXPECT_EQ(elements[5].GetSum(), 1);
  EXPECT_EQ(elements[6].GetSum(), 2);

  joins = CreateJoins(elements, {
                                    {2, 3},
                                    {4, 5},
                                    {5, 6},
                                });
  Element::BatchJoin(joins);
  EXPECT_EQ(elements[0].GetSum(), 8);
  EXPECT_EQ(elements[7].GetSum(), 8);
  EXPECT_EQ(Element::GetSubsequenceSum(&elements[0], &elements[7]), 8);
  EXPECT_EQ(Element::GetSubsequenceSum(&elements[2], &elements[5]), 4);
  EXPECT_EQ(Element::GetSubsequenceSum(&elements[1], &elements[2]), 2);
  EXPECT_EQ(Element::GetSubsequenceSum(&elements[6], &elements[6]), 1);

  parlay::sequence<Element*> splits = CreateSplits(elements, {0, 3, 4});
  Element::BatchSplit(splits);
  EXPECT_EQ(elements[0].GetSum(), 1);
  EXPECT_EQ(elements[3].GetSum(), 3);
  EXPECT_EQ(elements[4].GetSum(), 1);
  EXPECT_EQ(elements[6].GetSum(), 3);
}

// Check that lists have the correct sum after joins and splits to form circular
// lists.
TEST(AugmentedParallelSkipListTest, SumCircularList) {
  auto elements = CreateElements(6);
  parlay::sequence<std::pair<Element*, Element*>> joins =
      CreateJoins(elements, {
                                {0, 1},
                                {1, 2},
                                {2, 3},
                                {3, 4},
                                {4, 5},
                                {5, 0},
                            });
  Element::BatchJoin(joins);
  EXPECT_EQ(elements[4].GetSum(), 6);
  EXPECT_EQ(Element::GetSubsequenceSum(&elements[3], &elements[2]), 6);
  EXPECT_EQ(Element::GetSubsequenceSum(&elements[5], &elements[1]), 3);

  parlay::sequence<Element*> splits = CreateSplits(elements, {1, 2, 5});
  Element::BatchSplit(splits);
  EXPECT_EQ(elements[0].GetSum(), 2);
  EXPECT_EQ(elements[2].GetSum(), 1);
  EXPECT_EQ(elements[3].GetSum(), 3);
}

}  // namespace augmented_skip_list
}  // namespace elektra::testing

#pragma once

#include "../../elektra/parallel_skip_list/augmented_skip_list.h"

#include "gtest/gtest.h"

namespace elektra::testing {
namespace augmented_skip_list {

using Element = parallel_skip_list::AugmentedElement;

class AugmentedParallelSkipListTest : public ::testing::Test {
 protected:
  AugmentedParallelSkipListTest() {
    Element::Initialize();
  }

  ~AugmentedParallelSkipListTest() override {
    Element::Finish();
  }
};

// Create `n` skip list elements.
parlay::sequence<Element> CreateElements(size_t n, size_t random_seed = 0) {
  auto elements_seq = parlay::sequence<Element>::uninitialized(n);
  Element* elements = elements_seq.data();
  parlay::random r{random_seed};
  parlay::parallel_for(0, n, [&](int i) {
    new (&elements[i]) Element(r.ith_rand(i));
  });
  return elements_seq;
}

TEST_F(AugmentedParallelSkipListTest, IsDisconnectedInitially) {
  auto elements = CreateElements(3);
  EXPECT_NE(elements[0].FindRepresentative(), elements[1].FindRepresentative());
  EXPECT_NE(elements[0].FindRepresentative(), elements[2].FindRepresentative());
  EXPECT_NE(elements[1].FindRepresentative(), elements[2].FindRepresentative());
}

}  // namespace augmented_skip_list
}  // namespace elektra::testing

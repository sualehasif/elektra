#pragma once

#include "../../elektra/parallel_skip_list/skip_list.h"

#include "gtest/gtest.h"
#include <parlay/parallel.h>
#include <parlay/random.h>
#include <parlay/sequence.h>

namespace elektra::testing {
namespace {

using Element = parallel_skip_list::Element;

class ParallelSkipListTest : public ::testing::Test {
 protected:
  ParallelSkipListTest() {
    Element::Initialize();
  }

  ~ParallelSkipListTest() override {
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

TEST_F(ParallelSkipListTest, DisjointElements) {
  auto elements = CreateElements(2);
  EXPECT_NE(elements[0].FindRepresentative(), elements[1].FindRepresentative());
}

// Join all elements together in a long chain.
TEST_F(ParallelSkipListTest, JoinIntoChain) {
  const int n = 5;
  auto elements = CreateElements(n);

  parlay::parallel_for(0, n - 1, [&](int i) {
    Element::Join(&elements[i], &elements[i + 1]);
  });

  Element* representative_0{elements[0].FindRepresentative()};
  for (int i = 1; i < n; i++) {
    EXPECT_EQ(representative_0, elements[i].FindRepresentative());
  }
}

// Join all elements together in a cycle.
TEST_F(ParallelSkipListTest, JoinIntoCycle) {
  const int n = 5;
  auto elements = CreateElements(n);

  parlay::parallel_for(0, n, [&](int i) {
    Element::Join(&elements[i], &elements[(i + 1) % n]);
  });

  Element* representative_0{elements[0].FindRepresentative()};
  for (int i = 1; i < n; i++) {
    EXPECT_EQ(representative_0, elements[i].FindRepresentative());
  }
}

// Join all elements together in a cycle and then split everything back into
// individual elements.
TEST_F(ParallelSkipListTest, SplitIntoLengthOneLists) {
  const int n = 10;
  auto elements = CreateElements(n);

  parlay::parallel_for(0, n, [&](int i) {
    Element::Join(&elements[i], &elements[(i + 1) % n]);
  });
  parlay::parallel_for(0, n, [&](int i) {
    elements[i].Split();
  });

  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      EXPECT_NE(elements[i].FindRepresentative(), elements[j].FindRepresentative());
    }
  }
}

// Join all elements together in a cycle and then split everything into length-2
// lists.
TEST_F(ParallelSkipListTest, SplitIntoLengthTwoLists) {
  const int n = 10;
  ASSERT_EQ(n % 2, 0);
  auto elements = CreateElements(n);

  parlay::parallel_for(0, n, [&](int i) {
    Element::Join(&elements[i], &elements[(i + 1) % n]);
  });
  parlay::parallel_for(0, n / 2, [&](int i) {
    elements[2 * i + 1].Split();
  });

  for (int i = 0; i < n; i += 2) {
    EXPECT_EQ(elements[i].FindRepresentative(), elements[i + 1].FindRepresentative());
    for (int j = i + 2; j < n; j += 2) {
      EXPECT_NE(elements[i].FindRepresentative(), elements[j].FindRepresentative());
    }
  }
}

}  // namespace
}  // namespace elektra::testing

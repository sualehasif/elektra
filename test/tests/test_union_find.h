#include <parlay/parallel.h>
#include <parlay/random.h>
#include <parlay/sequence.h>

#include <random>

#include "../../elektra/union_find.h"
#include "gtest/gtest.h"

namespace elektra::testing {
namespace {

class UnionFindTest : public ::testing::Test {
 protected:
  UnionFindTest() {
    std::cout << "Testing on " << parlay::num_workers() << " Parlay workers."
              << std::endl;
    rng.seed(0);
  }

  std::mt19937 rng{};
};

// Test that the union-find algorithm works on a small set of elements.
TEST_F(UnionFindTest, SmallSet) {
  constexpr int num_elements = 10;
  constexpr auto find{find_variants::find_compress};
  auto unite{unite_variants::Unite<decltype(find)>{find}};

  parlay::sequence<uintE> elements{num_elements};

  parlay::parallel_for(0, num_elements, [&](int i) { elements[i] = i; });

  unite(0, 1, elements);
  unite(0, 4, elements);
  unite(0, 5, elements);

  unite(1, 6, elements);
  unite(1, 7, elements);

  unite(2, 8, elements);
  unite(2, 9, elements);

  ASSERT_EQ(find(0, elements), 0);
  ASSERT_EQ(find(1, elements), 0);
  ASSERT_EQ(find(2, elements), 2);
  ASSERT_EQ(find(3, elements), 3);
  ASSERT_EQ(find(4, elements), 0);
  ASSERT_EQ(find(5, elements), 0);
  ASSERT_EQ(find(6, elements), 0);
  ASSERT_EQ(find(7, elements), 0);
  ASSERT_EQ(find(8, elements), 2);
  ASSERT_EQ(find(9, elements), 2);
}

// Test that the union-find algorithm works on a large set of elements.
TEST_F(UnionFindTest, LargeSet) {
  constexpr int num_elements = 1000;
  constexpr auto find{find_variants::find_compress};
  auto unite{unite_variants::Unite<decltype(find)>{find}};

  parlay::sequence<uintE> elements{num_elements};

  parlay::parallel_for(0, num_elements, [&](int i) { elements[i] = i; });

  // for(int i = 0; i < num_elements; ++i) {
  parlay::parallel_for(0, num_elements, [&](int i) {
    for (int j = i + 1; j < num_elements; ++j) {
      if (i != j) {
        unite(i, j, elements);
      }
    }
  });

  // parlay::parallel_for(int i = 0; i < num_elements; ++i) {
  parlay::parallel_for(0, num_elements,
                       [&](int i) { ASSERT_EQ(find(i, elements), 0); });
}

}  // namespace
}  // namespace elektra::testing
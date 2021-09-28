#include <parlay/parallel.h>
#include <parlay/random.h>
#include <parlay/sequence.h>

#include <string>

#include "../../elektra/parallel_skip_list/skip_list.h"
#include "gtest/gtest.h"

// TODO: parallelize the loops in this based on the original.

namespace elektra::testing {
namespace {
using Element = parallel_skip_list::Element;
#define kNumElements 100

// The fixture for testing class Foo.
class ParallelSkipListTest : public ::testing::Test {
 protected:
  ParallelSkipListTest() {
    // print starting the test
    std::cout << "ParallelSkipListTest: Number of Elements: " << kNumElements
              << " elements" << std::endl;

    Element::Initialize();
    parlay::random r;

    elements_seq = parlay::sequence<Element>::uninitialized(kNumElements);
    elements = elements_seq.data();

    // parlay::parallel_for(0, kNumElements, [&](int i) {
    //   new (&elements[i]) Element(r.ith_rand(i));
    // });

    for (int i = 0; i < kNumElements; i++) {
      auto rand = r.ith_rand(i);
      new (&elements[i]) Element(rand);
    }
    PrimeSieve();
  }

  ~ParallelSkipListTest() override { Element::Finish(); }

  // Class members declared here can be used by all tests
  parlay::sequence<Element> elements_seq{};
  Element* elements;
  // const int kNumElements{1000};

  bool split_points[kNumElements];
  int start_index_of_list[kNumElements];

  // Mark that we will split at prime-numbered indices.
  inline void PrimeSieve() {
    split_points[2] = true;
    for (int i = 3; i < kNumElements; i += 2) split_points[i] = true;
    for (int64_t i = 3; i * i < kNumElements; i += 2) {
      if (split_points[i]) {
        for (int64_t j = i * i; j < kNumElements; j += 2 * i) {
          split_points[j] = false;
        }
      }
    }
  }
};

// Tests
// TEST_F(ParallelSkipListTest, ReperesentativeTest1_DisjointElements) {
//   int start_index{0};
//   for (int i = 0; i < kNumElements; i++) {
//     start_index_of_list[i] = start_index;
//     if (split_points[i]) {
//       start_index = i + 1;
//     }
//   }
//   start_index_of_list[0] = start_index_of_list[1] = start_index_of_list[2] =
//       start_index % kNumElements;

//   for (int i = 0; i < kNumElements; i++) {
//     Element* representative_i{elements[i].FindRepresentative()};

//     for (int j = i + 1; j < kNumElements; j++) {
//       std::cout << "here 4" << std::endl;
//       Element* representative_j{elements[j].FindRepresentative()};

//       bool assertion = (representative_j != representative_i);
//       std::cout << representative_j << " " << representative_i << std::endl;
//       std::cout << assertion << std::endl;
//       EXPECT_TRUE(assertion);
//       std::cout << assertion << std::endl;
//     }
//   }

//   // parlay::parallel_for(0, kNumElements, [&](int i) {
//   //   Element* representative_i{elements[i].FindRepresentative()};

//   //   for (int j = i + 1; j < kNumElements; j++) {
//   //     ASSERT_TRUE(representative_i != elements[j].FindRepresentative());
//   //   }
//   // });
// }

TEST_F(ParallelSkipListTest, ReperesentativeTest2_JoiningAllElements) {
  int start_index{0};
  for (int i = 0; i < kNumElements; i++) {
    start_index_of_list[i] = start_index;
    if (split_points[i]) {
      start_index = i + 1;
    }
  }
  start_index_of_list[0] = start_index_of_list[1] = start_index_of_list[2] =
      start_index % kNumElements;

  // could be parallelized
  // Join all elements together

  for (int i = 0; i < kNumElements - 1; i++) {
    Element::Join(&elements[i], &elements[i + 1]);
  }

  // could be parallelized
  Element* representative_0{elements[0].FindRepresentative()};
  for (int i = 0; i < kNumElements; i++) {
    assert(representative_0 == elements[i].FindRepresentative());
  }
}

// TEST_F(ParallelSkipListTest, ReperesentativeTest3_JoiningAllIntoOneBigCycle)
// {
//   int start_index{0};
//   // could be parallelized
//   for (int i = 0; i < kNumElements; i++) {
//     start_index_of_list[i] = start_index;
//     if (split_points[i]) {
//       start_index = i + 1;
//     }
//   }
//   start_index_of_list[0] = start_index_of_list[1] = start_index_of_list[2] =
//       start_index % kNumElements;

//   // Join all elements together
//   for (int i = 0; i < kNumElements - 1; i++) {
//     Element::Join(&elements[i], &elements[i + 1]);
//   }

//   // Join into one big cycle
//   Element::Join(&elements[kNumElements - 1], &elements[0]);

//   // could be parallelized
//   Element* representative_0 = elements[0].FindRepresentative();
//   for (int i = 0; i < kNumElements; i++) {
//     ASSERT_TRUE(representative_0 == elements[i].FindRepresentative());
//   }
// }

// TEST_F(ParallelSkipListTest,
//        ReperesentativeTest4_JoinIndividualListsIntoIndividualCycles) {
//   int start_index{0};
//   for (int i = 0; i < kNumElements; i++) {
//     start_index_of_list[i] = start_index;
//     if (split_points[i]) {
//       start_index = i + 1;
//     }
//   }
//   start_index_of_list[0] = start_index_of_list[1] = start_index_of_list[2] =
//       start_index % kNumElements;

//   // could be parallelized
//   // Join all elements together
//   for (int i = 0; i < kNumElements - 1; i++) {
//     Element::Join(&elements[i], &elements[i + 1]);
//   }

//   // Join into one big cycle
//   Element::Join(&elements[kNumElements - 1], &elements[0]);

//   // could be parallelized
//   // Join individual lists into individual cycles
//   for (int i = 0; i < kNumElements; i++) {
//     if (split_points[i]) {
//       Element::Join(&elements[i], &elements[start_index_of_list[i]]);
//     }
//   }

//   // could be parallelized
//   for (int i = 0; i < kNumElements; i++) {
//     const int start{start_index_of_list[i]};
//     ASSERT_TRUE(elements[start].FindRepresentative() ==
//                 elements[i].FindRepresentative());
//     if (start > 0) {
//       ASSERT_TRUE(elements[start - 1].FindRepresentative() !=
//                   elements[i].FindRepresentative());
//     }
//   }
// }

// TEST_F(ParallelSkipListTest,
//        ReperesentativeTest5_BreakCyclesIntoIndividualLists) {
//   int start_index{0};
//   for (int i = 0; i < kNumElements; i++) {
//     start_index_of_list[i] = start_index;
//     if (split_points[i]) {
//       start_index = i + 1;
//     }
//   }
//   start_index_of_list[0] = start_index_of_list[1] = start_index_of_list[2] =
//       start_index % kNumElements;

//   // could be parallelized
//   // Join all elements together
//   for (int i = 0; i < kNumElements - 1; i++) {
//     Element::Join(&elements[i], &elements[i + 1]);
//   }

//   // Join into one big cycle
//   Element::Join(&elements[kNumElements - 1], &elements[0]);

//   // could be parallelized
//   // Join individual lists into individual cycles
//   for (int i = 0; i < kNumElements; i++) {
//     if (split_points[i]) {
//       Element::Join(&elements[i], &elements[start_index_of_list[i]]);
//     }
//   }

//   // could be parallelized
//   // Break cycles back into lists
//   for (int i = 0; i < kNumElements; i++) {
//     if (split_points[i]) {
//       elements[i].Split();
//     }
//   }

//   // could be parallelized
//   // Join lists together
//   for (int i = 0; i < kNumElements; i++) {
//     if (split_points[i]) {
//       Element::Join(&elements[i], &elements[(i + 1) % kNumElements]);
//     }
//   }

//   Element* representative_0 = elements[0].FindRepresentative();
//   // could be parallelized
//   for (int i = 0; i < kNumElements; i++) {
//     ASSERT_TRUE(representative_0 == elements[i].FindRepresentative());
//   }
// }

}  // namespace
}  // namespace elektra::testing

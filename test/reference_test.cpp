#include "gtest/gtest.h"
#include "reference_tests/reference_connectivity.hpp"

auto main(int argc, char **argv) -> int {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

#include <iostream>
#include <random>

#include "benchmarks/benchmark.h"
#include "elektra/batch_dynamic_connectivity/dynamic_connectivity.h"
#include "elektra/utilities/simple_forest_connectivity.h"

using BatchDynamicConnectivity =
    bdcty::BatchDynamicConnectivity;
auto main(int argc, char** argv) -> int {
  std::cout << "Hello, from elektra!\n\n";
  elektra::benchmark::RunBenchmark<BatchDynamicConnectivity>(
      argc, argv, "Dynamic Connectivity");
}

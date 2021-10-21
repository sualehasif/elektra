#include <iostream>
#include <random>

#include "benchmarks/benchmark.h"
#include "elektra/batch_dynamic_connectivity/dynamic_connectivity.h"
#include "elektra/utilities/simple_forest_connectivity.h"

using BatchDynamicConnectivity =
    batchDynamicConnectivity::BatchDynamicConnectivity;
int main(int argc, char** argv) {
  std::cout << "Hello, from elektra!\n";
  elektra::benchmark::RunBenchmark<BatchDynamicConnectivity>(
      argc, argv, "Dynamic Connectivity");
}

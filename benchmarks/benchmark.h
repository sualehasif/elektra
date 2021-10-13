#pragma once

#include "get_time.h"
#include "io.h"
#include "parse_command_line.h"

namespace elektra {
namespace benchmark {

typedef std::pair<uintE, uintE> intPair;
typedef parlay::sequence<intPair> edgeList;
using timer = elektra::timer;

// For `num_iters` iterations, construct a forest from the `m` edges in `edges`,
// then batch cut and batch link the first `batch_size` edges in `edges`.
// Report the median batch cut and batch link time.
template <typename Connectivity>
void IncrementallUpdateConnectivity(edgeList edges, int batch_size,
                                    int num_iters, int n, int m) {
  vector<double> link_times(num_iters);

  Connectivity connect = Connectivity(n);

  std::random_device rd;
  std::mt19937 gen(rd());  // Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> distrib(1, m - batch_size);

  for (int j = 0; j < num_iters; j++) {
    timer link_t;

    // construct a slice of edges to insert
    auto start_slice = distrib(gen);
    edgeList edges_to_link{
        edges.slice(start_slice, start_slice + batch_size - 1)};

    // batch link and time
    link_t.start();
    forest->BatchAddEdges(edges_to_link);
    link_times[j] = link_t.stop();
  }
  const string batch_str{to_string(batch_size)};
  timer::report_time("link-" + batch_str, median(link_times));
}

template <typename Connectivity>
inline void RunBenchmark(int argc, char** argv) {
  commandLine P{argc, argv, "[-iters] graph_filename"};
  char* graph_filename{P.getArgument(0)};

  int num_iters{P.getOptionIntValue("-iters", 4)};
  int nworkers{P.getOptionIntValue("-workers", 1)};

  std::cout << "Running with " << nworkers << " workers" << std::endl;
  auto vertex_and_edge_list = io::read_unweighted_edge_list(graph_filename);

  const int n = vertex_and_edge_list.first;
  const int m = vertex_and_edge_list.size() - 1;

  sequence<intPair> edges = vertex_and_edge_list.second;

  std::random_device rd;
  std::mt19937 gen(rd());  // Standard mersenne_twister_engine seeded with rd()

  // shuffle the edges
  std::shuffle(edges.begin(), edges.end(), generator);

  //   Connectivity* connect = new Connectivity(n);

  for (int batch_size = 100; batch_size < m; batch_size *= 10) {
    IncrementallUpdateConnectivity<Connectivity>(edges, batch_size, num_iters,
                                                 n, m);
  }
  //   UpdateForest(&forest, edges, m, num_iters, m);
}

}  // namespace benchmark
}  // namespace elektra
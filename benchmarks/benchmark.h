#pragma once

#include <cstdlib>
#include <vector>

#include "../elektra/batch_dynamic_connectivity/connectivity-helpers.h"
#include "get_time.h"
#include "io.h"
#include "parse_command_line.h"

namespace elektra::benchmark {

using V = uint32_t;
using IntPair = std::pair<uintE, uintE>;
using std::vector;
using EdgeList = parlay::sequence<IntPair>;
using Timer = elektra::timer;
using E = std::pair<V, V>;

// TODO: move this to a utility file
template <typename T> auto median(std::vector<T> v) -> T {
  const size_t len = v.size();
  if (len == 0) {
    std::cerr << "median(): empty vector" << std::endl;
    return 0;
  }
  std::sort(v.begin(), v.end());
  if (len % 2) {
    return v[len / 2];
  }

  return v[len / 2 - 1] + (v[len / 2] - v[len / 2 - 1]) / 2;
}

auto IntPairBatchToEdgeArray(parlay::sequence<IntPair> &se)
    -> parlay::sequence<E> {
  // turns a sequence of edges to an array of pairs
  // useful for interfacing with EulerTourTrees
  parlay::sequence<E> array;

  for (auto &i : se) {
    array.push_back({i.first, i.second});
  }
  return array;
}

// For `num_iters` iterations, construct a forest from the `m` edges in `edges`,
// then batch cut and batch link the first `batch_size` edges in `edges`.
// Report the median batch cut and batch link time.
template <typename Connectivity>
void incrementallUpdateConnectivity(EdgeList edges, int batch_size,
                                    int num_iters, int n, int m) {
  vector<double> link_times(num_iters);

  Connectivity connect = Connectivity(n);

  std::random_device rd;
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> distrib(1, m - batch_size);

  for (int j = 0; j < num_iters; j++) {
    // print the iteration number
    // std::cout << "iteration " << j << std::endl;

    Timer link_t;

    // construct a slice of edges to insert (because Connectivity doesnt support
    // slices yet.)
    auto start_slice = distrib(gen);
    auto slice = edges.cut(start_slice, start_slice + batch_size - 1);
    EdgeList el{slice.begin(), slice.end()};

    // print the edges to be inserted
    // std::cout << "inserting edges: " << std::endl;
    // for (auto e : el) {
    //   std::cout << e.first << " " << e.second << std::endl;
    // }

    auto edges_to_link = IntPairBatchToEdgeArray(el);

    // batch link and time
    link_t.start();
    connect.BatchAddEdges(edges_to_link);

    link_times[j] = link_t.stop();
  }
  const std::string batch_str{std::to_string(batch_size)};
  Timer::report_time("  link-" + batch_str, median(link_times));
}

// [Experiment: Overhead in an Insertion-Only Scenario]
// Take an input graph
// Construct the graph
// Insert p% of the edges in the start
// Then incrementally add batches.
// (^ batch size varied, measurements done)
// Measure the throughput=edges/second as a function of batch size. Can plot
// incremental-only algorithm (UF) against our dynamic connectivity
// implementations
template <typename Connectivity>
void insertionOnly(EdgeList edges, int batch_size, int n, int m) {
  // construct a forest from the edges
  Connectivity connect = Connectivity(n);

  // insert p% of the edges in the start
  const int p = 40;
  const int p_edges = m * p / 100;
  auto start_slice = edges.cut(0, p_edges - 1);
  EdgeList el{start_slice.begin(), start_slice.end()};
  connect.BatchAddEdges(IntPairBatchToEdgeArray(el));

  // then incrementally add batches

  auto start_slice_idx = p_edges;

  const int num_iters = (m - p_edges) / batch_size - 2;

  vector<double> link_times(num_iters);

  for (int i = 0; i < num_iters; i++) {
    Timer link_t;

    auto el_batch =
        parlay::make_slice(edges.begin() + start_slice_idx,
                           edges.begin() + start_slice_idx + batch_size - 1);

    auto edges_to_link =
        parlay::sequence<E>::from_function(batch_size, [&](size_t i) {
          return E(el_batch[i].first, el_batch[i].second);
        });

    // batch link and time
    link_t.start();
    connect.BatchAddEdges(edges_to_link);

    link_times[i] = link_t.stop();

    start_slice_idx += batch_size;
  }
  const std::string kBatchStr{"  link-insertion-only" +
                              std::to_string(batch_size)};
  Timer::report_time(kBatchStr, median(link_times));
}

template <typename Connectivity>
inline void RunBenchmark(int argc, char **argv, const std::string &name) {
  commandLine P{argc, argv, "[-iters] [-workers] graph_filename"};
  char *graph_filename{P.getArgument(0)};

  int num_iters{P.getOptionIntValue("-iters", 5)};
  int nworkers{P.getOptionIntValue("-workers", 1)};

  // We start a benchmark
  std::cout << " ---- Starting benchmark ---- " << std::endl;
  std::cout << "  Benchmarks Name: " << name << std::endl;
  std::cout << "  num_iters: " << num_iters << std::endl;
  std::cout << "  Running with " << nworkers << " workers" << std::endl;
  std::cout << "  graph_filename: " << graph_filename << std::endl;

  auto vertex_and_edge_list = io::read_unweighted_edge_list(graph_filename);

  const int n = vertex_and_edge_list.first;

  sequence<IntPair> edges = vertex_and_edge_list.second;
  const int m = edges.size();

  // print the number of vertices and edges
  std::cout << "  n: " << n << std::endl;
  std::cout << "  m: " << m << std::endl;

  std::random_device rd;
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

  // shuffle the edges
  std::shuffle(edges.begin(), edges.end(), gen);

  //   Connectivity* connect = new Connectivity(n);

  // for (int batch_size = 2; batch_size < m; batch_size *= 10) {
  //   std::cout << "batch_size: " << batch_size << std::endl;
  //   incrementallUpdateConnectivity<Connectivity>(edges, batch_size,
  //   num_iters,
  //                                                n, m);
  // }
  insertionOnly<Connectivity>(edges, 100, n, m);
}

} // namespace elektra::benchmark
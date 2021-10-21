#pragma once

#include <vector>

#include "../elektra/batch_dynamic_connectivity/graph.h"
#include "get_time.h"
#include "io.h"
#include "parse_command_line.h"

namespace elektra {
namespace benchmark {

typedef std::pair<uintE, uintE> intPair;
using std::vector;
typedef parlay::sequence<intPair> edgeList;
using timer = elektra::timer;
using UndirectedEdge = dynamicGraph::UndirectedEdge;

// TODO: move this to a utility file
template <typename T>
T median(std::vector<T> v) {
  const size_t len = v.size();
  if (len == 0) {
    std::cerr << "median(): empty vector" << std::endl;
    abort();
  }
  std::sort(v.begin(), v.end());
  if (len % 2) {
    return v[len / 2];
  } else {
    return v[len / 2 - 1] + (v[len / 2] - v[len / 2 - 1]) / 2;
  }
}

parlay::sequence<UndirectedEdge> intPairBatchToEdgeArray(
    parlay::sequence<intPair>& se) {
  // turns a sequence of edges to an array of pairs
  // useful for interfacing with EulerTourTrees
  parlay::sequence<UndirectedEdge> array;

  for (size_t i = 0; i < se.size(); i++) {
    array.push_back(UndirectedEdge(se[i].first, se[i].second));
  }
  return array;
}

// For `num_iters` iterations, construct a forest from the `m` edges in `edges`,
// then batch cut and batch link the first `batch_size` edges in `edges`.
// Report the median batch cut and batch link time.
template <typename Connectivity>
void incrementallUpdateConnectivity(edgeList edges, int batch_size,
                                    int num_iters, int n, int m) {
  vector<double> link_times(num_iters);

  Connectivity connect = Connectivity(n);

  std::random_device rd;
  std::mt19937 gen(rd());  // Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> distrib(1, m - batch_size);

  for (int j = 0; j < num_iters; j++) {
    // print the iteration number
    std::cout << "iteration " << j << std::endl;

    timer link_t;

    // construct a slice of edges to insert (because Connectivity doesnt support
    // slices yet.)
    auto start_slice = distrib(gen);
    auto slice = edges.cut(start_slice, start_slice + batch_size - 1);
    edgeList el{slice.begin(), slice.end()};

    // print the edges to be inserted
    std::cout << "inserting edges: " << std::endl;
    for (auto e : el) {
      std::cout << e.first << " " << e.second << std::endl;
    }

    auto edges_to_link = intPairBatchToEdgeArray(el);

    std::cout << "here" << std::endl;

    // batch link and time
    link_t.start();
    connect.BatchAddEdges(edges_to_link);

    link_times[j] = link_t.stop();
  }
  const std::string batch_str{std::to_string(batch_size)};
  timer::report_time("link-" + batch_str, median(link_times));
}

template <typename Connectivity>
inline void RunBenchmark(int argc, char** argv, std::string name) {
  commandLine P{argc, argv, "[-iters] [-workers] graph_filename"};
  char* graph_filename{P.getArgument(0)};

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

  sequence<intPair> edges = vertex_and_edge_list.second;
  const int m = edges.size();

  // print the number of vertices and edges
  std::cout << "  n: " << n << std::endl;
  std::cout << "  m: " << m << std::endl;

  std::random_device rd;
  std::mt19937 gen(rd());  // Standard mersenne_twister_engine seeded with rd()

  // shuffle the edges
  std::shuffle(edges.begin(), edges.end(), gen);

  //   Connectivity* connect = new Connectivity(n);

  for (int batch_size = 2; batch_size < m; batch_size *= 10) {
    incrementallUpdateConnectivity<Connectivity>(edges, batch_size, num_iters,
                                                 n, m);
  }
  //   UpdateForest(&forest, edges, m, num_iters, m);
}

}  // namespace benchmark
}  // namespace elektra
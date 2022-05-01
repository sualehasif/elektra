#include <parlay/parallel.h>
#include <parlay/random.h>
#include <parlay/sequence.h>

#include <dynamic_graph/dynamic_connectivity.hpp>
#include <random>
#include <string>

#include "../../benchmarks/io.h"
#include "../../elektra/batch_dynamic_connectivity/dynamic_connectivity.h"
#include "gtest/gtest.h"

auto get_benchmark_path(const std::string& name) -> std::string {
  // get path using "git rev-parse --show-toplevel"
  char cwd[PATH_MAX];
  getcwd(cwd, sizeof(cwd));

  // find where "elektra" lies
  std::string path(cwd);
  size_t pos = path.find("elektra");
  path = path.substr(0, pos) + "elektra/";

  return path + "benchmarks/graphs/" + name;
}

namespace elektra::testing {
namespace {
using std::make_pair;

using E = pair<V, V>;
using BatchDynamicConnectivity = bdcty::BatchDynamicConnectivity;

class ReferenceConnectivityTest : public ::testing::Test {
 protected:
  ReferenceConnectivityTest() {
    std::cout << "Testing on " << parlay::num_workers() << " Parlay workers."
              << std::endl;
    rng.seed(0);
  }

  std::mt19937 rng{};
};

TEST_F(ReferenceConnectivityTest, SimpleGraph) {
  const auto filename = get_benchmark_path("email-Eu-core.txt");
  auto vertex_and_edge_list = io::read_unweighted_edge_list(filename.c_str());

  const int num_vertices = vertex_and_edge_list.first;
  const auto edges_base = vertex_and_edge_list.second;
  const int num_edges = edges_base.size();

  std::cout << "num_vertices: " << num_vertices << std::endl;
  std::cout << "num_edges: " << num_edges << std::endl;

  DynamicConnectivity reference(num_vertices);
  BatchDynamicConnectivity batch_connect(num_vertices);

  // remove self-loops
  auto edges = parlay::sequence<E>::from_function(num_edges, [&](size_t i) {
    if (edges_base[i].first < edges_base[i].second) {
      return E(edges_base[i].first, edges_base[i].second);
    }
    return E(UINT32_MAX - 1, UINT32_MAX - 1);
  });
  edges = parlay::filter(edges,
                         [&](const E& e) { return e.first != UINT32_MAX - 1; });

  for (const auto& edge : edges) {
    const auto e = UndirectedEdgeR(edge.first, edge.second);
    reference.AddEdge(e);
  }

  batch_connect.BatchAddEdges(edges_base);

  // // generate random queries
  // std::uniform_int_distribution<int> dist(0, num_vertices - 1);
  // const int num_queries = 100;
  // auto queries = sequence<E>::from_function(
  //     num_queries, [&]() { return make_pair(dist(rng), dist(rng)); });

  // // run queries
  // auto batch_results = batch_connect.BatchConnected(queries);
  // for (int i = 0; i < num_queries; i++) {
  //   EXPECT_EQ(batch_results[i], static_cast<char>(reference.IsConnected(
  //                                   queries[i].first, queries[i].second)));
  // }
}

}  // namespace
}  // namespace elektra::testing
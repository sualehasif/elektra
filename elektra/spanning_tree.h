#pragma once

#include <parlay/sequence.h>

#include <unordered_set>

#include "hash_pair.h"
#include "macros.h"
#include "union_find.h"

namespace elektra {
using std::make_pair;

template<typename T>
class SpanningTree {
 private:
  using E = std::pair<T, T>;
  using Ehash = HashIntPairStruct;
  static constexpr auto empty_edge_ = make_pair(UINT_E_MAX, UINT_E_MAX);

  parlay::sequence<uintE> parents_;

 public:
  // return the tree
  auto get_spanning_tree();

  SpanningTree(const parlay::sequence<E> &edges, parlay::sequence<E> &tree) {
    constexpr auto find{find_variants::find_compress};
    auto unite{unite_variants::Unite<decltype(find)>{find}};

    // count the number of vertices
    auto left_vertices = parlay::sequence<uintE>::from_function(
        edges.size(), [&](uintE i) { return edges[i].first; });
    auto right_vertices = parlay::sequence<uintE>::from_function(
        edges.size(), [&](uintE i) { return edges[i].second; });

    // sort the vertices
    parlay::integer_sort_inplace(left_vertices);
    parlay::integer_sort_inplace(right_vertices);

    // merge the vertices
    auto vertices =
        parlay::unique(parlay::merge(left_vertices, right_vertices));

    // create a map from the vertices to their index in the sequence
    auto vertex_map =
        elektra::make_resizable_table<uintE, uintE, std::hash<uintE>>(
            vertices.size(), empty_edge_,
            make_tuple(UINT_E_MAX - 1, UINT_E_MAX - 1), std::hash<uintE>{});

    // fill the map
    parlay::parallel_for(0, vertices.size(), [&](uintE i) {
      vertex_map.insert(make_tuple(vertices[i], i));
    });

    // resize parents_ to the number of vertices
    // reserve enough space in parents_ for all elements
    auto parents = parlay::sequence<uintE>::from_function(
        vertices.size(), [](uintE i) { return i; });

    tree = parlay::sequence<E>::from_function(
        vertices.size(), [&](uintE i) { return empty_edge_; });

    // union the elements in the edges
    parlay::parallel_for(0, edges.size(), [&](uintE i) {
      auto e = edges[i];
      auto first = vertex_map.find(e.first);
      auto second = vertex_map.find(e.second);

      unite(first, second, parents, tree);
    });

    // map the tree to the original vertices
    parlay::parallel_for(0, tree.size(), [&](uintE i) {
      auto e = tree[i];
      auto first = ((e.first != UINT_E_MAX) ? vertices[e.first] : -1);
      auto second = ((e.first != UINT_E_MAX) ? vertices[e.second] : -1);

      tree[i] = make_pair(first, second);
    });

    tree = parlay::filter(tree, [&](const E &e) {
      return (e.first != static_cast<T>(-1)) && (e.second != static_cast<T>(-1));
    });
    // print the parents_ for each element in edges if you want
    // std::cout << "parents_: " << std::endl;
    // for (auto &e : edges) {
    //   auto first = vertex_map.find(e.first);
    //   auto second = vertex_map.find(e.second);

    //   std::cout << "( " << first << ", " << second << " ) "
    //             << ": "
    //             << "( " << parents_[first] << " " << parents_[second] << " ) "
    //             << std::endl;
    // }

    // print the tree if you want.
    // std::cout << "tree: " << std::endl;
    // for (auto &e : tree) {
    //   std::cout << e.first << " " << e.second << " -> "
    //             << ((e.first != UINT_E_MAX) ? vertices[e.first] : -1) << " "
    //             << ((e.second != UINT_E_MAX) ? vertices[e.second] : -1)
    //             << std::endl;
    // }
  }
};

} // namespace elektra

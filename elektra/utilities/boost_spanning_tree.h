//
// Created by Sualeh Asif on 12/25/21.
//

#ifndef ELEKTRA_ELEKTRA_UTILITIES_BOOST_SPANNING_TREE_H_
#define ELEKTRA_ELEKTRA_UTILITIES_BOOST_SPANNING_TREE_H_

#include <boost/pending/disjoint_sets.hpp>
#include <boost/property_map/property_map.hpp>
#include <parlay/primitives.h>
#include <parlay/sequence.h>
#include <utility>
#include <unordered_set>
#include <map>

#include "../hash_pair.h"

namespace elektra::deprecate {
using V = int;
using std::pair;
using E = pair<V, V>;
using EHash = HashIntPairStruct;

using TreeSet = std::unordered_set<E, EHash>;

template<typename Rank, typename Parent>
auto ConstructTree(Rank &r, Parent &p,
                   const parlay::sequence<E> &se)
-> TreeSet {
  // Given a sequence of edges, returns a set of a
  // spanning forest of the graph formed by them
  boost::disjoint_sets<Rank, Parent> dsu(r, p);
  TreeSet tree;

  // BUG POSSIBLE: check whether this causes a problem.
  for (auto v : se) {
    dsu.make_set(v.first);
    dsu.make_set(v.second);
  }

  for (auto v : se) {
    V first = v.first;
    V second = v.second;
    // TODO is there a race condition here if we parallelize this? How can we
    // resolve that

    if (dsu.find_set(first) != dsu.find_set(second)) {
      tree.insert(v);
      dsu.link(first, second);
    }
  }
  return tree;
}
auto GetSpanningTree(const parlay::sequence<E> &se)
-> TreeSet {
  // I am assuming the interface in
  // https://github.com/ParAlg/gbbs/blob/master/gbbs/union_find.h?fbclid=IwAR0U_Nbe1SpQF7mbmN0CEGLyF-5v362oy1q-9eQLvjQz916jhfTH69bMx9s
  // could be worth it to parallelize this

  using RankT = std::map<V, size_t>;
  using ParentT = std::map<V, V>;

  RankT rank_map;
  ParentT parent_map;

  boost::associative_property_map<RankT> rank_pmap(rank_map);
  boost::associative_property_map<ParentT> parent_pmap(parent_map);

  return ConstructTree(rank_pmap, parent_pmap, se);
}
} // namespace elektra::deprecate
#endif //ELEKTRA_ELEKTRA_UTILITIES_BOOST_SPANNING_TREE_H_

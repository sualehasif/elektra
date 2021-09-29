#include "connectivity.h"

namespace batchDynamicConnectivity {

// -------------
// PUBLIC METHODS
// -------------

// Batch Dynamic Connectivity: BatchConnected
// @input: A sequence of pairs of vertices
// @output: A sequence boolean values which are true if the corresponding
// vertices are connected
// @description: This method checks if the batch of vertices are connected in
// the graph
parlay::sequence<char> BatchDynamicConnectivity::BatchConnected(
    parlay::sequence<std::pair<Vertex, Vertex>> suv) const {
    parlay::sequence<char> s(suv.size(), 0);
  // check if they are connected in the highest level forest
  BatchDynamicET *pMaxLevelEulerTree =
      parallel_spanning_forests_[max_level_ - 1];

  parlay::parallel_for(0, suv.size(), [&](int i) {
    auto v1 = suv[i].first;
    auto v2 = suv[i].second;
    s[i] = pMaxLevelEulerTree->IsConnected(v1, v2);
  });

  return s;
}

// TODO: fix the getRepresentative function
// parlay::sequence<Vertex> BatchDynamicConnectivity::BatchFindRepr(
//     const parlay::sequence<Vertex> &sv) {
//   auto pMaxLevelEulerTree = parallel_spanning_forests_[max_level_ - 1];
//   // get representatives from the highest level foredt
//   return parlay::map(sv, [&](Vertex v) {
//     return (int64_t)pMaxLevelEulerTree->getRepresentative(v);
//   });
// }

// TODO: get the concurrent union find for this
template <typename Rank, typename Parent>
treeSet BatchDynamicConnectivity::constructTree(
    Rank &r, Parent &p, const parlay::sequence<UndirectedEdge> &se) {
  // Given a sequence of edges, returns a set of a
  // spanning forest of the graph formed by them
  boost::disjoint_sets<Rank, Parent> dsu(r, p);
  treeSet tree;

  // BUG POSSIBLE: check whether this causes a problem.
  for (auto v : se) {
    dsu.make_set(v.first);
    dsu.make_set(v.second);
  }

  for (auto v : se) {
    Vertex first = v.first;
    Vertex second = v.second;
    // TODO is there a race condition here if we paralize this? How can we
    // resolve that

    if (dsu.find_set(first) != dsu.find_set(second)) {
      tree.insert(v);
      dsu.link(first, second);
    }
  }
  return tree;
}

// TODO: add parallel DSU structure to implement this
treeSet BatchDynamicConnectivity::getSpanningTree(
    const parlay::sequence<UndirectedEdge> &se) {
  // I am assuming the interface in
  // https://github.com/ParAlg/gbbs/blob/master/gbbs/union_find.h?fbclid=IwAR0U_Nbe1SpQF7mbmN0CEGLyF-5v362oy1q-9eQLvjQz916jhfTH69bMx9s
  // could be worth paralelizing this

  typedef std::map<Vertex, size_t> rank_t;
  typedef std::map<Vertex, Vertex> parent_t;

  rank_t rank_map;
  parent_t parent_map;

  boost::associative_property_map<rank_t> rank_pmap(rank_map);
  boost::associative_property_map<parent_t> parent_pmap(parent_map);

  return constructTree(rank_pmap, parent_pmap, se);
}

parlay::sequence<std::pair<int, int>> edgeBatchToPairArray(
    parlay::sequence<UndirectedEdge> &se) {
  // turns a sequence of edges to an array of pairs
  // useful for interfacing with EulerTourTrees
  auto array = parlay::sequence<std::pair<int, int>>::uninitialized(se.size());

  parlay::parallel_for(0, se.size(), [&](int i) {
    array[i] = std::make_pair(se[i].first, se[i].second);
  });

  // parallel_for(int i = 0; i < se.size(); i++) {
  //   array[i].first = se[i].first;
  //   array[i].second = se[i].second;
  // }
  return array;
}

void BatchDynamicConnectivity::BatchAddEdges(
    const parlay::sequence<UndirectedEdge> &se) {
  auto maxLevelEulerTree = parallel_spanning_forests_[max_level_ - 1];

  parlay::sequence<UndirectedEdge> auxiliaryEdges =
      parlay::map(se, [&](UndirectedEdge e) {
        return UndirectedEdge(
            (Vertex)maxLevelEulerTree->getRepresentative(e.first),
            (Vertex)maxLevelEulerTree->getRepresentative(e.second));
      });
  auto tree = getSpanningTree(auxiliaryEdges);

  parlay::sequence<std::pair<int, int>> treeEdges;
  parlay::sequence<UndirectedEdge> nonTreeEdges;

  // update the tree and nonTree edges based on the ST computation
  parlay::parallel_for(0, se.size(), [&](int i) {
    if (tree.count(se[i])) {
      treeEdges.push_back(make_pair(se[i].first, se[i].second));
      detail::EdgeInfo ei = {(detail::Level)(max_level_ - 1),
                             detail::EdgeType::kTree};
      edges_[se[i]] = ei;
    } else {
      nonTreeEdges.push_back(se[i]);
      detail::EdgeInfo ei = {(detail::Level)(max_level_ - 1),
                             detail::EdgeType::kNonTree};
      edges_[se[i]] = ei;
    }
  });

  // parallel_for(int i = 0; i < se.size(); i++) {
  //   if (tree.count(se[i])) {
  //     treeEdges.push_back(se[i]);
  //     detail::EdgeInfo ei = {(detail::Level)(max_level_ - 1),
  //                            detail::EdgeType::kTree};
  //     edges_[se[i]] = ei;
  //   } else {
  //     nonTreeEdges.push_back(se[i]);
  //     detail::EdgeInfo ei = {(detail::Level)(max_level_ - 1),
  //                            detail::EdgeType::kNonTree};
  //     edges_[se[i]] = ei;
  //   }
  // }

  // add tree edges
  maxLevelEulerTree->BatchLink(treeEdges);

  // add to adjacancy list
  parlay::parallel_for(0, nonTreeEdges.size(), [&](int i) {
    non_tree_adjacency_lists_[max_level_ - 1][nonTreeEdges[i].first].insert(
        nonTreeEdges[i].second);
    non_tree_adjacency_lists_[max_level_ - 1][nonTreeEdges[i].second].insert(
        nonTreeEdges[i].first);
  });
}

// TODO implement semisort or any sort
auto BatchDynamicConnectivity::removeDuplicates(parlay::sequence<Vertex> &seq) {
  // TODO: possibly change this to use a not inplace sort
  // FIXME: Turn this into a integer sort by converting vertices to uints
  parlay::sort_inplace(seq);

  auto newSeq = parlay::unique(seq);
  return newSeq;
}

UndirectedEdge BatchDynamicConnectivity::componentSearch(int level, Vertex v) {
  auto levelEulerTree = parallel_spanning_forests_[level];
  // TODO: make sure there is a return in every case.
  auto cc = levelEulerTree->ConnectedComponent(v);
  for (int i = 0; i < cc.size(); i++) {
    auto u = cc[i];
    for (auto w : non_tree_adjacency_lists_[level][u]) {
      if (levelEulerTree->getRepresentative(u) !=
          levelEulerTree->getRepresentative(v)) {
        return UndirectedEdge(u, w);
      }
    }
  }

  // BUG: WE HAVE A BUG HERE
  assert(false);
  return UndirectedEdge(0, 0);
}

parlay::sequence<Vertex> BatchDynamicConnectivity::parallelLevelSearch(
    const parlay::sequence<UndirectedEdge> &se,
    parlay::sequence<Vertex> &components,
    parlay::sequence<std::pair<int, int>> &promotedEdges, int level) {
  auto levelEulerTree = parallel_spanning_forests_[level];
  levelEulerTree->BatchLink(promotedEdges);

  auto ncomponents = parlay::map(components, [&](Vertex v) {
    return (Vertex)levelEulerTree->getRepresentative(v);
  });
  components = ncomponents;
  components = removeDuplicates(components);

  parlay::sequence<Vertex> componentsToConsider;
  parlay::sequence<Vertex> largeComponents;

  parlay::parallel_for(0, components.size(), [&](int i) {
    if (levelEulerTree->ConnectedComponent(components[i]).size() <=
        1 << (level - 1)) {
      componentsToConsider.push_back(components[i]);
    } else {
      largeComponents.push_back(components[i]);
    }
  });
  // parallel_for(int i = 0; i < components.size(); i++) {
  //   if (levelEulerTree->ConnectedComponent(components[i]).size() <=
  //       1 << (level - 1)) {
  //     componentsToConsider.push_back(components[i]);
  //   } else {
  //     largeComponents.push_back(components[i]);
  //   }
  // }

  parlay::sequence<UndirectedEdge> R;

  while (componentsToConsider.size() != 0) {
    parlay::sequence<std::pair<int, int>> edgesToDropLevel;
    for (Vertex v : componentsToConsider) {
      // Move all the edges of small components down a level
      auto componentSearched = levelEulerTree->ConnectedComponent(v);

      // FIXME: convert to a parlay::sequence for now
      std::unordered_set<Vertex> cc;
      for (int i = 0; i < componentSearched.size(); ++i) {
        cc.insert(componentSearched[i]);
      }

      // BUG: fix the speed here
      for (const auto &ed : edges_) {  // edges loops over the edges of the tree
        auto e = ed.first;
        auto details = ed.second;
        if (details.level == level && details.type == detail::EdgeType::kTree &&
            (cc.find(e.first) != cc.end() || cc.find(e.second) != cc.end())) {
          edges_[e] = detail::EdgeInfo{
              static_cast<detail::Level>(details.level - 1), edges_[e].type};
          edgesToDropLevel.push_back(make_pair(e.first, e.second));
        }
      }
      R.push_back(componentSearch(level, v));
    }
    parallel_spanning_forests_[level - 1]->BatchLink(edgesToDropLevel);

    auto maxLevelEulerTree = parallel_spanning_forests_[max_level_];
    auto auxiliaryEdges = parlay::map(R, [&](UndirectedEdge e) {
      return UndirectedEdge(maxLevelEulerTree->getRepresentative(e.first),
                            maxLevelEulerTree->getRepresentative(e.second));
    });
    auto promIndices = getSpanningTree(auxiliaryEdges);

    parlay::sequence<std::pair<int, int>> newPromotedEdges;
    parlay::sequence<UndirectedEdge> notPromotedEdges;

    parlay::parallel_for(0, promIndices.size(), [&](int i) {
      if (promIndices.find(auxiliaryEdges[i]) != promIndices.end())
        newPromotedEdges.push_back(
            make_pair(auxiliaryEdges[i].first, auxiliaryEdges[i].second));
      else
        notPromotedEdges.push_back(auxiliaryEdges[i]);
    });

    // parallel_for(int i = 0; i < auxiliaryEdges.size(); i++) {
    //   if (promIndices.find(auxiliaryEdges[i]) != promIndices.end())
    //     newPromotedEdges.push_back(auxiliaryEdges[i]);
    //   else
    //     notPromotedEdges.push_back(auxiliaryEdges[i]);
    // }

    levelEulerTree->BatchLink(newPromotedEdges);

    parlay::parallel_for(0, newPromotedEdges.size(), [&](size_t i) {
      UndirectedEdge e = {newPromotedEdges[i].first,
                          newPromotedEdges[i].second};
      level = edges_[e].level;
      auto u = newPromotedEdges[i].first;
      auto v = newPromotedEdges[i].second;

      auto ul = non_tree_adjacency_lists_[level][u];
      auto vl = non_tree_adjacency_lists_[level][v];

      ul.erase(ul.find(u));
      vl.erase(ul.find(v));

      promotedEdges.push_back(newPromotedEdges[i]);
    });
    // parallel_for(int i = 0; i < newPromotedEdges.size(); i++) {
    //   level = edges_[newPromotedEdges[i]].level;
    //   auto u = newPromotedEdges[i].first;
    //   auto v = newPromotedEdges[i].second;

    //   auto ul = non_tree_adjacency_lists_[level][u];
    //   auto vl = non_tree_adjacency_lists_[level][v];

    //   ul.erase(ul.find(u));
    //   vl.erase(ul.find(v));

    //   promotedEdges.push_back(newPromotedEdges[i]);
    // }

    auto ccu = parlay::map(componentsToConsider, [&](Vertex v) {
      return (Vertex)levelEulerTree->getRepresentative(v);
    });

    componentsToConsider = removeDuplicates(ccu);

    parlay::sequence<Vertex> newComponentsToConsider;

    parlay::parallel_for(0, componentsToConsider.size(), [&](int i) {
      if (levelEulerTree->ConnectedComponent(componentsToConsider[i]).size() <=
          1 << (level - 1)) {
        newComponentsToConsider.push_back(componentsToConsider[i]);
      } else {
        largeComponents.push_back(componentsToConsider[i]);
      }
    });

    // parallel_for(int i = 0; i < componentsToConsider.size(); i++) {
    //   if (levelEulerTree->ConnectedComponent(componentsToConsider[i]).size()
    //   <=
    //       1 << (level - 1)) {
    //     newComponentsToConsider.push_back(componentsToConsider[i]);
    //   } else {
    //     largeComponents.push_back(componentsToConsider[i]);
    //   }
    // }
    componentsToConsider = newComponentsToConsider;
  }
  return largeComponents;
}

void BatchDynamicConnectivity::BatchDeleteEdges(
    const parlay::sequence<UndirectedEdge> &se) {
  // TODO: split se into tree and non tree edges
  // delete edges from adjacency list
  parlay::sequence<UndirectedEdge> treeEdges;

  auto min_tree_edge_level = max_level_;

  parlay::parallel_for(0, se.size(), [&](int i) {
    auto level = edges_[se[i]].level;
    auto u = se[i].first;
    auto v = se[i].second;

    auto ul = non_tree_adjacency_lists_[level][u];
    auto vl = non_tree_adjacency_lists_[level][v];

    if (ul.find(v) != ul.end()) {
      ul.erase(ul.find(v));
      vl.erase(vl.find(u));
    } else {
      treeEdges.push_back(se[i]);
      if (level < min_tree_edge_level) {
        min_tree_edge_level = level;
      }
    }
  });

  for (int level = min_tree_edge_level; level < max_level_; level++) {
    auto levelEulerTree = parallel_spanning_forests_[level];
    auto toDelete = parlay::filter(
        treeEdges, [&](UndirectedEdge e) { return edges_[e].level <= level; });

    parlay::sequence<std::pair<int, int>> toDeletePairSequence =
        edgeBatchToPairArray(toDelete);

    levelEulerTree->BatchCut(toDeletePairSequence);
  }

  parlay::sequence<Vertex> lcomponents =
      parlay::map(treeEdges, [](UndirectedEdge e) { return e.first; });
  parlay::sequence<Vertex> rcomponents =
      parlay::map(treeEdges, [](UndirectedEdge e) { return e.second; });
  lcomponents.append(rcomponents);

  auto components = removeDuplicates(lcomponents);

  parlay::sequence<std::pair<int, int>> promotedEdges;
  for (int i = min_tree_edge_level; i < max_level_; i++) {
    components = parallelLevelSearch(se, components, promotedEdges, i);
  }

  return;
}
}  // namespace batchDynamicConnectivity
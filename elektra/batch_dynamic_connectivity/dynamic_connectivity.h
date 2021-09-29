#include "connectivity.h"

namespace batchDynamicConnectivity {
template <typename T>
using sequence = parlay::sequence<T>;
using V = Vertex;

using std::pair;
// -------------
// PUBLIC METHODS
// -------------

// Batch Dynamic Connectivity: BatchConnected
// @input: A sequence of pairs of vertices
// @output: A sequence boolean values which are true if the corresponding
// vertices are connected
// @description: This method checks if the batch of vertices are connected in
// the graph
sequence<char> BatchDynamicConnectivity::BatchConnected(
    sequence<pair<V, V>> suv) const {
  sequence<char> s(suv.size(), 0);
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

// BUG (Possible): the tree count check needs to happen on auxillary edges.
// BUG (Possible): we need to check if the edges are already present in the
// graph.
void BatchDynamicConnectivity::BatchAddEdges(
    const sequence<UndirectedEdge> &se) {
  // Look at the max level Euler Tour Tree in the parallel spanning forests.
  auto maxLevelEulerTree = parallel_spanning_forests_[max_level_ - 1];

  sequence<UndirectedEdge> auxiliaryEdges =
      parlay::map(se, [&](UndirectedEdge e) {
        return UndirectedEdge(
            (V)maxLevelEulerTree->getRepresentative(e.first),
            (V)maxLevelEulerTree->getRepresentative(e.second));
      });
  auto tree = getSpanningTree(auxiliaryEdges);

  sequence<pair<int, int>> treeEdges;
  sequence<UndirectedEdge> nonTreeEdges;

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
auto BatchDynamicConnectivity::removeDuplicates(sequence<V> &seq) {
  // TODO: possibly change this to use a not inplace sort
  // FIXME: Turn this into a integer sort by converting vertices to uints
  parlay::sort_inplace(seq);

  auto newSeq = parlay::unique(seq);
  return newSeq;
}

UndirectedEdge BatchDynamicConnectivity::componentSearch(int level, V v) {
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

sequence<V> BatchDynamicConnectivity::parallelLevelSearch(
    const sequence<UndirectedEdge> &se, sequence<V> &components,
    sequence<pair<int, int>> &promotedEdges, int level) {
  auto levelEulerTree = parallel_spanning_forests_[level];
  levelEulerTree->BatchLink(promotedEdges);

  auto ncomponents = parlay::map(
      components, [&](V v) { return (V)levelEulerTree->getRepresentative(v); });
  components = ncomponents;
  components = removeDuplicates(components);

  sequence<V> componentsToConsider;
  sequence<V> largeComponents;

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

  sequence<UndirectedEdge> R;

  while (componentsToConsider.size() != 0) {
    sequence<pair<int, int>> edgesToDropLevel;
    for (V v : componentsToConsider) {
      // Move all the edges of small components down a level
      auto componentSearched = levelEulerTree->ConnectedComponent(v);

      // FIXME: convert to a sequence for now
      std::unordered_set<V> cc;
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

    sequence<pair<int, int>> newPromotedEdges;
    sequence<UndirectedEdge> notPromotedEdges;

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

    auto ccu = parlay::map(componentsToConsider, [&](V v) {
      return (V)levelEulerTree->getRepresentative(v);
    });

    componentsToConsider = removeDuplicates(ccu);

    sequence<V> newComponentsToConsider;

    parlay::parallel_for(0, componentsToConsider.size(), [&](int i) {
      if (levelEulerTree->ConnectedComponent(componentsToConsider[i]).size() <=
          1 << (level - 1)) {
        newComponentsToConsider.push_back(componentsToConsider[i]);
      } else {
        largeComponents.push_back(componentsToConsider[i]);
      }
    });

    // parallel_for(int i = 0; i < componentsToConsider.size(); i++) {
    //   if
    //   (levelEulerTree->ConnectedComponent(componentsToConsider[i]).size()
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
    const sequence<UndirectedEdge> &se) {
  // TODO: split se into tree and non tree edges
  // delete edges from adjacency list
  sequence<UndirectedEdge> treeEdges;

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

    sequence<pair<int, int>> toDeletePairSequence =
        edgeBatchToPairArray(toDelete);

    levelEulerTree->BatchCut(toDeletePairSequence);
  }

  sequence<V> lcomponents =
      parlay::map(treeEdges, [](UndirectedEdge e) { return e.first; });
  sequence<V> rcomponents =
      parlay::map(treeEdges, [](UndirectedEdge e) { return e.second; });
  lcomponents.append(rcomponents);

  auto components = removeDuplicates(lcomponents);

  sequence<pair<int, int>> promotedEdges;
  for (int i = min_tree_edge_level; i < max_level_; i++) {
    components = parallelLevelSearch(se, components, promotedEdges, i);
  }

  return;
}
}  // namespace batchDynamicConnectivity
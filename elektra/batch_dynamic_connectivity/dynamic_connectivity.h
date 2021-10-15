#include <parlay/primitives.h>

#include <atomic>
#include <mutex>

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

// Batch Dynamic Connectivity: BatchAddEdges
// @input: A sequence of pairs of vertices
// @output: void
// @description: Inserts the edges in the graph and updates the spanning
// forests(ETTs).
// BUG (Possible): we need to check if the edges are already
// present in the graph.
void BatchDynamicConnectivity::BatchAddEdges(
    const sequence<UndirectedEdge> &se) {
  // Look at the max level Euler Tour Tree in the parallel spanning forests.
  auto maxLevelEulerTree = parallel_spanning_forests_[max_level_ - 1];

  // you can print all the edges for debugging
  // std::cout << "Adding edges to the graph" << std::endl;
  // for (auto &e : se) {
  //   std::cout << "Edge: " << e.first << " " << e.second << std::endl;
  // }

  // Construct the auxillary edges.
  sequence<UndirectedEdge> auxiliaryEdges =
      parlay::map(se, [&](UndirectedEdge e) {
        return UndirectedEdge(
            (V)maxLevelEulerTree->getRepresentative(e.first),
            (V)maxLevelEulerTree->getRepresentative(e.second));
      });
  auto tree = getSpanningTree(auxiliaryEdges);

  // print all the tree edges for debugging if wanted.
  // std::cout << "Tree edges" << std::endl;
  // for (auto &e : tree) {
  //   std::cout << "Edge: " << e.first << " " << e.second << std::endl;
  // }

  auto num_edges_inserted = se.size();

  sequence<pair<int, int>> treeEdges;
  sequence<UndirectedEdge> nonTreeEdges;

  treeEdges.reserve(tree.size());
  nonTreeEdges.reserve(se.size() - tree.size());

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
// template <typename T>
auto BatchDynamicConnectivity::removeDuplicates(sequence<int> &seq) {
  // TODO: possibly change this to use a not inplace sort
  // FIXME: Turn this into a integer sort by converting vertices to uints
  parlay::integer_sort_inplace(seq, [](int x) { return (unsigned)x; });

  auto newSeq = parlay::unique(seq);
  return newSeq;
}

UndirectedEdge BatchDynamicConnectivity::componentSearch(int level, V v) {
  auto levelEulerTree = parallel_spanning_forests_[level];
  // TODO: make sure there is a return in every case.
  auto cc = levelEulerTree->ComponentEdges(v);
  // `ComponentEdges()` + `seen_vertices` is a silly way to iterate over
  // vertices of a component but we're going to refactor this anyway
  std::unordered_set<V> seen_vertices;
  for (int i = 0; i < (int)cc.size(); i++) {
    auto e = cc[i];
    for (auto u : {e.first, e.second}) {
      if (seen_vertices.count(u) > 0) {
        continue;
      }
      seen_vertices.insert(u);
      for (auto w : non_tree_adjacency_lists_[level][u]) {
        if (levelEulerTree->getRepresentative(u) !=
            levelEulerTree->getRepresentative(v)) {
          return UndirectedEdge(u, w);
        }
      }
    }
  }

  // BUG: WE HAVE A BUG HERE
  assert(false);
  return UndirectedEdge(0, 0);
}

// static inline auto getUnionFind() {}

void BatchDynamicConnectivity::replacementSearch(
    int level, sequence<int> components,
    sequence<pair<int, int>> &promoted_edges) {
  // set up a course search size and an overall stride for the component search
  int search_size = 256;
  int total_search_stride = 0;

  int critical_component_size = 1 << level;

  sequence<char> component_indicator(components.size(), 0);

  sequence<pair<int, int>> push_down_edges;

  auto ett = parallel_spanning_forests_[level];

  // TODO: move union find into it's own class.s
  // set up a union find data structure to keep track of the components
  typedef std::map<int, size_t> rank_t;
  typedef std::map<int, int> parent_t;

  rank_t rank_map;
  parent_t parent_map;

  typedef boost::associative_property_map<rank_t> Rank;
  typedef boost::associative_property_map<parent_t> Parent;

  Rank rank_pmap(rank_map);
  Parent parent_pmap(parent_map);
  boost::disjoint_sets<Rank, Parent> UF(rank_pmap, parent_pmap);

  // initialize the union find data structure
  for (int i = 0; i < (int)components.size(); i++) {
    UF.make_set(components[i]);
  }

  while (parlay::count(component_indicator, 1) < components.size()) {
    parlay::parallel_for(0, components.size(), [&](int i) {
      // check if the component is already marked
      if (component_indicator[i] == 1) return;

      auto c = components[i];
      auto cc = ett->ComponentEdges(c);
      const int ccSize = cc.size() + 1;

      // size check to ensure that we are only running over small components
      if (ccSize > critical_component_size) return;

      // doing a search over the components.
      parlay::parallel_for(0, search_size, [&](int j) {
        auto idx = j + total_search_stride;
        // if we have already searched this component, skip it.
        if (idx >= ccSize) {
          component_indicator[i] = 1;
          return;
        } else {
          auto e = UndirectedEdge(cc[idx]);
          if (edges_[e].type == detail::EdgeType::kNonTree) {
            auto u = e.first;
            auto v = e.second;

            auto u_rep = ett->getRepresentative(u);
            auto v_rep = ett->getRepresentative(v);

            auto u_parent = UF.find_set(u_rep);
            auto v_parent = UF.find_set(v_rep);

            if (u_parent != v_parent) {
              // link the two components together
              UF.link(u_parent, v_parent);
              // this should just be promoted.
              promoted_edges.push_back(make_pair(u, v));
            } else {
              // this failed and will be pushed down a level.
              push_down_edges.push_back(make_pair(u, v));
            }
          }
        }
      });

      // TODO: push down the edges and clear out the push_down_edges
      // Note we will promote the edges later.
    });

    // update the total search stride
    total_search_stride += search_size;

    // update the search size
    search_size *= 2;
  }
}

void BatchDynamicConnectivity::BatchDeleteEdges(
    const sequence<UndirectedEdge> &se) {
  // split se into tree and non tree edges
  sequence<UndirectedEdge> treeEdges;
  // reserve space for the sequence
  treeEdges.reserve(se.size());

  int8_t min_tree_edge_level = max_level_;
  // std::atomic<int8_t> min_tree_edge_level(max_level_);

  std::mutex m;
  std::unique_lock<std::mutex> lock(m, std::defer_lock);
  // delete edges from the non tree adjacency lists.
  parlay::parallel_for(0, se.size(), [&](int i) {
    auto level = edges_[se[i]].level;
    auto u = se[i].first;
    auto v = se[i].second;

    auto ul = non_tree_adjacency_lists_[level][u];
    auto vl = non_tree_adjacency_lists_[level][v];

    // if (u, v) is not a tree edge, then we need to remove it from the
    // non tree adjacency list
    // TODO: could this possibly do this with edges_[se[i]] == kTree?
    if (ul.find(v) != ul.end()) {
      ul.erase(ul.find(v));
      vl.erase(vl.find(u));
    } else {
      treeEdges.push_back(se[i]);

      // atomic compare and swap to find the minimum tree edge level
      // FIXME: Fix this. Find a way to do this without locking
      if (level < min_tree_edge_level) {
        lock.lock();
        if (level < min_tree_edge_level) {
          min_tree_edge_level = level;
        }
        lock.unlock();
      }
      // if (level < min_tree_edge_level) {
      //   std::atomic_compare_exchange_strong(
      //       &min_tree_edge_level, &min_tree_edge_level,
      //       std::min(atomic_level, min_tree_edge_level));
      // }
    }
  });

  // delete edges from the tree at each level from the minimum tree edge level
  // to the maximum tree edge level
  for (int l = min_tree_edge_level; l < max_level_; l++) {
    auto levelEulerTree = parallel_spanning_forests_[l];

    // get the edges to delete which have level at max l.
    auto toDelete = parlay::filter(
        treeEdges, [&](UndirectedEdge e) { return edges_[e].level <= l; });

    // FIXME: We do an extraneously expensive copy here because of API design.
    sequence<pair<int, int>> toDeletePairSequence =
        edgeBatchToPairArray(toDelete);

    levelEulerTree->BatchCut(toDeletePairSequence);
  }

  // We now try to find replacement edges starting from the min level
  // and working our way up.
  for (int l = min_tree_edge_level; l < max_level_; l++) {
    auto levelEulerTree = parallel_spanning_forests_[l];

    auto edgesToReplace = parlay::filter(
        treeEdges, [&](UndirectedEdge e) { return edges_[e].level == l; });

    sequence<int> lcomponents =
        parlay::map(edgesToReplace, [&](UndirectedEdge e) {
          return levelEulerTree->getRepresentative(e.first);
        });
    sequence<int> rcomponents =
        parlay::map(edgesToReplace, [&](UndirectedEdge e) {
          return levelEulerTree->getRepresentative(e.second);
        });
    lcomponents.append(rcomponents);

    auto components_to_consider = removeDuplicates(lcomponents);
    // FIXME: we might have to reserve space here.
    sequence<pair<int, int>> promoted_edges;

    // Doing a replacement search for the disconnected edge at level l.
    replacementSearch(l, components_to_consider, promoted_edges);

    levelEulerTree->BatchLink(promoted_edges);

    // TODO: possibly coursen the loop
    // update the tree edges with the promoted edges for all the higher levels
    parlay::parallel_for(l + 1, max_level_, [&](int i) {
      auto levelEulerTree = parallel_spanning_forests_[i];
      levelEulerTree->BatchLink(promoted_edges);
    });
  }

  return;
}

// sequence<V> BatchDynamicConnectivity::parallelLevelSearch(
//     const sequence<UndirectedEdge> &se, sequence<V> &components,
//     sequence<pair<int, int>> &promotedEdges, int level) {
//   auto levelEulerTree = parallel_spanning_forests_[level];
//   levelEulerTree->BatchLink(promotedEdges);

//   auto ncomponents = parlay::map(
//       components, [&](V v) { return (V)levelEulerTree->getRepresentative(v);
//       });
//   components = ncomponents;
//   // components = removeDuplicates(components);

//   sequence<V> components_to_consider;
//   sequence<V> largeComponents;

//   parlay::parallel_for(0, components.size(), [&](int i) {
//     if (levelEulerTree->ConnectedComponent(components[i]).size() <=
//         1 << (level - 1)) {
//       components_to_consider.push_back(components[i]);
//     } else {
//       largeComponents.push_back(components[i]);
//     }
//   });

//   sequence<UndirectedEdge> R;

//   while (components_to_consider.size() != 0) {
//     sequence<pair<int, int>> edgesToDropLevel;
//     for (V v : components_to_consider) {
//       // Move all the edges of small components down a level
//       auto componentSearched = levelEulerTree->ConnectedComponent(v);

//       // FIXME: convert to a sequence for now
//       std::unordered_set<V> cc;
//       for (int i = 0; i < (int)componentSearched.size(); ++i) {
//         cc.insert(componentSearched[i]);
//       }

//       // BUG: fix the speed here
//       for (const auto &ed : edges_) {  // edges loops over the edges of the
//       tree
//         auto e = ed.first;
//         auto details = ed.second;
//         if (details.level == level && details.type == detail::EdgeType::kTree
//         &&
//             (cc.find(e.first) != cc.end() || cc.find(e.second) != cc.end()))
//             {
//           edges_[e] = detail::EdgeInfo{
//               static_cast<detail::Level>(details.level - 1), edges_[e].type};
//           edgesToDropLevel.push_back(make_pair(e.first, e.second));
//         }
//       }
//       R.push_back(componentSearch(level, v));
//     }
//     parallel_spanning_forests_[level - 1]->BatchLink(edgesToDropLevel);

//     auto maxLevelEulerTree = parallel_spanning_forests_[max_level_];
//     auto auxiliaryEdges = parlay::map(R, [&](UndirectedEdge e) {
//       return UndirectedEdge(maxLevelEulerTree->getRepresentative(e.first),
//                             maxLevelEulerTree->getRepresentative(e.second));
//     });
//     auto promIndices = getSpanningTree(auxiliaryEdges);

//     sequence<pair<int, int>> newPromotedEdges;
//     sequence<UndirectedEdge> notPromotedEdges;

//     parlay::parallel_for(0, promIndices.size(), [&](int i) {
//       if (promIndices.find(auxiliaryEdges[i]) != promIndices.end())
//         newPromotedEdges.push_back(
//             make_pair(auxiliaryEdges[i].first, auxiliaryEdges[i].second));
//       else
//         notPromotedEdges.push_back(auxiliaryEdges[i]);
//     });

//     levelEulerTree->BatchLink(newPromotedEdges);

//     parlay::parallel_for(0, newPromotedEdges.size(), [&](size_t i) {
//       UndirectedEdge e = {newPromotedEdges[i].first,
//                           newPromotedEdges[i].second};
//       level = edges_[e].level;
//       auto u = newPromotedEdges[i].first;
//       auto v = newPromotedEdges[i].second;

//       auto ul = non_tree_adjacency_lists_[level][u];
//       auto vl = non_tree_adjacency_lists_[level][v];

//       ul.erase(ul.find(u));
//       vl.erase(ul.find(v));

//       promotedEdges.push_back(newPromotedEdges[i]);
//     });

//     auto ccu = parlay::map(components_to_consider, [&](V v) {
//       return (V)levelEulerTree->getRepresentative(v);
//     });

//     // components_to_consider = removeDuplicates(ccu);

//     sequence<V> newComponentsToConsider;

//     parlay::parallel_for(0, components_to_consider.size(), [&](int i) {
//       if (levelEulerTree->ConnectedComponent(components_to_consider[i])
//               .size() <= 1 << (level - 1)) {
//         newComponentsToConsider.push_back(components_to_consider[i]);
//       } else {
//         largeComponents.push_back(components_to_consider[i]);
//       }
//     });

//     components_to_consider = newComponentsToConsider;
//   }
//   return largeComponents;
// }
}  // namespace batchDynamicConnectivity

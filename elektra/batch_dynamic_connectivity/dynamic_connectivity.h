#include <parlay/primitives.h>

#include <atomic>
#include <mutex>
#include <set>

#include "connectivity.h"
namespace bdcty {
template <typename T> using sequence = parlay::sequence<T>;

// -------------
// PUBLIC METHODS
// -------------

// Batch Dynamic Connectivity: BatchConnected
// @input: A sequence of pairs of vertices
// @output: A sequence boolean values which are true if the corresponding
// vertices are connected
// @description: This method checks if the batch of vertices is connected in
// the graph
auto BatchDynamicConnectivity::BatchConnected(sequence<pair<V, V>> suv) const
    -> sequence<char> {
  sequence<char> s_connected(suv.size(), 0);
  // check if they are connected in the highest level forest
  BatchDynamicEtt *p_max_level_euler_tree =
      parallel_spanning_forests_[max_level_ - 1];

  parlay::parallel_for(0, suv.size(), [&](size_t i) {
    auto [v_1, v_2] = suv[i];
    s_connected[i] = p_max_level_euler_tree->IsConnected(v_1, v_2);
  });

  return s_connected;
}

// Batch Dynamic Connectivity: BatchAddEdges
// @input: A sequence of pairs of vertices
// @output: void
// @description: Inserts the edges in the graph and updates the spanning
// forests(ETTs).
// BUG (Possible): we need to check if the edges are already
// present in the graph.
void BatchDynamicConnectivity::BatchAddEdges(const sequence<E> &se) {
  // Look at the max level Euler Tour Tree in the parallel spanning forests.
  BatchDynamicEtt *max_level_euler_tree =
      parallel_spanning_forests_[max_level_ - 1];

#ifdef DEBUG
  // you can print all the edges for debugging
  std::cout << "Adding edges to the graph" << std::endl;
  for (auto &e : se) {
    std::cout << "Edge: " << e.first << " " << e.second << std::endl;
  }
#endif

  // Construct the auxiliary edges.
  sequence<E> auxiliary_edges = parlay::map(se, [&](E e) {
    return E(max_level_euler_tree->GetRepresentative(e.first),
             max_level_euler_tree->GetRepresentative(e.second));
  });
  // TODO(laxmand): change GetSpanningTree to use the concurrent
  // union-find code.
  // TODO(sualeh): Seems we can just make a TreeSet a parlay::sequence
  // of UndirectedEdges?
  auto tree = GetSpanningTree(auxiliary_edges);

  // sequence<pair<uintE, uintE>> aux_int_edges =
  //     parlay::map(se, [&](E e) {
  //       return
  //       make_pair((uintE)max_level_euler_tree->GetRepresentative(e.first),
  //                        (uintE)max_level_euler_tree->GetRepresentative(e.second));
  //     });

  // auto stree = sequence<pair<uintE,uintE>>::uninitialized (se.size());
  // elektra::SpanningTree(aux_int_edges, stree);

  // auto spanning_tree = stree.cut(1, stree.size());

// print all the tree edges for debugging if wanted.
#ifdef DEBUG
  std::cout << "Tree edges" << std::endl;
  for (auto &e : tree) {
    std::cout << "Edge: " << e.first << " " << e.second << std::endl;
  }
#endif

  auto tree_edges = sequence<pair<int, int>>();

  // FIXME(sualeh): This is a hack to get the tree edges.
  for (const auto &e : tree) {
    tree_edges.push_back(make_pair(e.first, e.second));
  }

#ifdef DEBUG
  printEdgeSequence(tree_edges, "Tree Edges from sequence");
#endif

  auto non_tree_edges =
      sequence<pair<int, int>>(se.size(), pair<int, int>(-1, -1));

  // This inserts the tree edges in to our hashtable.
  parlay::parallel_for(0, tree_edges.size(), [&](size_t i) {
    auto &e = tree_edges[i];
    bdcty::EInfo ei = {(bdcty::Level)(max_level_ - 1), bdcty::EType::K_TREE};
    edges_.insert(make_tuple(pair<V, V>(e.first, e.second), ei));

    // TODO(sualeh): Think about whether you can get away without
    // inserting the reverse edge add the reverse edge to the edges_ map
    bdcty::EInfo ei_rev = {(bdcty::Level)(max_level_ - 1),
                           bdcty::EType::K_TREE};
    edges_.insert(make_tuple(pair<V, V>(e.second, e.first), ei_rev));
  });

  // update the nonTree edges based on the ST computation
  parlay::parallel_for(0, se.size(), [&](size_t i) {
    if (tree.count(se[i]) == 0) {
      // std::string edge_str = "Non Tree Edge: " +
      // std::to_string(se[i].first)
      // +
      //                        " " + std::to_string(se[i].second) + "\n";
      // std::cout << edge_str;
      non_tree_edges[i] = (make_pair(se[i].first, se[i].second));
      bdcty::EInfo ei = {(bdcty::Level)(max_level_ - 1),
                         bdcty::EType::K_NON_TREE};
      edges_.insert(make_tuple(pair<V, V>(se[i].first, se[i].second), ei));

      // add the reverse edge to the edges_ map
      bdcty::EInfo ei_rev = {(bdcty::Level)(max_level_ - 1),
                             bdcty::EType::K_NON_TREE};
      edges_.insert(make_tuple(pair<V, V>(se[i].second, se[i].first), ei_rev));
    }
  });

  // add tree edges
  // std::cout << "Adding tree edges to the ETT" << std::endl;
  max_level_euler_tree->BatchLink(tree_edges);

  // filter and add non-tree edges
  auto filtered_non_tree_edges = parlay::filter(
      non_tree_edges, [&](auto e) { return e.first != -1 && e.second != -1; });

  // add to adjacency list
#ifdef DEBUG
  std::cout << "Adding the non-tree edges to the hashtable" << std::endl;
  std::cout << "non_tree_edges.size() = " << non_tree_edges.size() << std::endl;
  printEdgeSequence(non_tree_edges, "nonTreeEdgesUnfiltered");
  printEdgeSequence(filtered_non_tree_edges, "nonTreeEdgesFiltered");
#endif

  non_tree_adjacency_lists_.BatchAddEdgesToLevel(filtered_non_tree_edges,
                                                 max_level_ - 1);
}

// template <typename T>
auto BatchDynamicConnectivity::RemoveDuplicates(sequence<int> &seq)
    -> sequence<int> {
  // TODO: possibly change this to use a semisort and not inplace sort
  // sort the sequence
  parlay::integer_sort_inplace(seq, [](int x) { return (unsigned)x; });
  auto new_seq = parlay::unique(seq);

  return new_seq;
}

// static inline auto getUnionFind() {}

void BatchDynamicConnectivity::ReplacementSearch(
    int level, sequence<int> components,
    sequence<pair<int, int>> &promoted_edges) {
  // set up a course search size and an overall stride for the component
  // search
  int search_size = 256;
  int total_search_stride = 0;

#ifdef DEBUG
  // print out the components to consider
  std::cout << "Components to consider" << std::endl;
  for (auto &c : components) {
    std::cout << c << " ";
  }
  std::cout << std::endl;
#endif

  // TODO: I think this should be 1 << (level - 1)
  auto critical_component_size = 1 << level;

  // make a vector of the components to consider with size components.size()
  // and initialized to 0
  std::vector<char> component_indicator(components.size(), 0);

  sequence<pair<int, int>> push_down_edges;

  auto *ett = parallel_spanning_forests_[level];

  // TODO: move union find into it's own class.s
  // set up a union find data structure to keep track of the components
  using RankT = std::map<int, size_t>;
  using ParentT = std::map<int, int>;

  RankT rank_map;
  ParentT parent_map;

  using Rank = boost::associative_property_map<RankT>;
  using Parent = boost::associative_property_map<ParentT>;

  Rank rank_pmap(rank_map);
  Parent parent_pmap(parent_map);
  boost::disjoint_sets<Rank, Parent> uf(rank_pmap, parent_pmap);

  // initialize the union find data structure
  for (int component : components) {
    uf.make_set(component);
  }

  // TODO(tom): This needs to be supported as an augmentation.
  // set up non-tree edges to search through
  // The first sequence indexes by the component number
  // The second sequence indexes by the edges in the component
  vector<vector<E>> non_tree_edges;

  non_tree_edges.reserve((int)components.size());
  for (int i = 0; i < (int)components.size(); i++) {
    non_tree_edges.emplace_back();
  }

  for (int i = 0; i < (int)components.size(); i++) {
    auto comp_id = components[i];
    // get the ComponentVertices for the component
    auto cc = ett->ComponentVertices(comp_id);

    // for each vertex in the component
    for (int u : cc) {
      // for each edge in the component
      for (auto &kw : non_tree_adjacency_lists_[level][u].entries()) {
        auto w = std::get<0>(kw);
        // push the edge into the non_tree_edges vector
        non_tree_edges[i].push_back(E(u, w));
      }
    }
  }

  while (parlay::count(component_indicator, 1) < components.size()) {
    parlay::parallel_for(0, components.size(), [&](int i) {
      // check if the component is already marked
      if (component_indicator[i] == 1) {
        return;
      }

      auto non_tree_edges_i = non_tree_edges[i];
      const int kNumEdges = non_tree_edges_i.size();

      // size check to ensure that we are only running over small components
      if (kNumEdges > critical_component_size) {
        return;
      }

      // doing a search over the components.
      parlay::parallel_for(0, min(kNumEdges, search_size) + 1, [&](int j) {
        auto idx = j + total_search_stride;
        // if we have already searched this component, skip it.
        if (idx >= kNumEdges) {
          component_indicator[i] = 1;
          return;
        }
        auto e = E{non_tree_edges_i[idx].first, non_tree_edges_i[idx].second};
        auto u = e.first;
        auto v = e.second;

        auto u_rep = ett->GetRepresentative(u);
        auto v_rep = ett->GetRepresentative(v);

        auto u_parent = uf.find_set(u_rep);
        auto v_parent = uf.find_set(v_rep);

        if (u_parent != v_parent) {
          // link the two components together
          uf.link(u_parent, v_parent);
          // this should just be promoted.
          promoted_edges.push_back(make_pair(u, v));
        } else {
          // this failed and will be pushed down a level.
          push_down_edges.push_back(make_pair(u, v));
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

#ifdef DEBUG
  // print the promoted edges
  std::cout << "Promoted edges" << std::endl;
  for (auto &e : promoted_edges) {
    std::cout << e.first << " " << e.second << std::endl;
  }
  std::cout << std::endl;
#endif
}

void BatchDynamicConnectivity::BatchDeleteEdges(sequence<E> &se) {
  // first make sure all the edges are correctly oriented and remove all the
  // edges that are not in the graph

  parlay::parallel_for(0, (int)se.size(), [&](int i) {
    auto e = se[i];
    auto u = e.first;
    auto v = e.second;

    if (edges_.find(pair<V, V>(v, u)) != empty_info_) {
      // swap the edges
      se[i] = E(v, u);
    } else if (edges_.find(pair<V, V>(u, v)) == empty_info_) {
      // if the edge is not in the graph, skip it
      se[i] = E(-1, -1);
    }
  });

  // filter out the (-1, -1) edges
  parlay::filter(se, [&](E e) { return e.first != -1; });

  // split se into tree and non tree edges
  sequence<E> tree_edges;
  // reserve space for the sequence
  tree_edges.reserve(se.size());

  int min_tree_edge_level = static_cast<int>(max_level_);
  // std::atomic<int8_t> min_tree_edge_level(max_level_);

  std::mutex m;
  std::unique_lock<std::mutex> lock(m, std::defer_lock);

  // delete edges from the non tree adjacency lists.
  // TODO(sualeh): maybe collect the edges to be deleted and then delete them
  // in a batch
  parlay::parallel_for(0, se.size(), [&](int i) {
    auto level = edges_.find(pair<V, V>(se[i].first, se[i].second)).level;
    auto u = se[i].first;
    auto v = se[i].second;

    auto ul = non_tree_adjacency_lists_[level][u];
    auto vl = non_tree_adjacency_lists_[level][v];

#ifdef DEBUG
    // print ul
    std::cout << "ul" << std::endl;
    for (auto &e : ul.entries()) {
      std::cout << "(" << u << " " << std::get<0>(e) << "), ";
    }
    std::cout << std::endl;

    // print vl
    std::cout << "vl" << std::endl;
    for (auto &e : vl.entries()) {
      std::cout << "(" << v << " " << std::get<0>(e) << "), ";
    }
    std::cout << std::endl;
#endif

    // if (u, v) is not a tree edge, then we need to remove it from the
    // non tree adjacency list
    // TODO: could this possibly do this with edges_[se[i]] == K_TREE?
    if (ul.contains(v)) {
      non_tree_adjacency_lists_[level][u].deleteVal(v);
      non_tree_adjacency_lists_[level][v].deleteVal(u);

    } else {
      // if it is a tree edge, then we need to add it to the tree edges
      tree_edges.push_back(se[i]);

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

#ifdef DEBUG
  // print min and max level here
  std::cout << "min tree edge level: " << min_tree_edge_level << std::endl;
  std::cout << "max tree edge level: " << (int)max_level_ << std::endl;
#endif

  // delete edges from the tree at each level from the minimum tree edge level
  // to the maximum tree edge level
  for (int l = min_tree_edge_level; l < max_level_; l++) {
    auto *level_euler_tree = parallel_spanning_forests_[l];

    // get the edges to delete which have level at max l.
    auto to_delete = parlay::filter(tree_edges, [&](E e) {
      return edges_.find(pair<V, V>(e.first, e.second)).level <= l;
    });

#ifdef DEBUG
    std::cout << "Deleting edges from level " << l << std::endl;
    for (auto &e : to_delete) {
      std::cout << e.first << " " << e.second << std::endl;
    }
#endif

    // FIXME: We do an extraneously expensive copy here because of API design.
    sequence<pair<int, int>> to_delete_pair_sequence =
        EdgeBatchToPairArray(to_delete);

#ifdef DEBUG
    for (auto &e : to_delete_pair_sequence) {
      std::cout << e.first << " " << e.second << std::endl;
    }
#endif

    level_euler_tree->BatchCut(to_delete_pair_sequence);
  }

  // We now try to find replacement edges starting from the min level
  // and working our way up.
  for (int l = min_tree_edge_level; l < max_level_; l++) {
    auto *level_euler_tree = parallel_spanning_forests_[l];

    auto edges_to_replace = parlay::filter(tree_edges, [&](E e) {
      return edges_.find(pair<V, V>(e.first, e.second)).level == l;
    });

#ifdef DEBUG
    // print out the edges to replace
    std::cout << "Edges to replace" << std::endl;
    for (auto &e : edges_to_replace) {
      std::cout << "(" << e.first << ", " << e.second << "), ";
    }
    std::cout << std::endl << std::endl;
#endif

    sequence<int> lcomponents = parlay::map(edges_to_replace, [&](E e) {
      return level_euler_tree->GetRepresentative(e.first);
    });
    sequence<int> rcomponents = parlay::map(edges_to_replace, [&](E e) {
      return level_euler_tree->GetRepresentative(e.second);
    });
    lcomponents.append(rcomponents);

    auto components_to_consider = RemoveDuplicates(lcomponents);

// print out the components to consider if debug is defined
#ifdef DEBUG
    std::cout << " ----------------------" << std::endl;
    std::cout << "Components to consider:" << std::endl;

    for (auto &c : components_to_consider) {
      std::cout << "Component " << c << " contains the edges:";
      auto ett = parallel_spanning_forests_[l];
      auto cc = ett->ComponentEdges(c);
      for (auto &e : cc) {
        std::cout << "(" << e.first << ", " << e.second << "), ";
      }
      std::cout << std::endl;
    }
#endif

    // FIXME: we might have to reserve space here.
    sequence<pair<int, int>> promoted_edges;

    // Doing a replacement search for the disconnected edge at level l.
    // BUG: push down all the tree edges at level i.
    ReplacementSearch(l, components_to_consider, promoted_edges);

    if (!promoted_edges.empty()) {
      // We have some promoted edges.
      level_euler_tree->BatchLink(promoted_edges);
      // TODO: possibly coursen the loop
      // update the tree edges with the promoted edges for all the higher
      // levels
      parlay::parallel_for(l + 1, max_level_, [&](int i) {
        auto *level_euler_tree = parallel_spanning_forests_[i];
        level_euler_tree->BatchLink(promoted_edges);
      });

      // remove the promoted edges from the non-tree edge lists
      parlay::parallel_for(0, promoted_edges.size(), [&](int i) {
        auto e = promoted_edges[i];
        auto level = edges_.find(e).level;
        auto u = promoted_edges[i].first;
        auto v = promoted_edges[i].second;

// if DEBUG is defined
#ifdef DEBUG
        // print
        std::cout << "Removing edge " << e.first << " " << e.second
                  << std::endl;
        // level
        std::cout << "level: " << (int)level << std::endl;
#endif

        // update the edges_ map to reflect that it is a tree edge
        bdcty::EInfo ei = {level, bdcty::EType::K_TREE};
        edges_.insert(make_tuple(e, ei));
        edges_.insert(make_tuple(make_pair(e.second, e.first), ei));

        auto &ul = non_tree_adjacency_lists_[level][u];
        auto &vl = non_tree_adjacency_lists_[level][v];

        // if (u, v) is a tree edge, then we need to remove it from the
        // non tree adjacency list

        assert(ul.contains(v) || vl.contains(u));

        ul.deleteVal(v);
        vl.deleteVal(u);

        assert(!(ul.contains(v) || vl.contains(u)));
      });
    }
  }
}

// sequence<V> BatchDynamicConnectivity::parallelLevelSearch(
//     const sequence<E> &se, sequence<V> &components,
//     sequence<pair<int, int>> &promotedEdges, int level) {
//   auto levelEulerTree = parallel_spanning_forests_[level];
//   levelEulerTree->BatchLink(promotedEdges);

//   auto ncomponents = parlay::map(
//       components, [&](V v) { return
//       (V)levelEulerTree->GetRepresentative(v);
//       });
//   components = ncomponents;
//   // components = RemoveDuplicates(components);

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

//   sequence<E> R;

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
//         if (details.level == level && details.type ==
//         batchDynamicConnectivity::EType::K_TREE
//         &&
//             (cc.find(e.first) != cc.end() || cc.find(e.second) !=
//             cc.end()))
//             {
//           edges_[e] = batchDynamicConnectivity::EInfo{
//               static_cast<batchDynamicConnectivity::Level>(details.level -
//               1), edges_[e].type};
//           edgesToDropLevel.push_back(make_pair(e.first, e.second));
//         }
//       }
//       R.push_back(componentSearch(level, v));
//     }
//     parallel_spanning_forests_[level - 1]->BatchLink(edgesToDropLevel);

//     auto maxLevelEulerTree = parallel_spanning_forests_[max_level_];
//     auto auxiliaryEdges = parlay::map(R, [&](E e) {
//       return E(maxLevelEulerTree->GetRepresentative(e.first),
//                             maxLevelEulerTree->GetRepresentative(e.second));
//     });
//     auto promIndices = GetSpanningTree(auxiliaryEdges);

//     sequence<pair<int, int>> newPromotedEdges;
//     sequence<E> notPromotedEdges;

//     parlay::parallel_for(0, promIndices.size(), [&](int i) {
//       if (promIndices.find(auxiliaryEdges[i]) != promIndices.end())
//         newPromotedEdges.push_back(
//             make_pair(auxiliaryEdges[i].first, auxiliaryEdges[i].second));
//       else
//         notPromotedEdges.push_back(auxiliaryEdges[i]);
//     });

//     levelEulerTree->BatchLink(newPromotedEdges);

//     parlay::parallel_for(0, newPromotedEdges.size(), [&](size_t i) {
//       E e = {newPromotedEdges[i].first,
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
//       return (V)levelEulerTree->GetRepresentative(v);
//     });

//     // components_to_consider = RemoveDuplicates(ccu);

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
} // namespace bdcty

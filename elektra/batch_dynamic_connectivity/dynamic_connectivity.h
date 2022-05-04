#include <parlay/primitives.h>

#include <atomic>
#include <mutex>
#include <set>

#include "connectivity.h"
#include "utilities/sequence_utils.h"

namespace bdcty {
// -------------
// PUBLIC METHODS
// -------------
//#define DEBUG 0

// Batch Dynamic Connectivity: BatchConnected
// @input: A sequence of pairs of vertices
// @output: A sequence boolean values which are true if the corresponding
// vertices are connected
// @description: This method checks if the batch of vertices is connected in
// the graph
auto BatchDynamicConnectivity::BatchConnected(sequence<pair<V, V>> suv)
    -> sequence<char> {
  CheckRep();

  sequence<char> s_connected(suv.size(), 0);
  // check if they are connected in the highest level forest
  const auto &p_max_level_euler_tree =
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
  CheckRep();

  // Look at the max level Euler Tour Tree in the parallel spanning forests.
  const auto &max_level_euler_tree = parallel_spanning_forests_[max_level_ - 1];

  // We get the spanning tree_set of the representatives of the vertices in the
  // max_level_euler_tree. We get back a vector of type uintE.
  auto spanning_tree =
      RepresentativeSpanningTree<E, uintE>(se, max_level_euler_tree);
  auto tree_set =
      elektra::MakeSet<sequence<pair<uintE, uintE>>, E, EHash>(spanning_tree);
  // these are all the tree edges, so we insert them in to the top level.
  parlay::parallel_for(0, spanning_tree.size(), [&](size_t i) {
    pair<V, V> e = E(spanning_tree[i].first, spanning_tree[i].second);
    InsertIntoEdgeTable(e, EType::K_TREE, max_level_ - 1);
  });

  //  auto tree_edges = NewEdgeSequence<uint32_t>(spanning_tree);
  max_level_euler_tree->BatchLink(
      spanning_tree, parlay::delayed_seq<bool>(spanning_tree.size(),
                                               [](size_t) { return true; }));

  const auto non_tree_edges = parlay::filter(se, [&](auto e) {
    // TODO(sualeh): we could get rid of this check by just checking the
    //   EdgeTable which labels edges as Tree Edges. So we could get rid of that
    //   table altogether.
    return tree_set.count(e) == 0;
  });
  parlay::parallel_for(0, non_tree_edges.size(), [&](size_t i) {
    InsertIntoEdgeTable(non_tree_edges[i], EType::K_NON_TREE, max_level_ - 1);
  });

  non_tree_adjacency_lists_.BatchAddEdgesToLevel(non_tree_edges,
                                                 max_level_ - 1);
  parallel_euler_tour_tree::UpdateNontreeEdges(max_level_euler_tree.get(),
                                               true,  // is_insertion
                                               non_tree_edges);

  // add to adjacency list
  {
#ifdef DEBUG
    std::cout << "Tree edges" << std::endl;
    for (const auto &e : tree_set) {
      std::cout << "Edge: " << e.first << " " << e.second << std::endl;
    }
    PrintEdgeSequence(spanning_tree, "Tree Edges from sequence");
    std::cout << "Adding the non-tree_set edges to the hashtable" << std::endl;
    std::cout << "non_tree_edges.size() = " << non_tree_edges.size()
              << std::endl;
    PrintEdgeSequence(non_tree_edges, "nonTreeEdges");
#endif
  }
  CheckRep();
}

void BatchDynamicConnectivity::PushDownTreeEdgesFromComponents(
    Level l, sequence<V> &components) {
  auto &euler_tree = parallel_spanning_forests_[l];
  auto &lower_level_euler_tree = parallel_spanning_forests_[l - 1];

  parlay::parallel_for(0, components.size(), [&](size_t i) {
    auto c = components[i];
    auto tree_edges = euler_tree->GetAndClearLevelIEdges(c);

    // move all these edges to the euler_tree of level l-1
    lower_level_euler_tree->BatchLink(
        tree_edges, parlay::delayed_seq<bool>(tree_edges.size(),
                                              [](size_t) { return true; }));

    parlay::parallel_for(0, tree_edges.size(), [&](size_t j) {
      auto e = tree_edges[j];
      InsertIntoEdgeTable(e, EType::K_TREE, l - 1);
    });
  });
}

void BatchDynamicConnectivity::PushDownNonTreeEdges(
    Level l, const parlay::sequence<E> &non_tree_edges) {
#if DEBUG
  std::cout << "--- Pushing down non-tree edges ---" << std::endl;
  PrintEdgeSequence(non_tree_edges,
                    "Pushing down from level: " + std::to_string(l));
#endif
  auto &level_non_tree_adjacency_lists = non_tree_adjacency_lists_[l];

  // TODO(sualeh): move this function into the non_tree_adjacency_lists_ class.
  parlay::parallel_for(0, non_tree_edges.size(), [&](size_t i) {
    const auto &e = non_tree_edges[i];
    const auto &[u, v] = e;
    const auto &[kLevel, kType] = edges_.find(e);
    assert(kType == EType::K_NON_TREE && kLevel == l);
    level_non_tree_adjacency_lists[u].deleteVal(v);
    level_non_tree_adjacency_lists[v].deleteVal(u);
    assert(!level_non_tree_adjacency_lists[u].contains(v));
    assert(!level_non_tree_adjacency_lists[v].contains(u));

    // update the edges_ map to reflect the level change
    InsertIntoEdgeTable(e, EType::K_NON_TREE, l - 1);
  });

  non_tree_adjacency_lists_.BatchAddEdgesToLevel(non_tree_edges, l - 1);

  parallel_euler_tour_tree::UpdateNontreeEdges(
      parallel_spanning_forests_[l].get(),
      false,  // is_insertion
      non_tree_edges);
  parallel_euler_tour_tree::UpdateNontreeEdges(
      parallel_spanning_forests_[l - 1].get(),
      true,  // is_insertion
      non_tree_edges);
}

void BatchDynamicConnectivity::BatchDeleteEdges(
    const sequence<E> &se_unfiltered) {
  CheckRep();
  // first make sure all the edges are correctly oriented and remove all the
  // edges that are not in the graph
  const sequence<E> se = RemoveUnknownEdges(se_unfiltered, edges_);

  // split se into tree and non tree edges
  sequence<E> tree_edges;

  //  V min_tree_edge_level = static_cast<V>(max_level_);
  std::atomic<Level> min_tree_edge_level(max_level_);

  std::mutex m;
  std::unique_lock<std::mutex> lock(m, std::defer_lock);

  auto levels = parlay::sequence<std::pair<E, Level>>::uninitialized(se.size());
  // delete edges from the non tree adjacency lists.
  // TODO(sualeh): maybe collect the edges to be deleted and then delete them
  // in a batch
  parlay::parallel_for(0, se.size(), [&](V i) {
    const auto &[kU, kV] = se[i];
    const auto &[kLevel, kType] = edges_.find(se[i]);

    auto &ul = non_tree_adjacency_lists_[kLevel][kU];
    auto &vl = non_tree_adjacency_lists_[kLevel][kV];

    // if (kU, kV) is not a tree edge, then we need to remove it from the
    // non tree adjacency list. Otherwise, we have to do hard work to remove it
    // from the tree.
    if (kType == EType::K_NON_TREE) {
      levels[i] = make_pair(se[i], kLevel);
      ul.deleteVal(kV);
      vl.deleteVal(kU);
    } else {
      levels[i] = make_pair(E(kV_Max, kV_Max), kLevelMax);
      tree_edges.push_back(se[i]);
      parlay::write_min(&min_tree_edge_level, kLevel, std::less<>());
    }

    {
#ifdef DEBUG
      auto u_entries = ul.entries();
      PrintSequence(u_entries, "ul: edges for vertex kU");
      auto v_entries = vl.entries();
      PrintSequence(v_entries, "vl: edges for vertex kV");
#endif
    }
  });

  // Update non-tree edge counts in the ETTs. This requires grouping the
  // non-tree edges by level.
  levels = parlay::filter(levels,
                          [&](auto &elem) { return elem.second != kLevelMax; });
  parlay::integer_sort_inplace(levels, [&](auto &elem) {
    // we must cast because integer_sort asserts on the type being unsigned
    return static_cast<uint8_t>(elem.second);
  });
  const auto unique_level_indices = elektra::get_offsets(levels);
  parlay::parallel_for(0, unique_level_indices.size(), [&](const size_t i) {
    const size_t start = unique_level_indices[i];
    const size_t end = i == unique_level_indices.size() - 1
                           ? levels.size()
                           : unique_level_indices[i + 1];
    const Level level = levels[start].second;

    parallel_euler_tour_tree::UpdateNontreeEdges(
        parallel_spanning_forests_[level].get(),
        false,  // is_insertion
        parlay::delayed_seq<E>(end - start, [&](const size_t j) {
          return levels[start + j].first;
        }));
  });

  // delete edges from the tree at each level from the minimum tree edge level
  // to the maximum tree edge level
  for (int l = min_tree_edge_level; l < max_level_; l++) {
    const auto &level_euler_tree = parallel_spanning_forests_[l];

    // get the edges to delete which have level at max l.
    auto to_delete = parlay::filter(
        tree_edges, [&](E e) { return edges_.find(e).level <= l; });
    level_euler_tree->BatchCut(to_delete);

    {
#ifdef DEBUG
      std::cout << "Deleting edges from level " << l << std::endl;
      for (auto &e : to_delete) {
        std::cout << e.first << " " << e.second << std::endl;
      }
#endif
    }
  }

  // We now try to find replacement edges starting from the min level.
  for (int l = min_tree_edge_level; l < max_level_; l++) {
    auto &level_euler_tree = parallel_spanning_forests_[l];
    auto &lower_level_euler_tree = parallel_spanning_forests_[l - 1];

    // Is edges_to_replace correct? Should we also take into account
    // edges at levels below that weren't "replaced" below?
    auto edges_to_replace = parlay::filter(tree_edges, [&](E e) {
      return edges_.find(pair<V, V>(e.first, e.second)).level == l;
    });
    auto components = sequence<V>::uninitialized(2 * edges_to_replace.size());
    parlay::parallel_for(0, edges_to_replace.size(), [&](size_t i) {
      components[2 * i] =
          level_euler_tree->GetRepresentative(edges_to_replace[i].first);
      components[2 * i + 1] =
          level_euler_tree->GetRepresentative(edges_to_replace[i].second);
    });

    auto components_to_consider = RemoveDuplicates(std::move(components));
    {
#ifdef DEBUG
      std::cout << "----------------- Key loop -----------------" << std::endl;
      cout << "Level " << l << endl;
      PrintEdgeSequence(edges_to_replace, "Edges to replace");
      std::cout << " ------Components to consider:-------" << std::endl;
      for (auto &c : components_to_consider) {
        auto cc = parallel_spanning_forests_[l]->ComponentEdges(c);
        PrintEdgeSequence(cc, "Component " + std::to_string(c) + " edges");
      }
#endif
    }

    // push down all the small component edges to a lower level.
    const auto kCriticalComponentSize = 1U << (l - 1);
    auto small_components = parlay::filter(components_to_consider, [&](V c) {
      return parallel_spanning_forests_[l]->ComponentSize(c) <=
             kCriticalComponentSize;
    });
    PushDownTreeEdgesFromComponents(l, small_components);

    // Doing a replacement search for the disconnected edge at level l.
    // TODO: use ::uninitialized and don't use push_back.
    const auto replacement_search_output =
        ReplacementSearch(l, components_to_consider);
    const auto &promoted_edges = replacement_search_output.first;
    const auto &promoted_edges_to_push = replacement_search_output.second;
    const auto promoted_edges_to_push_seq = promoted_edges_to_push.keys();

    {
#ifdef DEBUG
      PrintEdgeSequence(promoted_edges, "Promoted Edges");
      PrintStructure();
#endif
    }

    if (!promoted_edges.empty()) {
      level_euler_tree->BatchLink(
          promoted_edges, parlay::map(promoted_edges, [&](const E e) -> bool {
            return !promoted_edges_to_push.contains(e);
          }));
      lower_level_euler_tree->BatchLink(
          promoted_edges_to_push_seq,
          parlay::delayed_seq<bool>(promoted_edges_to_push_seq.size(),
                                    [](size_t) { return true; }));
      parallel_euler_tour_tree::UpdateNontreeEdges(level_euler_tree.get(),
                                                   false,  // is_insertion
                                                   promoted_edges);

      // TODO: possibly coursen the loop
      // update all the higher levels with the promoted tree edges
      parlay::parallel_for(
          l + 1, max_level_,
          [&](int i) {
            auto &level_euler_tree = parallel_spanning_forests_[i];
            level_euler_tree->BatchLink(
                promoted_edges,
                parlay::delayed_seq<bool>(promoted_edges.size(),
                                          [](size_t) { return false; }));
          },
          1);
      // remove the promoted edges from the non-tree edge lists
      parlay::parallel_for(0, promoted_edges.size(), [&](auto i) {
        const auto kE = promoted_edges[i];
        const auto &[kU, kV] = kE;
        const Level new_level = promoted_edges_to_push.contains(kE) ? l - 1 : l;
        // update the edges_ map to reflect that it is a tree edge
        InsertIntoEdgeTable(kE, EType::K_TREE, new_level);
        // then we remove it from the non tree adjacency list.
        auto &ul = non_tree_adjacency_lists_[l][kU];
        auto &vl = non_tree_adjacency_lists_[l][kV];
        assert(ul.contains(kV) || vl.contains(kU));
        ul.deleteVal(kV);
        vl.deleteVal(kU);
        assert(!(ul.contains(kV) || vl.contains(kU)));
      });
    }
  }
  parlay::parallel_for(0, se.size(),
                       [&](const size_t i) { DeleteFromEdgeTable(se[i]); });

  CheckRep();
}

std::pair<sequence<E>,
          elektra::resizable_table<E, elektra::empty, HashIntPairStruct>>
BatchDynamicConnectivity::ReplacementSearch(Level level,
                                            sequence<V> components) {
  // pseudocode:
  //
  // M: map from component to super-component size
  // search_size = 256
  // edges_to_promote = []
  //
  // // edges to be pushed
  // tree_edges_to_push = []
  // nontree_edges_to_push = []
  //
  // while components > 1:
  //   filter out components s.t. M[...] > 2^(i-1)
  //
  //   for each component C:
  //     parent[C] = uf.find(C)
  //
  //   for each component C:
  //     for the first `search_size`-th level-i non-tree edges (u, v) of C:
  //       if uf.union(u, v):  // tree edge
  //         edges_to_promote.insert(u, v)
  //
  //   for each component C:
  //     current_parent = uf.find(C)
  //     if parent[C] != current_parent:
  //       M[current_parent] += M[parent[C]]
  //
  //   for each component C:
  //     if M[uf.find(C)] > 2^(i-1): // component grew too large this round
  //       continue
  //     for the first `search_size`-th level-i non-tree edges (u, v) of C:
  //       if (u, v) in edges_to_promote:
  //         tree_edges_to_push.insert(u, v)
  //       else:
  //         non_tree_edges_to_push.insert(u, v)
  //
  //   search_size *= 2
  //
  // update ETTs and other stuff:
  //   etts[level - 1].link(tree_edges_to_push)
  //   push non_tree_edges_to_push to level - 1 as non-tree edges
  //   etts[level].link(edges_to_promote)

  // set up a coarse search size and an overall stride for the search.
  constexpr auto kInitialSearchSize = 256;
  uintE search_size = kInitialSearchSize;
  uintE search_offset = 0;
  const auto kCriticalComponentSize = 1U << (level - 1);

  const auto &ett = parallel_spanning_forests_[level];
  auto &level_non_tree_adjacency_lists = non_tree_adjacency_lists_[level];

  // successful replacement edges that should be added
  auto edges_to_promote =
      elektra::resizable_table<E, elektra::empty, HashIntPairStruct>(
          2 * components.size(), edge_with_empty_struct_,
          tombstone_with_empty_struct_, HashIntPairStruct());
  // successful replacement edges (appearing in edges_to_promote) that should be
  // pushed to the next level
  auto tree_edges_to_push =
      elektra::resizable_table<E, elektra::empty, HashIntPairStruct>(
          2 * components.size(), edge_with_empty_struct_,
          tombstone_with_empty_struct_, HashIntPairStruct());
  // unsuccessful replacement edges that should be pushed to the next level
  auto nontree_edges_to_push =
      elektra::resizable_table<E, elektra::empty, HashIntPairStruct>(
          2 * components.size(), edge_with_empty_struct_,
          tombstone_with_empty_struct_, HashIntPairStruct());

  auto union_find = elektra::UnionFindCompress{
      components.size()};  // TODO(laxman): read this :-)
  // an indicator of whether we are done with that component.
  sequence<char> component_indicator(components.size(), 0);
  // a component -> idx map
  auto component_map = elektra::MakeIndexMap<V, V>(components);
  // A supercomponent is a root inside `union_find`, and
  // supercomponent_sizes[i] is the size of that supercomponent.
  auto supercomponent_sizes = sequence<std::atomic<V>>::from_function(
      components.size(),
      [&](size_t i) { return ett->ComponentSize(components[i]); });
  // component_parents[i] == cached value of union_find.find(i)
  auto component_parents = sequence<V>::uninitialized(components.size());
  auto edge_finders =
      sequence<parallel_euler_tour_tree::NontreeEdgeFinder>::from_function(
          components.size(), [&](size_t i) {
            return ett->CreateNontreeEdgeFinder(components[i]);
          });

  while (parlay::count(component_indicator, 1) < components.size()) {
    // save supercomponent IDs of the start of this while-loop iteration
    parlay::parallel_for(0, components.size(), [&](const size_t i) {
      component_parents[i] = union_find.find(i, union_find.parents);
    });

    // do search over the components
    parlay::parallel_for(0, components.size(), [&](const size_t i) {
      if (component_indicator[i] == 1) {
        return;
      }

      edge_finders[i].ForEachIncidentVertex(
          search_offset, search_offset + search_size,
          [&](V u, V adj_begin, V adj_end) {
            // TODO(sualeh) TODO(tom): entries() here is inefficient. What we
            // "should" do is only look at the `adj_begin`-th to `adj_end`-th
            // edges in level_non_tree_adjacency_lists[u]. But we don't have a
            // nice way of fetching such a range of edges out of
            // level_non_tree_adjacency_lists[u]. Possible fixes:
            //   - Implement a doubling/block search on the hash table
            //   - Maintain an implicit augmented tree (with high fanout so the
            //     overhead is small) over the hash table that counts elements
            //     over ranges of the hash table
            //   - Make level_non_tree_adjacency_lists into a parallel BST
            //     - seems annoying to port a parallelBST, and would be harder
            //       to use since it would no longer be concurrent
            //
            // The hacky thing we're doing here is just searching over all the
            // incident edges of a vertex the first time that edge_finder finds
            // it (i.e., when adj_begin == 0). The bad case that could still
            // happen here is if we hit a high degree vertex v early on in the
            // search, and .entries() costs Omega(deg(v)).
            if (adj_begin != 0) {
              return;
            }

            const auto incident_edges =
                level_non_tree_adjacency_lists[u].entries();

            const auto u_component =
                component_map.at(ett->GetRepresentative(u));
            parlay::parallel_for(0, incident_edges.size(), [&](const size_t j) {
              const V v = std::get<0>(incident_edges[j]);

              {
#ifdef DEBUG
                std::cout << "working on edge: " << u << " " << v << std::endl;
#endif
              }

              const auto v_component =
                  component_map.at(ett->GetRepresentative(v));
              if (union_find.unite(u_component, v_component,
                                   union_find.parents) != UINT_E_MAX) {
                // found a valid replacement edge
                const auto edge = make_pair(std::min(u, v), std::max(u, v));
                edges_to_promote.insert(make_tuple(edge, elektra::empty{}));
              }
            });
          });
    });

    // note(sualeh/tom): are the atomics slow here? if so we can semisort by
    // supercomponent and reduce over each supercomponent. or we can implement
    // the non-interleaved replacement search algorithm
    parlay::parallel_for(0, components.size(), [&](const size_t i) {
      const V supercomponent = union_find.find(i, union_find.parents);
      const V old_component = component_parents[i];
      if (supercomponent != old_component) {
        supercomponent_sizes[supercomponent] +=
            supercomponent_sizes[old_component];
      }
    });

    // Determine which edges to push down to the next level (amortizing the
    // cost of searching)
    parlay::parallel_for(0, components.size(), [&](const size_t i) {
      if (component_indicator[i] == 1) {
        return;
      }
      const V supercomponent = union_find.find(i, union_find.parents);
      if (supercomponent_sizes[supercomponent] > kCriticalComponentSize) {
        // The supercomponent grew too large this round. (We should not push
        // edges in this case.)
        component_indicator[i] = 1;
      } else {
        const size_t search_end = search_offset + search_size;
        edge_finders[i].ForEachIncidentVertex(
            search_offset, search_end, [&](V u, V adj_begin, V adj_end) {
              if (adj_begin != 0) {
                return;
              }
              const auto incident_edges =
                  level_non_tree_adjacency_lists[u].entries();
              parlay::parallel_for(0, incident_edges.size(), [&](size_t j) {
                const V v = std::get<0>(incident_edges[j]);
                const auto edge = make_pair(std::min(u, v), std::max(u, v));
                if (edges_to_promote.contains(edge)) {
                  tree_edges_to_push.insert(make_tuple(edge, elektra::empty{}));
                } else {
                  nontree_edges_to_push.insert(
                      make_tuple(edge, elektra::empty{}));
                }
              });
            });
        if (edge_finders[i].NumEdges() <=
            search_end) {  // no more edges to search
          component_indicator[i] = 1;
        }
      }
    });
    search_offset += search_size;
    search_size *= 2;
  }  // end while loop

  PushDownNonTreeEdges(level, nontree_edges_to_push.keys());
  return {edges_to_promote.keys(), tree_edges_to_push};
}

}  // namespace bdcty

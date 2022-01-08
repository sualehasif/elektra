#include <parlay/primitives.h>

#include <atomic>
#include <mutex>
#include <set>

#include "connectivity.h"
namespace bdcty {
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
  const auto &p_max_level_euler_tree =
      parallel_spanning_forests_[max_level_ - 1];

  parlay::parallel_for(0, suv.size(), [&](size_t i) {
    auto[v_1, v_2] = suv[i];
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
  const auto &max_level_euler_tree = parallel_spanning_forests_[max_level_ - 1];

  // We get the spanning tree_set of the representatives of the vertices in the
  // max_level_euler_tree. We get back a vector of type uintE.
  auto spanning_tree =
      RepresentativeSpanningTree<E, uintE>(se, max_level_euler_tree);
  auto tree_set =
      elektra::MakeSet<sequence<pair<uintE, uintE>>, E, EHash>(spanning_tree);
  parlay::parallel_for(0, spanning_tree.size(), [&](size_t i) {
    pair<V, V> e = E(spanning_tree[i].first, spanning_tree[i].second);
    InsertIntoEdgeTable(e, EType::K_TREE, max_level_ - 1);
  });
  // TODO(major: sualeh): get rid of this and just use the spanning tree_set
  // vector.
  auto tree_edges = NewEdgeSequence<int>(spanning_tree);
  max_level_euler_tree->BatchLink(tree_edges);

  auto non_tree_edges = sequence<E>(se.size(), E(-1, -1));
  // update the nonTree edges based on the ST computation
  parlay::parallel_for(0, se.size(), [&](size_t i) {
    if (tree_set.count(se[i]) == 0) {
      non_tree_edges[i] = (make_pair(se[i].first, se[i].second));
      InsertIntoEdgeTable(se[i], EType::K_NON_TREE, max_level_ - 1);
    }
  });
  // filter and add non-tree_set edges
  auto filtered_non_tree_edges = parlay::filter(
      non_tree_edges, [&](auto e) { return e.first != -1 && e.second != -1; });

  non_tree_adjacency_lists_.BatchAddEdgesToLevel(filtered_non_tree_edges,
                                                 max_level_ - 1);
  // add to adjacency list
#ifdef DEBUG
  std::cout << "Tree edges" << std::endl;
  for (auto &e : tree_set) {
    std::cout << "Edge: " << e.first << " " << e.second << std::endl;
  }
  PrintEdgeSequence(tree_edges, "Tree Edges from sequence");
  std::cout << "Adding the non-tree_set edges to the hashtable" << std::endl;
  std::cout << "non_tree_edges.size() = " << non_tree_edges.size() << std::endl;
  PrintEdgeSequence(non_tree_edges, "nonTreeEdgesUnfiltered");
  PrintEdgeSequence(filtered_non_tree_edges, "nonTreeEdgesFiltered");
#endif
}

void BatchDynamicConnectivity::BatchDeleteEdges(sequence<E> &se) {
  // first make sure all the edges are correctly oriented and remove all the
  // edges that are not in the graph
  RemoveUnknownEdges(se, edges_, kEmptyInfo);

  // split se into tree and non tree edges
  sequence<E> tree_edges;
  tree_edges.reserve(se.size());

  //  V min_tree_edge_level = static_cast<V>(max_level_);
  std::atomic<V> min_tree_edge_level(max_level_);

  std::mutex m;
  std::unique_lock<std::mutex> lock(m, std::defer_lock);

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
      ul.deleteVal(kV);
      vl.deleteVal(kU);
    } else {
      tree_edges.push_back(se[i]);
      parlay::write_min(&min_tree_edge_level, kLevel, std::less<>());
    }
#ifdef DEBUG
    auto u_entries = ul.entries();
    PrintSequence(u_entries, "ul: edges for vertex kU");
    auto v_entries = vl.entries();
    PrintSequence(v_entries, "vl: edges for vertex kV");
#endif
  });
  // delete edges from the tree at each level from the minimum tree edge level
  // to the maximum tree edge level
  for (int l = min_tree_edge_level; l < max_level_; l++) {
    const auto &level_euler_tree = parallel_spanning_forests_[l];

    // get the edges to delete which have level at max l.
    auto to_delete = parlay::filter(
        tree_edges, [&](E e) { return edges_.find(e).level <= l; });
    level_euler_tree->BatchCut(to_delete);

#ifdef DEBUG
    std::cout << "Deleting edges from level " << l << std::endl;
    for (auto &e : to_delete) {
      std::cout << e.first << " " << e.second << std::endl;
    }
#endif
  }

  // We now try to find replacement edges starting from the min level.
  for (int l = min_tree_edge_level; l < max_level_; l++) {
    const auto &level_euler_tree = parallel_spanning_forests_[l];

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
#ifdef DEBUG
    PrintEdgeSequence(edges_to_replace, "Edges to replace");
    std::cout << " ----------Components to consider:------------" << std::endl;
    for (auto &c : components_to_consider) {
      auto cc = parallel_spanning_forests_[l]->ComponentEdges(c);
      PrintEdgeSequence(cc, "Component" + std::to_string(c) + " edges");
    }
#endif
    // TODO(sualeh: BUG): push down all the tree edges at level i.

    // Doing a replacement search for the disconnected edge at level l.
    sequence<E> promoted_edges;
    promoted_edges.reserve(tree_edges.size());
    ReplacementSearch(l, components_to_consider, promoted_edges);

    if (!promoted_edges.empty()) {
      level_euler_tree->BatchLink(promoted_edges);
      // TODO: possibly coursen the loop
      // update all the higher levels with the promoted tree edges
      parlay::parallel_for(l + 1, max_level_, [&](int i) {
        auto &level_euler_tree = parallel_spanning_forests_[i];
        level_euler_tree->BatchLink(promoted_edges);
      });
      // remove the promoted edges from the non-tree edge lists
      parlay::parallel_for(0, promoted_edges.size(), [&](auto i) {
        const auto kE = promoted_edges[i];
        const auto &[kU, kV] = kE;
        const auto &[kLevel, kType] = edges_.find(kE);
        // update the edges_ map to reflect that it is a tree edge
        InsertIntoEdgeTable(kE, EType::K_TREE, kLevel);
        // then we remove it from the non tree adjacency list.
        auto &ul = non_tree_adjacency_lists_[kLevel][kU];
        auto &vl = non_tree_adjacency_lists_[kLevel][kV];
        assert(ul.contains(kV) || vl.contains(kU));
        ul.deleteVal(kV);
        vl.deleteVal(kU);
        assert(!(ul.contains(kV) || vl.contains(kU)));
#ifdef DEBUG
        std::cout << "Removing edge " << kU << " " << kV
                  << " at level: " << kLevel << std::endl;
#endif
      });
    }
  }
}

void BatchDynamicConnectivity::ReplacementSearch(Level level,
                                                 sequence<V> components,
                                                 sequence<E> &promoted_edges) {
  // set up a course search size and an overall stride for the search.
  constexpr auto kInitialSearchSize = 256;
  uintE k_search_size = kInitialSearchSize;
  uintE total_search_stride = 0;
  const auto kCriticalComponentSize = 1U << (level - 1);

  const auto &ett = parallel_spanning_forests_[level];
  sequence<pair<V, V>> push_down_edges;
  auto union_find = elektra::UnionFindCompress{components.size()};

  // an indicator of whether we are done with that component.
  std::vector<char> component_indicator(components.size(), 0);
  // a component -> idx map
  auto component_map = elektra::MakeIndexMap<V, V>(components);

  // set up non-tree edges to search through
  // The first sequence indexes by the component number
  // The second sequence indexes by the edges in the component
  auto search_edges = GetSearchEdges(level, components, ett);

  while (parlay::count(component_indicator, 1) < components.size()) {
    parlay::parallel_for(0, components.size(), [&](int i) {
      // early return if the component is already marked
      if (component_indicator[i] == 1) {
        return;
      }

      auto non_tree_edges_i = search_edges[i];
      const uintE kNumEdges = non_tree_edges_i.size();
      // early return if the component is big
      if (kNumEdges > kCriticalComponentSize) {
        return;
      }

      // doing a search over the components.
      parlay::parallel_for(0, min(kNumEdges, k_search_size) + 1, [&](int j) {
        auto idx = j + total_search_stride;
        // if we have already searched this component, skip it.
        if (idx >= kNumEdges) {
          component_indicator[i] = 1;
          return;
        }

        const auto &[kU, kV] = non_tree_edges_i[idx];
        auto u_component_rep = component_map.at(ett->GetRepresentative(kU));
        auto v_component_rep = component_map.at(ett->GetRepresentative(kV));
        auto u_parent = union_find.find(u_component_rep, union_find.parents);
        auto v_parent = union_find.find(v_component_rep, union_find.parents);

        if (u_parent != v_parent) {
          // link the two components together: we found an edge that connects
          // them
          union_find.unite(u_parent, v_parent, union_find.parents);
          promoted_edges.push_back(make_pair(kU, kV));
        } else {
          // this failed and will be pushed down a level.
          push_down_edges.push_back(make_pair(kU, kV));
        }
      });

      // TODO: push down the edges and clear out the push_down_edges
      // Note we will promote the edges later.
    });

    // update the total search stride and search size
    total_search_stride += k_search_size;
    k_search_size *= 2;
  }

#ifdef DEBUG
  PrintEdgeSequence(promoted_edges, "Promoted Edges");
#endif
}

}  // namespace bdcty

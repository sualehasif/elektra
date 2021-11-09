#include "elektra/hash_pair.h"
#include "elektra/resizable_table.h"

// write a main function

typedef int V;

typedef int8_t Level;

enum class EdgeType {
  // Edge is in the spanning forest of the graph.
  kNonTree,
  // Edge is not in the spanning forest of the graph.
  kTree,
};

// A struct that contains information about a particular edge.
struct EdgeInfo {
  Level level;
  EdgeType type;

  //   EdgeInfo() : level(-1), type(EdgeType::kNonTree) {}

  // equality operator
  bool operator==(const EdgeInfo &other) const {
    return level == other.level && type == other.type;
  }
};

using std::make_pair;
using std::pair;

int main() {
  elektra::resizable_table<pair<V, V>, EdgeInfo, HashIntPairStruct> table(
      15,
      std::make_tuple(std::make_pair(-1, -1), EdgeInfo{-1, EdgeType::kNonTree}),
      HashIntPairStruct());

  table.insert(std::make_tuple(make_pair(1, 1), EdgeInfo{1, EdgeType::kTree}));
  table.insert(std::make_tuple(make_pair(1, 2), EdgeInfo{4, EdgeType::kTree}));
  table.insert(
      std::make_tuple(make_pair(2, 4), EdgeInfo{2, EdgeType::kNonTree}));

  if (table.find(std::make_pair(1, 2)).type == EdgeType::kTree) {
    std::cout << "Edge (1, 2) is in the spanning forest of the graph."
              << std::endl;
  } else {
    std::cout << "Edge (1, 2) is not in the spanning forest of the graph."
              << std::endl;
  }

  if (table.find(std::make_pair(2, 4)).type == EdgeType::kTree) {
    std::cout << "Edge (1, 2) is in the spanning forest of the graph."
              << std::endl;
  } else {
    std::cout << "Edge (1, 2) is not in the spanning forest of the graph."
              << std::endl;
  }

  if (table.find(std::make_pair(2, 4)).level == 2) {
    std::cout << "Edge (1, 2) has weight 2" << std::endl;
  } else {
    std::cout << "Edge (1, 2) has incorrect weight" << std::endl;
  }

  return 0;
}
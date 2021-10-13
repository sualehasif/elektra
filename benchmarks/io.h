#pragma once
#include <parlay/sequence.h>

#include <unordered_set>

#include "elektra/macros.h"

template <typename T>
using sequence = parlay::sequence<T>;
using std::pair;

namespace elektra {
namespace io {

typedef std::pair<uintE, uintE> intPair;

namespace internal {

// Starting from the current position, skips all consecutive lines in the stream
// that start with '#' or are empty.
//
// The intent here is that lines starting with '#' are interpreted to be
// comments that should be ignored.
void skip_ifstream_comments(std::ifstream* stream) {
  std::string line;
  while (*stream) {
    std::streampos current_position = stream->tellg();
    std::getline(*stream, line);
    if (!(line.empty() || line[0] == '#')) {
      stream->seekg(current_position);
      return;
    }
  }
}

}  // namespace internal

// convert the graph into a list of pairs of edges in a sequence
// return as a pair of the number of vertices and the edge sequence
pair<intV, sequence<intPair>> read_unweighted_edge_list(const char* filename) {
  std::ifstream file{filename};
  if (!file.is_open()) {
    std::cout << "ERROR: Unable to open file: " << filename << '\n';
    std::terminate();
  }
  internal::skip_ifstream_comments(&file);

  intV num_vertices = 0;
  std::unordered_set<intV> vertices;

  sequence<intPair> edge_list;
  uintE from;
  uintE to;
  while (file >> from >> to) {
    if (vertices.find(from) == vertices.end()) {
      vertices.insert(from);
      num_vertices++;
    }
    if (vertices.find(to) == vertices.end()) {
      vertices.insert(to);
      num_vertices++;
    }
    // edge_list.push_back(std::make_pair(from, to));
    edge_list.emplace_back(from, to);
  }
  return std::make_pair(num_vertices, edge_list);
}

}  // namespace io
}  // namespace elektra
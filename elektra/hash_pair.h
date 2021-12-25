#pragma once

#include <cstdlib>
#include <utility>

// TODO(sualeh, urgent): very cleanly change this to use absl::hash_combine.
inline unsigned int hashInt(unsigned int a) {
  a = (a + 0x7ed55d16) + (a << 12);
  a = (a ^ 0xc761c23c) ^ (a >> 19);
  a = (a + 0x165667b1) + (a << 5);
  a = (a + 0xd3a2646c) ^ (a << 9);
  a = (a + 0xfd7046c5) + (a << 3);
  a = (a ^ 0xb55a4f09) ^ (a >> 16);
  return a;
}

inline auto HashIntPair(const std::pair<unsigned int, unsigned int> &p) -> unsigned {
  unsigned h{hashInt(p.first)};
  h ^= hashInt(p.second) + 0x9e3779b9 + (h << 6) + (h >> 2);
  return h;
}

// For use in hash containers. For instance:
//   std::unordered_map<std::pair<int, int>, std::string, HashIntPairStruct>
//     int_to_string_map;
struct HashIntPairStruct {
  auto operator()(const std::pair<int, int> &p) const -> size_t {
    return HashIntPair(p);
  }
  auto operator()(const std::pair<unsigned int, unsigned int> &p) const -> size_t {
    return HashIntPair(p);
  }
};

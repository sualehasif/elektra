//
// Created by Sualeh Asif on 12/24/21.
//
#pragma once

#include "hash_pair.h"
#include "resizable_table.h"
#include <unordered_map>
#include <unordered_set>

namespace elektra {
template<class T, class K, class KeyHash>
inline auto MakeSet(T &seq) -> std::unordered_set<K, KeyHash> {
  auto set = std::unordered_set<K, KeyHash>{};
  for (auto &item : seq) {
    set.insert(item);
  }
  return set;
}

template<class K, class V=size_t, class T>
inline auto MakeIndexMap(T &seq) -> std::unordered_map<K, V> {
  auto map = std::unordered_map<K, V>{};
  for (V i = 0; i < static_cast<V>(seq.size()); i++) {
    map.insert(std::make_pair(seq[i], i));
  }
  return map;
}

template<class K, class V=size_t, class KeyHash, class T>
inline auto MakeIndexMap(T &seq) -> std::unordered_map<K, V, KeyHash> {
  auto map = std::unordered_map<K, V, KeyHash>{};
  for (size_t i = 0; i < seq.size(); i++) {
    map.insert(std::make_pair(seq[i], i));
  }
  return map;
}
} // namespace elektra
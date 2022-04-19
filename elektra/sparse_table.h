#pragma once

#include <parlay/primitives.h>
#include <parlay/sequence.h>
#include <parlay/slice.h>

#include <tuple>
#include <unordered_set>
#include <vector>

#include "macros.h"

namespace elektra {

struct HashUint32Empty {
  static inline size_t hash(uint32_t key) { return parlay::hash32_3(key); }
  static constexpr uint32_t empty_key = std::numeric_limits<uint32_t>::max();
  static constexpr std::pair<uint32_t, elektra::empty> empty = {
      empty_key, elektra::empty()};
  static constexpr uint32_t tombstone =
      std::numeric_limits<uint32_t>::max() - 1;
  static constexpr double space_mult = 1.1;
};

// Overfilling the table could put it into an infinite loop.
template <class K, class V, class KVInfo>
class sparse_table {
 public:
  using T = std::tuple<K, V>;

  sequence<T> table;

  size_t capacity() const { return table.size(); }
  size_t size() const {
    return parlay::reduce(
        parlay::delayed_seq<size_t>(capacity(), [&](size_t i) {
          return std::get<0>(table[i]) != KVInfo::empty_key &&
                 std::get<0>(table[i]) != KVInfo::tombstone;
        }));
  }

  inline size_t hashToRange(size_t h) const { return h & (table.size() - 1); }
  inline size_t firstIndex(K& k) const { return hashToRange(KVInfo::hash(k)); }
  inline size_t incrementIndex(size_t h) const { return hashToRange(h + 1); }

  sparse_table() {}

  sparse_table(size_t _m) {
    size_t m =
        (size_t)1 << parlay::log2_up((size_t)(KVInfo::space_mult * _m) + 1);
    table = sequence<T>::uninitialized(m);
    clear_table();
  }

  bool insert(std::tuple<K, V> kv) {
    K k = std::get<0>(kv);
    size_t h = firstIndex(k);
    while (true) {
      if (std::get<0>(table[h]) == KVInfo::empty_key) {
        if (elektra::atomic_compare_and_swap(&std::get<0>(table[h]),
                                             KVInfo::empty_key, k)) {
          if
            constexpr(sizeof(V) > 0) {
              std::get<1>(table[h]) = std::get<1>(kv);
            }
          return true;
        }
      } else if (std::get<0>(table[h]) == KVInfo::tombstone) {
        if (elektra::atomic_compare_and_swap(&std::get<0>(table[h]),
                                             KVInfo::tombstone, k)) {
          if
            constexpr(sizeof(V) > 0) {
              std::get<1>(table[h]) = std::get<1>(kv);
            }
          return true;
        }
      } else if (std::get<0>(table[h]) == k) {
        return false;
      }
      h = incrementIndex(h);
    }
    return false;
  }

  bool deleteVal(K k) {
    size_t h = firstIndex(k);
    while (true) {
      if (std::get<0>(table[h]) == KVInfo::empty_key) {
        return false;
      } else if (std::get<0>(table[h]) == k) {
        return elektra::atomic_compare_and_swap(&std::get<0>(table[h]), k,
                                                KVInfo::tombstone);
      }
      h = incrementIndex(h);
    }
  }

  bool contains(K k) const {
    size_t h = firstIndex(k);
    while (true) {
      if (std::get<0>(table[h]) == k) {
        return true;
      } else if (std::get<0>(table[h]) == KVInfo::empty_key) {
        return false;
      }
      h = incrementIndex(h);
    }
    return false;
  }

  V find(K k, V default_value) const {
    size_t h = firstIndex(k);
    while (true) {
      if (std::get<0>(table[h]) == k) {
        return std::get<1>(table[h]);
      } else if (std::get<0>(table[h]) == KVInfo::empty_key) {
        return default_value;
      }
      h = incrementIndex(h);
    }
    return default_value;
  }

  sequence<T> entries() const {
    auto pred = [&](const T& t) {
      return std::get<0>(t) != KVInfo::empty_key &&
             std::get<0>(t) != KVInfo::tombstone;
    };
    return parlay::filter(table, pred);
  }

  // Not sure if there are any users of this function, currently.
  parlay::slice<T, T> get_slice(size_t start, size_t length) {
    return parlay::make_slice(
        table.begin() + start,
        table.begin() + std::min(start + length, capacity()));
  }

  // Incoming must be a power of two
  void maybe_resize(size_t incoming) {
    // TODO: use random sampling for large arrays
    size_t num_full = size();
    if ((num_full + incoming) * KVInfo::space_mult > capacity()) {
      size_t new_size = 1 << parlay::log2_up((size_t)(KVInfo::space_mult *
                                                      (num_full + incoming)));

      auto old_backing = std::move(table);
      table = sequence<T>(new_size, KVInfo::empty);

      parlay::parallel_for(0, old_backing.size(), [&](size_t i) {
        if (std::get<0>(old_backing[i]) != KVInfo::empty_key &&
            std::get<0>(old_backing[i]) != KVInfo::tombstone) {
          insert(old_backing[i]);
        }
      });
    }
  }

  void clear_table() {
    parlay::parallel_for(0, table.size(),
                         [&](size_t i) { table[i] = KVInfo::empty; });
  }
};

}  // namespace elektra
